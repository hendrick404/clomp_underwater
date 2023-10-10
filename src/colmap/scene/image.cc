// Copyright (c) 2023, ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include "colmap/scene/image.h"

#include "colmap/geometry/pose.h"
#include "colmap/scene/projection.h"

namespace colmap {
namespace {

static constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

}  // namespace

const int Image::kNumPoint3DVisibilityPyramidLevels = 6;

Image::Image()
    : image_id_(kInvalidImageId),
      name_(""),
      camera_id_(kInvalidCameraId),
      registered_(false),
      num_points3D_(0),
      num_observations_(0),
      num_correspondences_(0),
      num_visible_points3D_(0),
      cam_from_world_prior_(Eigen::Quaterniond(kNaN, kNaN, kNaN, kNaN),
                            Eigen::Vector3d(kNaN, kNaN, kNaN)),
      cam_from_world_prior_cov_(Eigen::Matrix7d::Zero()) {}

void Image::SetUp(const class Camera& camera) {
  CHECK_EQ(camera_id_, camera.CameraId());
  point3D_visibility_pyramid_ = VisibilityPyramid(
      kNumPoint3DVisibilityPyramidLevels, camera.Width(), camera.Height());
}

void Image::TearDown() {
  point3D_visibility_pyramid_ = VisibilityPyramid(0, 0, 0);
}

void Image::SetPoints2D(const std::vector<Eigen::Vector2d>& points) {
  CHECK(points2D_.empty());
  points2D_.resize(points.size());
  num_correspondences_have_point3D_.resize(points.size(), 0);
  for (point2D_t point2D_idx = 0; point2D_idx < points.size(); ++point2D_idx) {
    points2D_[point2D_idx].xy = points[point2D_idx];
  }
}

void Image::SetPoints2D(const std::vector<struct Point2D>& points) {
  CHECK(points2D_.empty());
  points2D_ = points;
  num_correspondences_have_point3D_.resize(points.size(), 0);
  num_points3D_ = 0;
  for (const auto& point2D : points2D_) {
    if (point2D.HasPoint3D()) {
      num_points3D_ += 1;
    }
  }
}

void Image::SetPoint3DForPoint2D(const point2D_t point2D_idx,
                                 const point3D_t point3D_id) {
  CHECK_NE(point3D_id, kInvalidPoint3DId);
  struct Point2D& point2D = points2D_.at(point2D_idx);
  if (!point2D.HasPoint3D()) {
    num_points3D_ += 1;
  }
  point2D.point3D_id = point3D_id;
}

void Image::ResetPoint3DForPoint2D(const point2D_t point2D_idx) {
  struct Point2D& point2D = points2D_.at(point2D_idx);
  if (point2D.HasPoint3D()) {
    point2D.point3D_id = kInvalidPoint3DId;
    num_points3D_ -= 1;
  }
}

bool Image::HasPoint3D(const point3D_t point3D_id) const {
  return std::find_if(points2D_.begin(),
                      points2D_.end(),
                      [point3D_id](const struct Point2D& point2D) {
                        return point2D.point3D_id == point3D_id;
                      }) != points2D_.end();
}

void Image::IncrementCorrespondenceHasPoint3D(const point2D_t point2D_idx) {
  const struct Point2D& point2D = points2D_.at(point2D_idx);

  num_correspondences_have_point3D_[point2D_idx] += 1;
  if (num_correspondences_have_point3D_[point2D_idx] == 1) {
    num_visible_points3D_ += 1;
  }

  point3D_visibility_pyramid_.SetPoint(point2D.xy(0), point2D.xy(1));

  assert(num_visible_points3D_ <= num_observations_);
}

void Image::DecrementCorrespondenceHasPoint3D(const point2D_t point2D_idx) {
  const struct Point2D& point2D = points2D_.at(point2D_idx);

  num_correspondences_have_point3D_[point2D_idx] -= 1;
  if (num_correspondences_have_point3D_[point2D_idx] == 0) {
    num_visible_points3D_ -= 1;
  }

  point3D_visibility_pyramid_.ResetPoint(point2D.xy(0), point2D.xy(1));

  assert(num_visible_points3D_ <= num_observations_);
}

Eigen::Vector3d Image::ProjectionCenter() const {
  return cam_from_world_.rotation.inverse() * -cam_from_world_.translation;
}

Eigen::Vector3d Image::ViewingDirection() const {
  return cam_from_world_.rotation.toRotationMatrix().row(2);
}

void Image::SetVirtualPoints2D(const std::vector<Eigen::Vector2d>& points) {
  CHECK(virtual_points2D_.empty());
  virtual_points2D_.resize(points.size());
  for (point2D_t point2D_idx = 0; point2D_idx < points.size(); ++point2D_idx) {
    points2D_[point2D_idx].xy = points[point2D_idx];
  }
}

void Image::ComputeVirtualTransformations(
    const Camera& camera,
    const std::vector<Eigen::Vector2d>& points2D,
    Eigen::Quaterniond& virtual_from_real_rotation,
    std::vector<Eigen::Vector3d>& virtual_from_real_translations,
    std::vector<Eigen::Vector2d>& virtual_points2D) const {
  // Make sure camera is refractive when calling this function from outside.
  CHECK(camera.IsCameraRefractive()) << "Camera is not refractive";

  const size_t num_points = points2D.size();

  virtual_from_real_translations.reserve(num_points);
  virtual_points2D.reserve(num_points);

  // Compute virtual camera intrinsics.
  Camera virtual_camera = camera.VirtualCamera();

  // Compute the Refraction Axis.
  Eigen::Vector3d refrac_axis;

  if (camera.RefracModelName() == "FLATPORT") {
    // For flatports, the refraction axis is the interface normal.
    refrac_axis.x() = camera.RefracParams()[0];
    refrac_axis.y() = camera.RefracParams()[1];
    refrac_axis.z() = camera.RefracParams()[2];
    refrac_axis.normalize();
  } else if (camera.RefracModelName() == "DOMEPORT") {
    // For domeports, the refraction axis is the decentering direction.
    refrac_axis = Eigen::Vector3d(camera.RefracParams()[0],
                                  camera.RefracParams()[1],
                                  camera.RefracParams()[2])
                      .normalized();
  } else {
    LOG(FATAL) << "Refractive camera model does not exist";
  }

  virtual_from_real_rotation =
      Eigen::Quaterniond::FromTwoVectors(refrac_axis, Eigen::Vector3d::UnitZ())
          .normalized();

  // Go through every point and compute virtual camera centers.
  for (size_t i = 0; i < num_points; ++i) {
    const colmap::Ray3D ray = camera.CamFromImgRefrac(points2D[i]);
    Eigen::Vector3d virtual_camera_center;

    CHECK(IntersectLinesWithTolerance<double>(Eigen::Vector3d::Zero(),
                                              refrac_axis,
                                              ray.ori,
                                              -ray.dir,
                                              virtual_camera_center))
        << "ERROR: No intersection with the Refraction Axis";

    Eigen::Vector3d virtual_from_real_translation =
        virtual_from_real_rotation * -virtual_camera_center;
    colmap::Rigid3d virtual_from_real(virtual_from_real_rotation,
                                      virtual_from_real_translation);

    // Compute virtual image points.

    Eigen::Vector3d point3D = ray.At(1.0);

    // Now project the 3D point onto the virtual camera
    const Eigen::Vector3d point3D_v = virtual_from_real * point3D;
    const Eigen::Vector2d point2D_v =
        virtual_camera.ImgFromCam(point3D_v.hnormalized());

    // Finish up.
    virtual_points2D.push_back(point2D_v);
    virtual_from_real_translations.push_back(virtual_from_real_translation);
  }
}

void Image::ComputeVirtualTransformations(const Camera& camera) {
  std::vector<Eigen::Vector2d> points2D_xy;
  points2D_xy.reserve(NumPoints2D());
  for (point2D_t point2D_idx = 0; point2D_idx; ++point2D_idx) {
    points2D_xy.emplace_back(points2D_[point2D_idx].xy);
  }

  std::vector<Eigen::Vector2d> virtual_points2D_xy;
  ComputeVirtualTransformations(camera,
                                points2D_xy,
                                virtual_from_real_rotation_,
                                virtual_from_real_translations_,
                                virtual_points2D_xy);

  SetVirtualPoints2D(virtual_points2D_xy);
}

}  // namespace colmap
