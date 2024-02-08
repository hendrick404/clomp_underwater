#include "colmap/estimators/cost_functions.h"
#include "colmap/estimators/generalized_pose.h"
#include "colmap/estimators/two_view_geometry.h"
#include "colmap/geometry/essential_matrix.h"
#include "colmap/geometry/pose.h"
#include "colmap/geometry/rigid3.h"
#include "colmap/geometry/triangulation.h"
#include "colmap/math/math.h"
#include "colmap/math/random.h"
#include "colmap/scene/camera.h"
#include "colmap/scene/projection.h"

#include <fstream>

#include <Eigen/LU>
#include <Eigen/SVD>
#include <unsupported/Eigen/KroneckerProduct>

using namespace colmap;

struct PointsData {
  std::vector<Eigen::Vector2d> points2D1;
  std::vector<Eigen::Vector2d> points2D1_refrac;
  std::vector<Eigen::Vector2d> points2D2;
  std::vector<Eigen::Vector2d> points2D2_refrac;

  // Store virtual cameras here.
  std::vector<Camera> virtual_cameras1;
  std::vector<Camera> virtual_cameras2;
  std::vector<Rigid3d> virtual_from_reals1;
  std::vector<Rigid3d> virtual_from_reals2;
  Camera best_fit_camera;

  colmap::Rigid3d cam2_from_cam1_gt;
};

enum class RelTwoViewMethod {
  kNonRefrac = -1,
  kIgnore = 0,
  kIgnoreRefine = 1,
  kGR6P = 2,
  kGR6PRefine = 3,
  kBestFit = 4,
  kBestFitRefine = 5,
};

void GenerateRandomSecondViewPose(const Eigen::Vector3d& proj_center,
                                  const double distance,
                                  Rigid3d& cam2_from_cam1) {
  Eigen::Vector3d target(0.0, 0.0, distance);
  Eigen::Vector3d zaxis = (target - proj_center).normalized();
  Eigen::Vector3d yaxis(0, -1, 0);
  Eigen::Vector3d xaxis = (zaxis.cross(yaxis)).normalized();
  yaxis = (zaxis.cross(xaxis)).normalized();

  Eigen::Matrix3d R;
  R.col(0) = xaxis;
  R.col(1) = yaxis;
  R.col(2) = zaxis;
  Rigid3d cam1_from_cam2(Eigen::Quaterniond(R), proj_center);
  cam2_from_cam1 = Inverse(cam1_from_cam2);
}

void GenerateRandom2D2DPoints(const Camera& camera,
                              size_t num_points,
                              const Rigid3d& cam2_from_cam1_gt,
                              PointsData& points_data,
                              double noise_level,
                              double inlier_ratio) {
  // Generate refractive image points first, because flatport reduces FOV.

  points_data.points2D1.clear();
  points_data.points2D1_refrac.clear();
  points_data.points2D2.clear();
  points_data.points2D2_refrac.clear();

  points_data.points2D1.reserve(num_points);
  points_data.points2D1_refrac.reserve(num_points);
  points_data.points2D2.reserve(num_points);
  points_data.points2D2_refrac.reserve(num_points);

  points_data.cam2_from_cam1_gt = cam2_from_cam1_gt;

  size_t num_inliers =
      static_cast<size_t>(static_cast<double>(num_points) * inlier_ratio);
  size_t cnt = 0;
  while (true) {
    if (cnt >= num_points) {
      break;
    }
    Eigen::Vector2d point2D1_refrac;
    point2D1_refrac.x() =
        RandomUniformReal(0.5, static_cast<double>(camera.width) - 0.5);
    point2D1_refrac.y() =
        RandomUniformReal(0.5, static_cast<double>(camera.height) - 0.5);

    Ray3D ray_refrac = camera.CamFromImgRefrac(point2D1_refrac);

    const double d = RandomUniformReal(6.0, 8.0);

    // Now, do projection.
    Eigen::Vector3d point3D1 = ray_refrac.At(d);
    Eigen::Vector3d point3D2 = cam2_from_cam1_gt * point3D1;

    Eigen::Vector2d point2D2_refrac = camera.ImgFromCamRefrac(point3D2);

    if (std::isnan(point2D2_refrac.x()) || std::isnan(point2D2_refrac.y())) {
      continue;
    }

    if (point2D2_refrac.x() < 0 || point2D2_refrac.x() > camera.width ||
        point2D2_refrac.y() < 0 || point2D2_refrac.y() > camera.height) {
      continue;
    }

    Eigen::Vector2d point2D1 = camera.ImgFromCam(point3D1.hnormalized());
    Eigen::Vector2d point2D2 = camera.ImgFromCam(point3D2.hnormalized());

    if (cnt < num_inliers) {
      // Add noise to the points.
      if (noise_level > 0) {
        point2D1.x() += RandomGaussian(0.0, noise_level);
        point2D1.y() += RandomGaussian(0.0, noise_level);
        point2D1_refrac.x() += RandomGaussian(0.0, noise_level);
        point2D1_refrac.y() += RandomGaussian(0.0, noise_level);

        point2D2.x() += RandomGaussian(0.0, noise_level);
        point2D2.y() += RandomGaussian(0.0, noise_level);
        point2D2_refrac.x() += RandomGaussian(0.0, noise_level);
        point2D2_refrac.y() += RandomGaussian(0.0, noise_level);
      }
    } else {
      // Add huge noise to the points, this should be an outlier point.
      point2D1.x() += RandomGaussian(0.0, 200.0);
      point2D1.y() += RandomGaussian(0.0, 200.0);
      point2D1_refrac.x() += RandomGaussian(0.0, 200.0);
      point2D1_refrac.y() += RandomGaussian(0.0, 200.0);

      point2D2.x() += RandomGaussian(0.0, 200.0);
      point2D2.y() += RandomGaussian(0.0, 200.0);
      point2D2_refrac.x() += RandomGaussian(0.0, 200.0);
      point2D2_refrac.y() += RandomGaussian(0.0, 200.0);
    }

    points_data.points2D1.push_back(point2D1);
    points_data.points2D2.push_back(point2D2);
    points_data.points2D1_refrac.push_back(point2D1_refrac);
    points_data.points2D2_refrac.push_back(point2D2_refrac);

    cnt++;
  }

  camera.ComputeVirtuals(points_data.points2D1_refrac,
                         points_data.virtual_cameras1,
                         points_data.virtual_from_reals1);
  camera.ComputeVirtuals(points_data.points2D2_refrac,
                         points_data.virtual_cameras2,
                         points_data.virtual_from_reals2);

  const double kApproxDepth = 5.0;
  points_data.best_fit_camera =
      BestFitNonRefracCamera(CameraModelId::kOpenCV, camera, kApproxDepth);
}

TwoViewGeometry EstimateRefractiveTwoViewIgnoreRefraction(
    const Camera& camera1,
    const std::vector<Eigen::Vector2d>& points1,
    const std::vector<Camera>& virtual_cameras1,
    const std::vector<Rigid3d>& virtual_from_reals1,
    const Camera& camera2,
    const std::vector<Eigen::Vector2d>& points2,
    const std::vector<Camera>& virtual_cameras2,
    const std::vector<Rigid3d>& virtual_from_reals2,
    const FeatureMatches& matches,
    const TwoViewGeometryOptions& options,
    bool refine = false) {
  // Perform essential matrix estimation
  TwoViewGeometry geometry = EstimateCalibratedTwoViewGeometry(
      camera1, points1, camera2, points2, matches, options);

  if (refine) {
    // Perform refinement afterwards:
    std::vector<Eigen::Vector2d> inlier_points1_normalized;
    std::vector<Eigen::Vector2d> inlier_points2_normalized;

    std::vector<Rigid3d> inlier_virtual_from_reals1;
    std::vector<Rigid3d> inlier_virtual_from_reals2;

    inlier_points1_normalized.reserve(geometry.inlier_matches.size());
    inlier_points2_normalized.reserve(geometry.inlier_matches.size());

    inlier_virtual_from_reals1.reserve(geometry.inlier_matches.size());
    inlier_virtual_from_reals2.reserve(geometry.inlier_matches.size());

    for (const auto& match : geometry.inlier_matches) {
      const Camera& virtual_camera1 = virtual_cameras1[match.point2D_idx1];
      const Camera& virtual_camera2 = virtual_cameras2[match.point2D_idx2];

      inlier_points1_normalized.push_back(
          virtual_camera1.CamFromImg(points1[match.point2D_idx1]));
      inlier_points2_normalized.push_back(
          virtual_camera2.CamFromImg(points2[match.point2D_idx2]));

      inlier_virtual_from_reals1.push_back(
          virtual_from_reals1[match.point2D_idx1]);
      inlier_virtual_from_reals2.push_back(
          virtual_from_reals2[match.point2D_idx2]);
    }

    if (!RefineRefractiveTwoViewGeometry(inlier_points1_normalized,
                                         inlier_virtual_from_reals1,
                                         inlier_points2_normalized,
                                         inlier_virtual_from_reals2,
                                         &geometry.cam2_from_cam1)) {
      // Optimization failed, directly return and clean up the inlier matches.
      geometry.inlier_matches.clear();
      return geometry;
    };
  }
  return geometry;
}

size_t EstimateRelativePose(Camera& camera,
                            const PointsData& points_data,
                            Rigid3d& cam2_from_cam1,
                            RelTwoViewMethod method_id) {
  size_t num_points = points_data.points2D1.size();
  size_t num_inliers = 0;

  TwoViewGeometryOptions two_view_geometry_options;
  two_view_geometry_options.compute_relative_pose = true;
  two_view_geometry_options.ransac_options.max_error = 4.0;

  FeatureMatches matches;
  matches.reserve(num_points);

  for (size_t i = 0; i < num_points; i++) {
    matches.emplace_back(i, i);
  }

  TwoViewGeometry two_view_geometry;

  switch (method_id) {
    case RelTwoViewMethod::kNonRefrac:
      two_view_geometry =
          EstimateCalibratedTwoViewGeometry(camera,
                                            points_data.points2D1,
                                            camera,
                                            points_data.points2D2,
                                            matches,
                                            two_view_geometry_options);
      break;
    case RelTwoViewMethod::kGR6P:
      two_view_geometry_options.compute_relative_pose = true;
      two_view_geometry =
          EstimateRefractiveTwoViewGeometry(points_data.points2D1_refrac,
                                            points_data.virtual_cameras1,
                                            points_data.virtual_from_reals1,
                                            points_data.points2D2_refrac,
                                            points_data.virtual_cameras2,
                                            points_data.virtual_from_reals2,
                                            matches,
                                            two_view_geometry_options,
                                            false);
      break;
    case RelTwoViewMethod::kGR6PRefine:
      two_view_geometry_options.compute_relative_pose = true;
      two_view_geometry =
          EstimateRefractiveTwoViewGeometry(points_data.points2D1_refrac,
                                            points_data.virtual_cameras1,
                                            points_data.virtual_from_reals1,
                                            points_data.points2D2_refrac,
                                            points_data.virtual_cameras2,
                                            points_data.virtual_from_reals2,
                                            matches,
                                            two_view_geometry_options,
                                            true);
      break;
    case RelTwoViewMethod::kBestFit:

      two_view_geometry = EstimateRefractiveTwoViewGeometryUseBestFit(
          points_data.best_fit_camera,
          points_data.points2D1_refrac,
          points_data.virtual_cameras1,
          points_data.virtual_from_reals1,
          points_data.best_fit_camera,
          points_data.points2D2_refrac,
          points_data.virtual_cameras2,
          points_data.virtual_from_reals2,
          matches,
          two_view_geometry_options,
          false);
      break;
    case RelTwoViewMethod::kBestFitRefine:
      two_view_geometry = EstimateRefractiveTwoViewGeometryUseBestFit(
          points_data.best_fit_camera,
          points_data.points2D1_refrac,
          points_data.virtual_cameras1,
          points_data.virtual_from_reals1,
          points_data.best_fit_camera,
          points_data.points2D2_refrac,
          points_data.virtual_cameras2,
          points_data.virtual_from_reals2,
          matches,
          two_view_geometry_options,
          true);
      break;
    case RelTwoViewMethod::kIgnore:
      two_view_geometry = EstimateRefractiveTwoViewIgnoreRefraction(
          camera,
          points_data.points2D1_refrac,
          points_data.virtual_cameras1,
          points_data.virtual_from_reals1,
          camera,
          points_data.points2D2_refrac,
          points_data.virtual_cameras2,
          points_data.virtual_from_reals2,
          matches,
          two_view_geometry_options,
          false);
      break;
    case RelTwoViewMethod::kIgnoreRefine:
      two_view_geometry = EstimateRefractiveTwoViewIgnoreRefraction(
          camera,
          points_data.points2D1_refrac,
          points_data.virtual_cameras1,
          points_data.virtual_from_reals1,
          camera,
          points_data.points2D2_refrac,
          points_data.virtual_cameras2,
          points_data.virtual_from_reals2,
          matches,
          two_view_geometry_options,
          true);
      break;
    default:
      LOG(ERROR) << "Relative two-view method does not exist!";
      break;
  }

  cam2_from_cam1 = two_view_geometry.cam2_from_cam1;
  num_inliers = two_view_geometry.inlier_matches.size();
  return num_inliers;
}

void RelativePoseError(const colmap::Rigid3d& cam2_from_cam1_gt,
                       const colmap::Rigid3d& cam2_from_cam1_est,
                       double& rotation_error,
                       double& angular_error) {
  colmap::Rigid3d cam2_from_cam1_gt_norm = cam2_from_cam1_gt;
  cam2_from_cam1_gt_norm.translation.normalize();
  colmap::Rigid3d cam2_from_cam1_est_norm = cam2_from_cam1_est;
  cam2_from_cam1_est_norm.translation.normalize();

  Eigen::Quaterniond rotation_diff =
      cam2_from_cam1_gt.rotation * cam2_from_cam1_est.rotation.inverse();
  rotation_error = colmap::RadToDeg(Eigen::AngleAxisd(rotation_diff).angle());

  double cos_theta = cam2_from_cam1_gt_norm.translation.dot(
      cam2_from_cam1_est_norm.translation);
  if (cos_theta < 0) {
    cos_theta = -cos_theta;
  }
  if (cos_theta > 1) {
    cos_theta = 1.0;
  }
  angular_error = RadToDeg(acos(cos_theta));
  // if (is_refractive) {
  //   const double baseline_est = (cam2_from_cam1_est.rotation.inverse() *
  //                                -cam2_from_cam1_est.translation)
  //                                   .norm();
  //   const double baseline_gt =
  //       (cam2_from_cam1_gt.rotation.inverse() *
  //       -cam2_from_cam1_gt.translation)
  //           .norm();
  //   scale_error = (1 - baseline_est / baseline_gt);
  // } else {
  //   scale_error = 0.0;
  // }
}

void EvaluateMultipleMethods(colmap::Camera& camera,
                             size_t num_points,
                             size_t num_exps,
                             double inlier_ratio,
                             const std::string& output_path) {
  std::vector<double> noise_levels = {0.0, 0.2, 0.5, 0.8, 1.2, 1.5, 1.8, 2.0};
  // std::vector<double> noise_levels = {0.0};

  std::ofstream file(output_path, std::ios::out);

  for (const double& noise : noise_levels) {
    std::cout << "Noise level: " << noise << std::endl;

    // Generate random datasets first.
    std::vector<PointsData> datasets;
    datasets.reserve(num_exps);

    std::cout << "Generating random data ..." << std::endl;
    for (size_t i = 0; i < num_exps; i++) {
      // Create a random GT pose.
      const double tx = colmap::RandomUniformReal(8.0, 10.0);
      const double ty = colmap::RandomUniformReal(-2.5, 2.5);
      const double tz = colmap::RandomUniformReal(-2.5, 2.5);
      const double distance = colmap::RandomUniformReal(6.0, 8.0);
      Eigen::Vector3d proj_center(tx, ty, tz);

      colmap::Rigid3d cam2_from_cam1;

      GenerateRandomSecondViewPose(proj_center, distance, cam2_from_cam1);
      PointsData points_data;
      GenerateRandom2D2DPoints(
          camera, num_points, cam2_from_cam1, points_data, noise, inlier_ratio);
      datasets.push_back(points_data);
    }

    std::vector<RelTwoViewMethod> methods = {RelTwoViewMethod::kNonRefrac,
                                             RelTwoViewMethod::kBestFit};

    // std::vector<RelTwoViewMethod> methods = {RelTwoViewMethod::kNonRefrac,
    //                                          RelTwoViewMethod::kBestFitRefine};

    // Evaluate random dataset.
    // Errors of all methods:
    std::vector<double> rot_error_mean_all_methods;
    std::vector<double> rot_error_std_all_methods;
    std::vector<double> ang_error_mean_all_methods;
    std::vector<double> ang_error_std_all_methods;
    std::vector<double> inlier_ratio_all_methods;

    for (const RelTwoViewMethod& method_id : methods) {
      // Rotation error (degrees).
      std::vector<double> rotation_errors;
      // Angular error (error between the estimated translation direction and
      // ground truth direction)/
      std::vector<double> angular_errors;
      // Inlier ratio.
      std::vector<double> inlier_ratios;

      double time = 0.0;

      std::cout << "Evaluating Method: " << static_cast<int>(method_id)
                << std::endl;

      Timer timer;

      timer.Start();
      for (size_t i = 0; i < num_exps; i++) {
        const PointsData& points_data = datasets[i];
        colmap::Rigid3d cam2_from_cam1_est;
        size_t num_inliers = EstimateRelativePose(
            camera, points_data, cam2_from_cam1_est, method_id);

        double rotation_error, angular_error;
        RelativePoseError(points_data.cam2_from_cam1_gt,
                          cam2_from_cam1_est,
                          rotation_error,
                          angular_error);
        rotation_errors.push_back(rotation_error);
        angular_errors.push_back(angular_error);
        inlier_ratios.push_back(static_cast<double>(num_inliers) /
                                static_cast<double>(num_points));
      }
      timer.Pause();
      time = timer.ElapsedSeconds();

      const double rot_error_mean = Mean(rotation_errors);
      const double rot_error_std = StdDev(rotation_errors);
      const double ang_error_mean = Mean(angular_errors);
      const double ang_error_std = StdDev(angular_errors);
      const double inlier_ratio_mean = Mean(inlier_ratios);

      rot_error_mean_all_methods.push_back(rot_error_mean);
      rot_error_std_all_methods.push_back(rot_error_std);
      ang_error_mean_all_methods.push_back(ang_error_mean);
      ang_error_std_all_methods.push_back(ang_error_std);
      inlier_ratio_all_methods.push_back(inlier_ratio_mean);

      std::cout << "Relative pose error : Rotation: " << rot_error_mean
                << " +/- " << rot_error_std << " -- Angular: " << ang_error_mean
                << " +/- " << ang_error_std
                << " -- inlier ratio: " << inlier_ratio_mean
                << " GT inlier ratio: " << inlier_ratio << std::endl;
    }

    file << noise;
    for (size_t i = 0; i < methods.size(); i++) {
      file << " " << rot_error_mean_all_methods[i] << " "
           << rot_error_std_all_methods[i] << " "
           << ang_error_mean_all_methods[i] << " "
           << ang_error_std_all_methods[i] << " "
           << inlier_ratio_all_methods[i];
    }
    file << std::endl;
  }

  file.close();
}

/// Xiaohu's relative pose init ///

namespace {
// method for calculating the pseudo-Inverse as recommended by Eigen developers
template <typename _Matrix_Type_>
_Matrix_Type_ pseudoInverse(
    const _Matrix_Type_& a,
    double epsilon = std::numeric_limits<double>::epsilon()) {
  Eigen::JacobiSVD<_Matrix_Type_> svd(
      a, Eigen::ComputeFullU | Eigen::ComputeFullV);
  // For a non-square matrix
  // Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeThinU |
  // Eigen::ComputeThinV);
  double tolerance = epsilon * std::max(a.cols(), a.rows()) *
                     svd.singularValues().array().abs()(0);
  return svd.matrixV() *
         (svd.singularValues().array().abs() > tolerance)
             .select(svd.singularValues().array().inverse(), 0)
             .matrix()
             .asDiagonal() *
         svd.matrixU().adjoint();
}
}  // namespace

TwoViewGeometry XiaoHuInit(const Camera& camera,
                           const PointsData& points_data) {
  LOG(INFO) << "Evaluating XiaoHu's approach";

  TwoViewGeometry two_view_geometry;

  /// prepare the rays.
  const size_t kNumPoints = points_data.points2D1_refrac.size();

  std::vector<colmap::Ray3D> rays1;
  std::vector<colmap::Ray3D> rays2;

  rays1.reserve(kNumPoints);
  rays2.reserve(kNumPoints);

  for (const auto& point : points_data.points2D1_refrac) {
    rays1.push_back(camera.CamFromImgRefrac(point));
  }
  for (const auto& point : points_data.points2D2_refrac) {
    rays2.push_back(camera.CamFromImgRefrac(point));
  }

  // virtual camera centers.
  std::vector<Eigen::Vector3d> vc1, vc2;
  vc1.reserve(kNumPoints);
  vc2.reserve(kNumPoints);

  for (size_t i = 0; i < kNumPoints; ++i) {
    vc1.push_back(points_data.virtual_from_reals1[i].rotation.inverse() *
                  -points_data.virtual_from_reals1[i].translation);
    vc2.push_back(points_data.virtual_from_reals2[i].rotation.inverse() *
                  -points_data.virtual_from_reals2[i].translation);
  }

  Eigen::MatrixXd A(kNumPoints, 18);
  A.setZero();

  for (size_t i = 0; i < kNumPoints; ++i) {
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    const Eigen::Vector3d r1 = rays1[i].dir;
    const Eigen::Matrix<double, 1, 3> r1_t = rays1[i].dir.transpose();
    const Eigen::Matrix<double, 1, 3> r2_t = rays2[i].dir.transpose();

    const Eigen::Vector3d& cv1 = vc1[i];
    const Eigen::Vector3d& cv2 = vc2[i];
    Eigen::Matrix3d cv1_hat = CrossProductMatrix(cv1);
    Eigen::Matrix3d cv2_hat = CrossProductMatrix(cv2);

    // Compute the first kronecker product.
    Eigen::KroneckerProduct<Eigen::Matrix<double, 1, 3>,
                            Eigen::Matrix<double, 3, 3>>
        kron1 = Eigen::KroneckerProduct<Eigen::Matrix<double, 1, 3>,
                                        Eigen::Matrix<double, 3, 3>>(r1_t, I);

    // Compute the second kronecker product.
    Eigen::Matrix<double, 1, 3> cv1_u = (cv1_hat * r1).transpose();
    Eigen::KroneckerProduct<Eigen::Matrix<double, 1, 3>,
                            Eigen::Matrix<double, 3, 3>>
        kron2 = Eigen::KroneckerProduct<Eigen::Matrix<double, 1, 3>,
                                        Eigen::Matrix<double, 3, 3>>(cv1_u, I);

    A.block<1, 9>(i, 0) = r2_t * cv2_hat * kron1 - r2_t * kron2;
    A.block<1, 9>(i, 9) = -r2_t * kron1;
  }

  Eigen::MatrixXd AR = A.block(0, 0, kNumPoints, 9);
  Eigen::MatrixXd AE = A.block(0, 9, kNumPoints, 9);
  Eigen::MatrixXd ARpinv = AR.completeOrthogonalDecomposition().pseudoInverse();
  // Eigen::MatrixXd ARpinv = pseudoInverse(AR);
  Eigen::MatrixXd eye_N(kNumPoints, kNumPoints);
  eye_N.setIdentity();

  Eigen::MatrixXd B = (AR * ARpinv - eye_N) * AE;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(B, Eigen::ComputeFullV);

  Eigen::MatrixXd sol = svd.matrixV().col(8);
  const Eigen::Map<const Eigen::Matrix3d> E_raw(sol.data());

  // Enforcing the internal constraint that two singular values must be equal
  // and one must be zero.
  Eigen::JacobiSVD<Eigen::Matrix3d> E_raw_svd(
      E_raw, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::Vector3d singular_values = E_raw_svd.singularValues();
  singular_values(0) = (singular_values(0) + singular_values(1)) / 2.0;
  singular_values(1) = singular_values(0);
  singular_values(2) = 0.0;
  const Eigen::Matrix3d E = E_raw_svd.matrixU() * singular_values.asDiagonal() *
                            E_raw_svd.matrixV().transpose();

  // Solve kernel part (Basically repeat PoseFromEssentialMatrix, but replace
  // proj_center by virtual_proj_center).
  Eigen::Matrix3d R_cam1_from_cam2;
  Eigen::Vector3d t_cam1_from_cam2;
  std::vector<Eigen::Vector3d> points3D;
  Eigen::Matrix3d R1;
  Eigen::Matrix3d R2;
  Eigen::Vector3d t;
  DecomposeEssentialMatrix(E, &R1, &R2, &t);
  points3D.clear();

  // Generate all possible projection matrix combinations.
  const std::array<Eigen::Matrix3d, 4> R_cmbs{{R1, R2, R1, R2}};
  const std::array<Eigen::Vector3d, 4> t_cmbs{{t, t, -t, -t}};

  auto CalculateDepth = [](const Eigen::Matrix3x4d& cam_from_world,
                           const Eigen::Vector3d& point3D) {
    const double proj_z = cam_from_world.row(2).dot(point3D.homogeneous());
    return proj_z * cam_from_world.col(2).norm();
  };

  for (size_t i = 0; i < R_cmbs.size(); ++i) {
    std::vector<Eigen::Vector3d> points3D_cmb;
    // Check cheriality here

    const double kMinDepth = std::numeric_limits<double>::epsilon();
    const double max_depth =
        1000.0f * (R_cmbs[i].transpose() * t_cmbs[i]).norm();

    for (size_t j = 0; j < kNumPoints; ++j) {
      const Rigid3d& virtual_from_real1 = points_data.virtual_from_reals1[j];
      const Rigid3d& virtual_from_real2 = points_data.virtual_from_reals2[j];

      const Rigid3d world_from_cam2(Eigen::Quaterniond(R_cmbs[i].transpose()),
                                    t_cmbs[i]);
      const Rigid3d cam2_from_world = Inverse(world_from_cam2);
      const Rigid3d virtual2_from_world = virtual_from_real2 * cam2_from_world;

      const Eigen::Matrix3x4d virtual_proj_matrix1 =
          virtual_from_real1.ToMatrix();
      const Eigen::Matrix3x4d virtual_proj_matrix2 =
          virtual2_from_world.ToMatrix();

      const Eigen::Vector3d point3D =
          TriangulatePoint(virtual_proj_matrix1,
                           virtual_proj_matrix2,
                           rays1[j].dir.hnormalized(),
                           rays2[j].dir.hnormalized());

      const double depth1 = CalculateDepth(virtual_proj_matrix1, point3D);
      if (depth1 > kMinDepth && depth1 < max_depth) {
        const double depth2 = CalculateDepth(virtual_proj_matrix2, point3D);
        if (depth2 > kMinDepth && depth2 < max_depth) {
          points3D_cmb.push_back(point3D);
        }
      }
    }

    if (points3D_cmb.size() >= points3D.size()) {
      R_cam1_from_cam2 = R_cmbs[i].transpose();
      t_cam1_from_cam2 = t_cmbs[i];
      points3D = points3D_cmb;
    }
  }

  // Solve for t

  A.resize(kNumPoints, 3);
  A.setZero();

  Eigen::VectorXd b(kNumPoints);
  b.setZero();

  for (size_t i = 0; i < kNumPoints; ++i) {
    A.row(i) = rays2[i].dir.transpose() * R_cam1_from_cam2.transpose() *
               CrossProductMatrix(rays1[i].dir);
    // b(i,:) =
    // r22(:,i)'*skewm(t2(:,i))*R'*r12(:,i)-r22(:,i)'*R'*skewm(t1(:,i))*r12(:,i);
    b.row(i) = rays2[i].dir.transpose() * CrossProductMatrix(vc2[i]) *
                   R_cam1_from_cam2.transpose() * rays1[i].dir -
               rays2[i].dir.transpose() * R_cam1_from_cam2.transpose() *
                   CrossProductMatrix(vc1[i]) * rays1[i].dir;
  }

  Eigen::JacobiSVD<Eigen::MatrixXd> t_solver_svd(
      A, Eigen::ComputeFullU | Eigen::ComputeFullV);
  t_cam1_from_cam2 = t_solver_svd.solve(b);

  // LOG(INFO) << "R1: \n" << R1;
  // LOG(INFO) << "R2: \n" << R2;
  // LOG(INFO) << "t: \n" << R2;

  Rigid3d cam1_from_cam2(Eigen::Quaterniond(R_cam1_from_cam2),
                         t_cam1_from_cam2);
  two_view_geometry.cam2_from_cam1 = Inverse(cam1_from_cam2);

  return two_view_geometry;
}
/// Xiaohu's relative pose init ///

void EvaluateXiaoHu(const Camera& camera) {
  // Create a random GT pose.
  const double tx = colmap::RandomUniformReal(8.0, 10.0);
  const double ty = colmap::RandomUniformReal(-2.5, 2.5);
  const double tz = colmap::RandomUniformReal(-2.5, 2.5);
  const double distance = colmap::RandomUniformReal(6.0, 8.0);
  Eigen::Vector3d proj_center(tx, ty, tz);

  colmap::Rigid3d cam2_from_cam1;

  GenerateRandomSecondViewPose(proj_center, distance, cam2_from_cam1);
  PointsData points_data;
  GenerateRandom2D2DPoints(camera, 100, cam2_from_cam1, points_data, 0.5, 1.0);

  TwoViewGeometry geometry = XiaoHuInit(camera, points_data);

  LOG(INFO) << "GT  : "
            << points_data.cam2_from_cam1_gt.rotation.coeffs().transpose()
            << " , " << points_data.cam2_from_cam1_gt.translation.transpose();
  LOG(INFO) << "EST : " << geometry.cam2_from_cam1.rotation.coeffs().transpose()
            << " , " << geometry.cam2_from_cam1.translation.transpose();
}

int main(int argc, char* argv[]) {
  SetPRNGSeed(time(NULL));

  // Camera parameters coming from Anton 131 map.
  // Camera camera;
  // camera.SetWidth(4104);
  // camera.SetHeight(3006);
  // camera.SetModelIdFromName("METASHAPE_FISHEYE");
  // std::vector<double> params = {1000.7964068878537,
  //                               1000.6679248258547,
  //                               2097.4832274550317,
  //                               1641.1207545881762,
  //                               0,
  //                               0,
  //                               0,
  //                               0,
  //                               0,
  //                               0};
  Camera camera;
  camera.width = 1113;
  camera.height = 835;
  camera.model_id = CameraModelId::kPinhole;
  std::vector<double> params = {400.5, 400.5, 556.5, 417.5};
  camera.params = params;

  // Flatport setup.
  camera.refrac_model_id = CameraRefracModelId::kFlatPort;
  Eigen::Vector3d int_normal;
  int_normal[0] = RandomUniformReal(-0.1, 0.1);
  int_normal[1] = RandomUniformReal(-0.1, 0.1);
  int_normal[2] = RandomUniformReal(0.9, 1.1);

  int_normal = Eigen::Vector3d::UnitZ();
  // int_normal.normalize();

  std::vector<double> flatport_params = {int_normal[0],
                                         int_normal[1],
                                         int_normal[2],
                                         0.05,
                                         0.02,
                                         1.0,
                                         1.52,
                                         1.334};
  camera.refrac_params = flatport_params;

  // camera.refrac_model_id = CameraRefracModelId::kDomePort;
  // Eigen::Vector3d decentering;
  // // decentering[0] = RandomUniformReal(-0.001, 0.001);
  // // decentering[1] = RandomUniformReal(-0.001, 0.001);
  // // decentering[2] = RandomUniformReal(-0.002, 0.002);

  // decentering[0] = 0;
  // decentering[1] = 0;
  // decentering[2] = 0;

  // std::vector<double> domeport_params = {decentering[0],
  //                                        decentering[1],
  //                                        decentering[2],
  //                                        0.05,
  //                                        0.007,
  //                                        1.0,
  //                                        1.52,
  //                                        1.334};
  // camera.refrac_params = domeport_params;

  // Generate simulated point data.
  const size_t num_points = 200;
  const double inlier_ratio = 1.0;

  std::string output_dir =
      "/home/mshe/workspace/omv_src/colmap-project/refrac_sfm_eval/plots/"
      "rel_pose/compare_methods/";
  std::stringstream ss;
  ss << output_dir
     << "/rel_pose_flat_non_ortho_far_num_points_eval_different_distances_"
        "debug_"
     << num_points << "_inlier_ratio_" << inlier_ratio << ".txt";
  std::string output_path = ss.str();

  EvaluateXiaoHu(camera);
  // EvaluateMultipleMethods(camera, num_points, 10, inlier_ratio, output_path);

  return true;
}