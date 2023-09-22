#pragma once

#include "colmap/geometry/rigid3.h"
#include "colmap/scene/reconstruction.h"

#include <memory>
#include <unordered_set>

#include <Eigen/Core>
#include <ceres/ceres.h>

namespace colmap {

class PoseGraphOptimizer {
 public:
  PoseGraphOptimizer(const std::shared_ptr<Reconstruction>& reconstruction);

  void AddAbsolutePose(image_t image_id,
                       const Rigid3d& tform_measured,
                       const Eigen::Matrix6d& information,
                       ceres::LossFunction* loss_function);

  void AddRelativePose(image_t image_id_a,
                       image_t image_id_b,
                       const Rigid3d& b_from_a_measured,
                       const Eigen::Matrix6d& information,
                       ceres::LossFunction* loss_function);

  void AddSmoothMotion(image_t image_id_a,
                       image_t image_id_b,
                       image_t image_id_c,
                       const Eigen::Matrix6d& information,
                       ceres::LossFunction* loss_function);
  bool Solve();

 protected:
  std::unique_ptr<ceres::Problem> problem_;
  ceres::Solver::Summary summary_;
  std::shared_ptr<Reconstruction> reconstruction_;
};

}  // namespace colmap