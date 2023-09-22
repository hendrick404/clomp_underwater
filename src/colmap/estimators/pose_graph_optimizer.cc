#include "colmap/estimators/pose_graph_optimizer.h"

#include "colmap/estimators/cost_functions.h"
#include "colmap/util/misc.h"
#include "colmap/util/threading.h"

namespace colmap {

PoseGraphOptimizer::PoseGraphOptimizer(
    const std::shared_ptr<Reconstruction>& reconstruction)
    : reconstruction_(reconstruction) {
  CHECK_NOTNULL(reconstruction_);
  CHECK(!problem_) << "Cannot use the same PoseGraphOptimizer multiple times";

  ceres::Problem::Options problem_options;
  problem_options.loss_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
  problem_ = std::make_unique<ceres::Problem>(problem_options);
}

void PoseGraphOptimizer::AddAbsolutePose(const image_t image_id,
                                         const Rigid3d& tform_measured,
                                         const Eigen::Matrix6d& information,
                                         ceres::LossFunction* loss_function) {
  CHECK(reconstruction_->IsImageRegistered(image_id));

  const Eigen::Matrix6d sqrt_information = information.llt().matrixL();
  Image& image = reconstruction_->Image(image_id);

  // CostFunction assumes unit quaternions.
  image.CamFromWorld().rotation.normalize();

  double* cam_from_world_rotation =
      image.CamFromWorld().rotation.coeffs().data();
  double* cam_from_world_translation = image.CamFromWorld().translation.data();

  ceres::CostFunction* cost_function =
      AbsolutePoseErrorCostFunction::Create(tform_measured, sqrt_information);
  problem_->AddResidualBlock(cost_function,
                             loss_function,
                             cam_from_world_rotation,
                             cam_from_world_translation);

  SetQuaternionManifold(problem_.get(), cam_from_world_rotation);
}

void PoseGraphOptimizer::AddRelativePose(const image_t image_id_a,
                                         const image_t image_id_b,
                                         const Rigid3d& b_from_a_measured,
                                         const Eigen::Matrix6d& information,
                                         ceres::LossFunction* loss_function) {
  CHECK(reconstruction_->IsImageRegistered(image_id_a) &&
        reconstruction_->IsImageRegistered(image_id_b));

  const Eigen::Matrix6d sqrt_information = information.llt().matrixL();

  Image& image_a = reconstruction_->Image(image_id_a);
  Image& image_b = reconstruction_->Image(image_id_b);

  // CostFunction assumes unit quaternions.
  image_a.CamFromWorld().rotation.normalize();
  image_b.CamFromWorld().rotation.normalize();

  double* a_from_world_rotation =
      image_a.CamFromWorld().rotation.coeffs().data();
  double* a_from_world_translation = image_a.CamFromWorld().translation.data();
  double* b_from_world_rotation =
      image_b.CamFromWorld().rotation.coeffs().data();
  double* b_from_world_translation = image_b.CamFromWorld().translation.data();

  ceres::CostFunction* cost_function =
      RelativePoseError6DoFCostFunction::Create(b_from_a_measured,
                                                sqrt_information);

  problem_->AddResidualBlock(cost_function,
                             loss_function,
                             a_from_world_rotation,
                             a_from_world_translation,
                             b_from_world_rotation,
                             b_from_world_translation);

  SetQuaternionManifold(problem_.get(), a_from_world_rotation);
  SetQuaternionManifold(problem_.get(), b_from_world_rotation);
}

void PoseGraphOptimizer::AddSmoothMotion(const image_t image_id_a,
                                         const image_t image_id_b,
                                         const image_t image_id_c,
                                         const Eigen::Matrix6d& information,
                                         ceres::LossFunction* loss_function) {
  CHECK(reconstruction_->IsImageRegistered(image_id_a) &&
        reconstruction_->IsImageRegistered(image_id_b) &&
        reconstruction_->IsImageRegistered(image_id_c));

  const Eigen::Matrix6d sqrt_information = information.llt().matrixL();

  Image& image_a = reconstruction_->Image(image_id_a);
  Image& image_b = reconstruction_->Image(image_id_b);
  Image& image_c = reconstruction_->Image(image_id_c);

  // CostFunction assumes unit quaternions.

  image_a.CamFromWorld().rotation.normalize();
  image_b.CamFromWorld().rotation.normalize();
  image_c.CamFromWorld().rotation.normalize();

  double* a_from_world_rotation =
      image_a.CamFromWorld().rotation.coeffs().data();
  double* a_from_world_translation = image_a.CamFromWorld().translation.data();
  double* b_from_world_rotation =
      image_b.CamFromWorld().rotation.coeffs().data();
  double* b_from_world_translation = image_b.CamFromWorld().translation.data();
  double* c_from_world_rotation =
      image_c.CamFromWorld().rotation.coeffs().data();
  double* c_from_world_translation = image_c.CamFromWorld().translation.data();

  ceres::CostFunction* cost_function =
      SmoothMotionCostFunction::Create(sqrt_information);

  problem_->AddResidualBlock(cost_function,
                             loss_function,
                             a_from_world_rotation,
                             a_from_world_translation,
                             b_from_world_rotation,
                             b_from_world_translation,
                             c_from_world_rotation,
                             c_from_world_translation);

  SetQuaternionManifold(problem_.get(), a_from_world_rotation);
  SetQuaternionManifold(problem_.get(), b_from_world_rotation);
  SetQuaternionManifold(problem_.get(), c_from_world_rotation);
}

bool PoseGraphOptimizer::Solve() {
  if (problem_->NumResiduals() == 0) {
    return false;
  }

  const double kEpsilon = 1e-10;
  const size_t kMaxNumIterations = 300;

  ceres::Solver::Options solver_options;
  solver_options.minimizer_progress_to_stdout = true;
  solver_options.max_num_iterations = kMaxNumIterations;
  solver_options.function_tolerance = 1e-2 * kEpsilon;
  solver_options.gradient_tolerance = 1e-2 * kEpsilon;
  solver_options.parameter_tolerance = 1e-2 * kEpsilon;
  solver_options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;

  const int kMinNumResidualsForMultiThreading = 500;

  if (problem_->NumResiduals() < kMinNumResidualsForMultiThreading) {
    solver_options.num_threads = 1;
#if CERES_VERSION_MAJOR < 2
    solver_options.num_linear_solver_threads = 1;
#endif  // CERES_VERSION_MAJOR
  } else {
    solver_options.num_threads = GetEffectiveNumThreads(-1);
#if CERES_VERSION_MAJOR < 2
    solver_options.num_linear_solver_threads = GetEffectiveNumThreads(-1);
#endif  // CERES_VERSION_MAJOR
  }

  std::string solver_error;
  CHECK(solver_options.IsValid(&solver_error)) << solver_error;

  ceres::Solve(solver_options, problem_.get(), &summary_);

  if (solver_options.minimizer_progress_to_stdout) {
    std::cout << std::endl;
  }

  PrintHeading2("Pose graph optimizer report");
  std::cout << summary_.BriefReport() << std::endl;

  return summary_.IsSolutionUsable();
}

}  // namespace colmap