#include "colmap/controllers/hybrid_mapper.h"

#include "colmap/util/misc.h"

namespace colmap {
namespace {
void AdjustGlobalBundle(const IncrementalMapperOptions& options,
                        IncrementalMapper* mapper) {
  BundleAdjustmentOptions custom_ba_options = options.GlobalBundleAdjustment();

  const size_t num_reg_images = mapper->GetReconstruction().NumRegImages();

  // Use stricter convergence criteria for first registered images.
  const size_t kMinNumRegImagesForFastBA = 10;
  if (num_reg_images < kMinNumRegImagesForFastBA) {
    custom_ba_options.solver_options.function_tolerance /= 10;
    custom_ba_options.solver_options.gradient_tolerance /= 10;
    custom_ba_options.solver_options.parameter_tolerance /= 10;
    custom_ba_options.solver_options.max_num_iterations *= 2;
    custom_ba_options.solver_options.max_linear_solver_iterations = 200;
  }

  if (options.ba_refine_intrin_after_num_images > 0 &&
      num_reg_images >
          static_cast<size_t>(options.ba_refine_intrin_after_num_images)) {
    custom_ba_options.refine_focal_length = true;
    custom_ba_options.refine_principal_point = true;
    custom_ba_options.refine_extra_params = true;
  }

  if (num_reg_images <
      static_cast<size_t>(options.ba_refine_prior_from_cam_after_num_images)) {
    custom_ba_options.refine_prior_from_cam = false;
  }

  PrintHeading1("Global bundle adjustment");
  mapper->AdjustGlobalBundle(options.Mapper(), custom_ba_options);
}

void IterativeGlobalRefinement(const IncrementalMapperOptions& options,
                               IncrementalMapper* mapper) {
  PrintHeading1("Retriangulation");
  CompleteAndMergeTracks(options, mapper);
  std::cout << "  => Retriangulated observations: "
            << mapper->Retriangulate(options.Triangulation()) << std::endl;

  for (int i = 0; i < options.ba_global_max_refinements; ++i) {
    const size_t num_observations =
        mapper->GetReconstruction().ComputeNumObservations();
    size_t num_changed_observations = 0;
    AdjustGlobalBundle(options, mapper);
    num_changed_observations += CompleteAndMergeTracks(options, mapper);
    num_changed_observations += FilterPoints(options, mapper);
    const double changed =
        num_observations == 0
            ? 0
            : static_cast<double>(num_changed_observations) / num_observations;
    std::cout << StringPrintf("  => Changed observations: %.6f", changed)
              << std::endl;
    if (changed < options.ba_global_max_refinement_change) {
      break;
    }
  }

  FilterImages(options, mapper);
}

}  // namespace

bool HybridMapperController::Options::Check() const {
  CHECK_OPTION_GE(num_workers, -1);
  CHECK_OPTION_GE(max_num_weak_area_revisit, 0);
  CHECK_OPTION_GT(re_max_num_images, 0);
  CHECK_OPTION_GT(re_max_distance, 0);
  CHECK_OPTION_GE(pgo_rel_pose_multi, 0);
  CHECK_OPTION_GE(pgo_abs_pose_multi, 0);
  CHECK_OPTION_GE(pgo_smooth_multi, 0);
  clustering_options.Check();
  CHECK_EQ(clustering_options.branching, 2);
  incremental_options.Check();
  return true;
}

HybridMapperController::HybridMapperController(
    const Options& options,
    std::shared_ptr<ReconstructionManager> reconstruction_manager)
    : options_(options),
      reconstruction_manager_(std::move(reconstruction_manager)) {
  CHECK(options_.Check());
}

void HybridMapperController::Run() {
  if (!LoadDatabase()) {
    return;
  }

  // Initialize a global reconstruction by the pose priors.
}

bool HybridMapperController::LoadDatabase() {
  PrintHeading1("Loading database");

  std::unordered_set<std::string> image_names;
  Database database(database_path_);
  Timer timer;
  timer.Start();
  const size_t min_num_matches = static_cast<size_t>(options_->min_num_matches);
  database_cache_ = DatabaseCache::Create(
      database, min_num_matches, options_->ignore_watermarks, image_names);
  std::cout << std::endl;
  timer.PrintMinutes();

  std::cout << std::endl;

  if (database_cache_->NumImages() == 0) {
    std::cout << "WARNING: No images with matches found in the database."
              << std::endl
              << std::endl;
    return false;
  }

  return true;
}

}  // namespace colmap