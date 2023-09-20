#include "colmap/sfm/hybrid_mapper.h"

#include "colmap/estimators/pose_graph_optimizer.h"
#include "colmap/estimators/triangulation.h"
#include "colmap/util/misc.h"
#include "colmap/util/threading.h"

namespace colmap {

bool HybridMapper::Options::Check() const {
  CHECK_OPTION_GT(re_max_num_images, 0);
  CHECK_OPTION_GT(re_max_distance, 0);
  CHECK_OPTION_GE(pgo_rel_pose_multi, 0);
  CHECK_OPTION_GE(pgo_abs_pose_multi, 0);
  CHECK_OPTION_GE(pgo_smooth_multi, 0);
  return true;
}

HybridMapper::HybridMapper(
    std::shared_ptr<const IncrementalMapperOptions> incremental_options,
    std::shared_ptr<const DatabaseCache> database_cache,
    const std::string& database_path,
    const std::string& image_path)
    : incremental_options_(std::move(incremental_options)),
      database_cache_(std::move(database_cache)),
      database_path_(database_path),
      image_path_(image_path),
      reconstruction_(nullptr),
      scene_clustering_(nullptr) {}

void HybridMapper::BeginReconstruction(
    const std::shared_ptr<Reconstruction>& reconstruction) {
  CHECK(reconstruction_ == nullptr);
  reconstruction_ = reconstruction;
  reconstruction_->Load(*database_cache_);
  reconstruction_->SetUp(database_cache_->CorrespondenceGraph());

  // Initialize all camera poses as their pose priors.
  const Rigid3d cam_from_prior = Inverse(reconstruction_->PriorFromCam());
  for (const auto& image_el : reconstruction_->Images()) {
    Image& image = reconstruction_->Image(image_el.first);
    const Rigid3d cam_from_world_prior =
        cam_from_prior * image.CamFromWorldPrior();
    image.CamFromWorld() = cam_from_world_prior;

    reconstruction_->RegisterImage(image_el.first);
    num_registrations_.emplace(image_el.first, 0);
  }

  // Read view graph stats from database.
  upgraded_image_pair_stats_.clear();
  image_pair_stats_.clear();

  for (const auto& image_pair : reconstruction_->ImagePairs()) {
    image_pair_stats_.emplace(image_pair.first,
                              image_pair.second.num_total_corrs);
  }
}

void HybridMapper::EndReconstruction() {
  CHECK_NOTNULL(reconstruction_);

  reconstruction_->TearDown();
  reconstruction_ = nullptr;

  num_registrations_.clear();
  image_pair_stats_.clear();
  upgraded_image_pair_stats_.clear();
  image_id_to_name_.clear();
}

const std::shared_ptr<const Reconstruction>& HybridMapper::GetReconstruction()
    const {
  CHECK_NOTNULL(reconstruction_);
  return reconstruction_;
}

void HybridMapper::PartitionScene(
    const SceneClustering::Options& clustering_options) {
  const Database database(database_path_);

  std::cout << "Reading images..." << std::endl;
  const auto images = database.ReadAllImages();
  std::unordered_map<image_t, std::string> image_id_to_name;
  image_id_to_name.reserve(images.size());
  for (const auto& image : images) {
    image_id_to_name.emplace(image.ImageId(), image.Name());
  }

  scene_clustering_ = std::make_unique<SceneClustering>(
      SceneClustering::Create(clustering_options, database));
}
}  // namespace colmap