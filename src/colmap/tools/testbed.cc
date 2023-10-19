#include "colmap/geometry/rigid3.h"
#include "colmap/scene/reconstruction.h"
#include "colmap/util/logging.h"

using namespace colmap;

int main(int argc, char* argv[]) {
  if (true) {
    // For David, add additional images to the reconstruction.
    const std::string input_path =
        "/data2/mshe/omv_src/colmap-project/dataset/2023-08_AL-Daycruise/"
        "2023-08-10_Alkor_0001_GMR_PFM-109_AUV-LUISE_Mission-305/"
        "reconstruct_last_100/result/exp1/sparse/0/";
    const std::string database_path =
        "/data2/mshe/omv_src/colmap-project/dataset/2023-08_AL-Daycruise/"
        "2023-08-10_Alkor_0001_GMR_PFM-109_AUV-LUISE_Mission-305/"
        "reconstruct_last_100/result/database.db";
    const std::string output_path =
        "/data2/mshe/omv_src/colmap-project/dataset/2023-08_AL-Daycruise/"
        "2023-08-10_Alkor_0001_GMR_PFM-109_AUV-LUISE_Mission-305/"
        "reconstruct_last_100/result/exp1/for_david/sparse/";

    Rigid3d prior_from_cam(
        Eigen::Quaterniond(0.711987, -0.00218027, -0.00757204, 0.702149),
        Eigen::Vector3d(0.347714, 0.0330715, -0.529309));

    Rigid3d cam_from_prior = Inverse(prior_from_cam);

    Reconstruction recon;
    recon.Read(input_path);

    Database database(database_path);

    std::vector<Image> images = database.ReadAllImages();

    std::cout << "Current reconstruction: " << std::endl;
    std::cout << "Number of registered images: " << recon.NumRegImages()
              << std::endl;

    // Extra images I want to manually add:
    std::vector<image_t> extra_image_ids;
    extra_image_ids.resize(42);
    std::iota(extra_image_ids.begin(), extra_image_ids.end(), 40);

    std::cout << "Extra image ids: ";
    for (image_t id : extra_image_ids) {
      std::cout << id << " ";
    }
    std::cout << std::endl;

    for (const auto image_id : extra_image_ids) {
      Image image_db = database.ReadImage(image_id);

      const Rigid3d prior_from_world = image_db.CamFromWorldPrior();
      const Rigid3d cam_from_world_prior = cam_from_prior * prior_from_world;

      image_db.CamFromWorld() = cam_from_world_prior;
      image_db.SetRegistered(true);
      recon.AddImage(image_db);
    }

    recon.WriteText(output_path);
  }

  return true;
}