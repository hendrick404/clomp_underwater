#include "colmap/geometry/rigid3.h"
#include "colmap/math/random.h"
#include "colmap/scene/reconstruction.h"
#include "colmap/util/logging.h"

using namespace colmap;

int main(int argc, char* argv[]) {
  if (false) {
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
  if (true) {
    Camera camera;
    camera.SetWidth(2048);
    camera.SetHeight(1536);
    camera.SetModelIdFromName("PINHOLE");
    std::vector<double> params = {
        1300.900000, 1300.900000, 1024.000000, 768.000000};
    camera.SetParams(params);

    // Flatport setup.
    camera.SetRefracModelIdFromName("FLATPORT");
    Eigen::Vector3d int_normal;
    int_normal[0] = RandomUniformReal(-0.3, 0.3);
    int_normal[1] = RandomUniformReal(-0.3, 0.3);
    int_normal[2] = RandomUniformReal(0.7, 1.3);

    int_normal.normalize();

    // int_normal = Eigen::Vector3d::UnitZ();

    std::vector<double> flatport_params = {int_normal[0],
                                           int_normal[1],
                                           int_normal[2],
                                           0.01,
                                           0.014,
                                           1.0,
                                           1.52,
                                           1.334};
    camera.SetRefracParams(flatport_params);

    for (int i = 0; i < 10; i++) {
      const double x = RandomUniformReal(
          0.5, 10.0 /*static_cast<double>(camera.Width()) - 0.5*/);
      const double y = RandomUniformReal(
          0.5, 10.0 /*static_cast<double>(camera.Height()) - 0.5*/);

      Eigen::Vector2d point2D(x, y);
      Ray3D ray_refrac = camera.CamFromImgRefrac(point2D);
      const Eigen::Vector3d virtual_cam_center =
          camera.VirtualCameraCenter(ray_refrac);

      std::cout << "point: " << point2D.transpose()
                << ", virtual cam center: " << virtual_cam_center.transpose()
                << std::endl;
    }
  }
  return true;
}