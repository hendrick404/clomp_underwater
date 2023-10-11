#include "colmap/geometry/rigid3.h"
#include "colmap/math/random.h"
#include "colmap/scene/camera.h"

using namespace colmap;

int main(int argc, char* argv[]) {
  if (true) {
    SetPRNGSeed(time(NULL));

    std::cout << "Setup some realistic camera model" << std::endl;

    // Camera parameters coming from Anton 131 map.
    colmap::Camera camera;
    camera.SetWidth(4104);
    camera.SetHeight(3006);
    camera.SetModelIdFromName("METASHAPE_FISHEYE");
    std::vector<double> params = {2115.7964068878537,
                                  2115.6679248258547,
                                  2097.4832274550317,
                                  1641.1207545881762,
                                  -0.0047743457655361719,
                                  0.053835104379748221,
                                  -0.086202869606942054,
                                  0.041314184080814435,
                                  0.00012198372904429465,
                                  0.00051849746014895923};
    camera.SetParams(params);

    // Flatport setup.
    camera.SetRefracModelIdFromName("FLATPORT");
    Eigen::Vector3d int_normal;
    // int_normal[0] = colmap::RandomUniformReal(-0.3, 0.3);
    // int_normal[1] = colmap::RandomUniformReal(-0.3, 0.3);
    // int_normal[2] = colmap::RandomUniformReal(0.7, 1.3);
    // int_normal.normalize();

    int_normal = Eigen::Vector3d::UnitZ();

    std::vector<double> flatport_params = {int_normal[0],
                                           int_normal[1],
                                           int_normal[2],
                                           0.01,
                                           0.007,
                                           1.0,
                                           1.47,
                                           1.334};
    camera.SetRefracParams(flatport_params);

    const size_t kNumPoints = 200;
    
    for (size_t i = 0; i < num_points; i++) {
      Eigen::Vector2d point2D_refrac;
      point2D_refrac.x() = colmap::RandomUniformReal(
          0.5, static_cast<double>(camera.Width()) - 0.5);
      point2D_refrac.y() = colmap::RandomUniformReal(
          0.5, static_cast<double>(camera.Height()) - 0.5);

      colmap::Ray3D ray_refrac = camera.CamFromImgRefrac(point2D_refrac);
      const double d = colmap::RandomUniformReal(0.5, 10.0);

      Eigen::Vector3d point3D = ray_refrac.At(d);

      Eigen::Quaterniond virtual_from_real_rotation =
          camera.VirtualFromRealRotation();
      Eigen::Vector3d virtual_cam_center =
          camera.VirtualCameraCenter(point2D_refrac);

      Eigen::Vector3d virtual_from_real_translation =
          virtual_from_real_rotation * -virtual_cam_center;

      Rigid3d virtual_from_real(virtual_from_real_rotation,
                                virtual_from_real_translation);

      Eigen::Vector3d point3D_v = virtual_from_real * point3D;

      Camera virtual_camera =
          camera.VirtualCamera(point2D_refrac, point3D_v.hnormalized());

      std::cout
          << "image point: " << point2D_refrac.transpose()
          << " -- virtual point2D: "
          << virtual_camera.ImgFromCam(point3D_v.hnormalized()).transpose()
          << std::endl;
    }
  }
  return true;
}