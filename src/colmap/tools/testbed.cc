#include "colmap/estimators/two_view_geometry.h"
#include "colmap/geometry/rigid3.h"
#include "colmap/math/random.h"
#include "colmap/scene/camera.h"

using namespace colmap;

int main(int argc, char* argv[]) {
  if (true) {
    SetPRNGSeed(time(NULL));
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

    // colmap::Camera camera;
    // camera.SetWidth(2048);
    // camera.SetHeight(1536);
    // camera.SetModelIdFromName("PINHOLE");
    // std::vector<double> params = {
    //     500.900000, 500.900000, 1024.000000, 768.000000};

    camera.SetParams(params);

    // Flatport setup.
    camera.SetRefracModelIdFromName("FLATPORT");
    Eigen::Vector3d int_normal;
    int_normal[0] = colmap::RandomUniformReal(-0.3, 0.3);
    int_normal[1] = colmap::RandomUniformReal(-0.3, 0.3);
    int_normal[2] = colmap::RandomUniformReal(0.7, 1.3);
    int_normal.normalize();

    // int_normal = Eigen::Vector3d::UnitZ();

    std::vector<double> flatport_params = {
        int_normal[0], int_normal[1], int_normal[2], 2.0, 1.0, 1.0, 5.3, 3.7};
    camera.SetRefracParams(flatport_params);

    const size_t kNumPoints = 200;

    // NOLINTNEXTLINE(clang-analyzer-security.FloatLoopCounter)
    for (double qx = 0; qx < 0.4; qx += 0.1) {
      // NOLINTNEXTLINE(clang-analyzer-security.FloatLoopCounter)
      for (double tx = 0; tx < 0.5; tx += 0.1) {
        Rigid3d cam2_from_cam1_gt(Eigen::Quaterniond(1, qx, 0, 0).normalized(),
                                  Eigen::Vector3d(tx, 0.1, 0));

        std::vector<Eigen::Vector2d> points2D1_refrac;
        std::vector<Eigen::Vector2d> points2D2_refrac;

        points2D1_refrac.reserve(kNumPoints);
        points2D2_refrac.reserve(kNumPoints);

        std::cout << "Generating data ... " << std::endl;
        size_t cnt = 0;
        while (cnt < kNumPoints) {
          // Randomly generate 3D points.
          const double u = RandomUniformReal(-2.0, 2.0);
          const double v = RandomUniformReal(-2.0, 2.0);
          const double w = RandomUniformReal(0.5, 5.5);

          Eigen::Vector3d point3D(u, v, w);
          // check for projection.
          Eigen::Vector2d point2D1 = camera.ImgFromCamRefrac(point3D);
          Eigen::Vector2d point2D2 =
              camera.ImgFromCamRefrac(cam2_from_cam1_gt * point3D);

          // Points must be visible by the camera.
          if (point2D1.x() < 0 || point2D1.x() > camera.Width() ||
              point2D1.y() < 0 || point2D1.y() > camera.Height()) {
            continue;
          }
          if (point2D2.x() < 0 || point2D2.x() > camera.Width() ||
              point2D2.y() < 0 || point2D2.y() > camera.Height()) {
            continue;
          }

          points2D1_refrac.push_back(point2D1);
          points2D2_refrac.push_back(point2D2);
          cnt++;
        }

        // Construct virtual cameras for two-view geometry estimation.
        std::vector<Camera> virtual_cameras1;
        Eigen::Quaterniond virtual_from_real_rotation1;
        std::vector<Rigid3d> virtual_from_reals1;

        std::vector<Camera> virtual_cameras2;
        Eigen::Quaterniond virtual_from_real_rotation2;
        std::vector<Rigid3d> virtual_from_reals2;

        virtual_cameras1.reserve(kNumPoints);
        virtual_from_reals1.reserve(kNumPoints);

        virtual_cameras2.reserve(kNumPoints);
        virtual_from_reals2.reserve(kNumPoints);

        virtual_from_real_rotation1 = camera.VirtualFromRealRotation();
        virtual_from_real_rotation2 = camera.VirtualFromRealRotation();

        for (const Eigen::Vector2d& point : points2D1_refrac) {
          const Ray3D ray_refrac = camera.CamFromImgRefrac(point);
          const Eigen::Vector3d virtual_cam_center =
              camera.VirtualCameraCenter(ray_refrac);
          virtual_from_reals1.push_back(
              Rigid3d(virtual_from_real_rotation1,
                      virtual_from_real_rotation1 * -virtual_cam_center));
          virtual_cameras1.push_back(
              camera.VirtualCamera(point, ray_refrac.dir.hnormalized()));
        }

        for (const Eigen::Vector2d& point : points2D2_refrac) {
          const Ray3D ray_refrac = camera.CamFromImgRefrac(point);
          const Eigen::Vector3d virtual_cam_center =
              camera.VirtualCameraCenter(ray_refrac);
          virtual_from_reals2.push_back(
              Rigid3d(virtual_from_real_rotation2,
                      virtual_from_real_rotation2 * -virtual_cam_center));
          virtual_cameras2.push_back(
              camera.VirtualCamera(point, ray_refrac.dir.hnormalized()));
        }

        FeatureMatches matches;
        matches.reserve(kNumPoints);

        for (point2D_t i = 0; i < kNumPoints; ++i) {
          matches.emplace_back(i, i);
        }
        std::cout << "Estimate ... " << std::endl;

        TwoViewGeometryOptions two_view_geometry_options;

        TwoViewGeometry two_view_geometry =
            EstimateRefractiveTwoViewGeometry(points2D1_refrac,
                                              virtual_cameras1,
                                              virtual_from_reals1,
                                              points2D2_refrac,
                                              virtual_cameras2,
                                              virtual_from_reals2,
                                              matches,
                                              two_view_geometry_options);

        std::cout << "Matrix norm: "
                  << (cam2_from_cam1_gt.ToMatrix() -
                      two_view_geometry.cam2_from_cam1.ToMatrix())
                         .norm()
                  << std::endl;
        std::cout << "Inlier ratio: "
                  << static_cast<double>(
                         two_view_geometry.inlier_matches.size()) /
                         static_cast<double>(kNumPoints)
                  << std::endl;
      }
    }
  }
  return true;
}