#pragma once

#include "colmap/sensor/models.h"
#include "colmap/sensor/ray3d.h"

namespace colmap {

// This file defines several different Non-Single-View-Point (non-SVP) camera
// models and arbitrary new non-SVP camera models can be added by the following
// steps:
//
//  1. Add a new struct in this file which implements all the necessary methods.
//  2. Define an unique non_svp_model_name and non_svp_model_id for the camera
//  model.
//  3. Add camera model to `CAMERA_NON_SVP_MODEL_CASES` macro in this file.
//  4. Add new template specialization of test case for camera model to
//     `camera_non_svp_models_test.cc`.
//
// A camera model can have three different types of camera parameters: focal
// length, principal point, extra parameters (distortion parameters). The
// parameter array is split into different groups, so that we can enable or
// disable the refinement of the individual groups during bundle adjustment. It
// is up to the camera model to access the parameters correctly (it is free to
// do so in an arbitrary manner) - the parameters are not accessed from outside.
//
// A camera model must have the following methods:
//
//  - `WorldToImage`: transform 3D world point in camera coordinates frame to
//  image
//    coordinates (the inverse of `ImageToWorld`). Assumes that the world
//    coordinates are given as (x, y, z).
//  - `ImageToWorld`: transform image coordinates to 3D ray in camera
//  coordinates frame
// (the inverse of `WorldToImage`). Produces the ray as ray origin and ray unit
// direction
//
//  - `ImageToWorldPoint`: transform image coordinates to 3D world point given
//  depth value
// of this pixel position. Produces the 3D world point in the camera coordinate
// frame (x, y, z).
//
// Whenever you specify the camera parameters in a list, they must appear
// exactly in the order as they are accessed in the defined model struct.
//
// The camera models follow the convention that the upper left image corner has
// the coordinate (0, 0), the lower right corner (width, height), i.e. that
// the upper left pixel center has coordinate (0.5, 0.5) and the lower right
// pixel center has the coordinate (width - 0.5, height - 0.5).

static const int kInvalidRefractiveCameraModelId = -1;

#ifndef CAMERA_REFRAC_MODEL_DEFINITIONS
#define CAMERA_REFRAC_MODEL_DEFINITIONS(                                       \
    refrac_model_id_value, refrac_model_name_value, num_params_value)          \
  static const int kRefracModelId = refrac_model_id_value;                     \
  static const size_t kNumParams = num_params_value;                           \
  static const int refrac_model_id;                                            \
  static const std::string refrac_model_name;                                  \
  static const size_t num_params;                                              \
  static const std::string params_info;                                        \
                                                                               \
  static inline int InitializeRefracModelId() {                                \
    return refrac_model_id_value;                                              \
  };                                                                           \
  static inline std::string InitializeRefracModelName() {                      \
    return refrac_model_name_value;                                            \
  }                                                                            \
  static inline size_t InitializeNumParams() { return num_params_value; }      \
  static inline std::string InitializeRefracModelParamsInfo();                 \
                                                                               \
  template <typename CameraModel, typename T>                                  \
  static void ImgFromCam(                                                      \
      const T* cam_params, const T* refrac_params, T u, T v, T w, T* x, T* y); \
                                                                               \
  template <typename CameraModel, typename T>                                  \
  static void CamFromImg(const T* cam_params,                                  \
                         const T* refrac_params,                               \
                         T x,                                                  \
                         T y,                                                  \
                         Eigen::Matrix<T, 3, 1>* ori,                          \
                         Eigen::Matrix<T, 3, 1>* dir);                         \
  template <typename CameraModel, typename T>                                  \
  static void CamFromImgPoint(const T* cam_params,                             \
                              const T* refrac_params,                          \
                              T x,                                             \
                              T y,                                             \
                              T d,                                             \
                              Eigen::Matrix<T, 3, 1>* uvw);
#endif

#ifndef CAMERA_REFRAC_MODEL_CASES
#define CAMERA_REFRAC_MODEL_CASES    \
  CAMERA_REFRAC_MODEL_CASE(FlatPort) \
  CAMERA_REFRAC_MODEL_CASE(DomePort)
#endif

#ifndef CAMERA_COMBINATION_MODEL_CASES
#define CAMERA_COMBINATION_MODEL_CASES                                    \
  CAMERA_COMBINATION_MODEL_CASE(FlatPort, SimplePinholeCameraModel)       \
  CAMERA_COMBINATION_MODEL_CASE(FlatPort, PinholeCameraModel)             \
  CAMERA_COMBINATION_MODEL_CASE(FlatPort, SimpleRadialCameraModel)        \
  CAMERA_COMBINATION_MODEL_CASE(FlatPort, SimpleRadialFisheyeCameraModel) \
  CAMERA_COMBINATION_MODEL_CASE(FlatPort, RadialCameraModel)              \
  CAMERA_COMBINATION_MODEL_CASE(FlatPort, RadialFisheyeCameraModel)       \
  CAMERA_COMBINATION_MODEL_CASE(FlatPort, OpenCVCameraModel)              \
  CAMERA_COMBINATION_MODEL_CASE(FlatPort, OpenCVFisheyeCameraModel)       \
  CAMERA_COMBINATION_MODEL_CASE(FlatPort, FullOpenCVCameraModel)          \
  CAMERA_COMBINATION_MODEL_CASE(FlatPort, FOVCameraModel)                 \
  CAMERA_COMBINATION_MODEL_CASE(FlatPort, ThinPrismFisheyeCameraModel)    \
  CAMERA_COMBINATION_MODEL_CASE(FlatPort, MetashapeFisheyeCameraModel)    \
  CAMERA_COMBINATION_MODEL_CASE(DomePort, SimplePinholeCameraModel)       \
  CAMERA_COMBINATION_MODEL_CASE(DomePort, PinholeCameraModel)             \
  CAMERA_COMBINATION_MODEL_CASE(DomePort, SimpleRadialCameraModel)        \
  CAMERA_COMBINATION_MODEL_CASE(DomePort, SimpleRadialFisheyeCameraModel) \
  CAMERA_COMBINATION_MODEL_CASE(DomePort, RadialCameraModel)              \
  CAMERA_COMBINATION_MODEL_CASE(DomePort, RadialFisheyeCameraModel)       \
  CAMERA_COMBINATION_MODEL_CASE(DomePort, OpenCVCameraModel)              \
  CAMERA_COMBINATION_MODEL_CASE(DomePort, OpenCVFisheyeCameraModel)       \
  CAMERA_COMBINATION_MODEL_CASE(DomePort, FullOpenCVCameraModel)          \
  CAMERA_COMBINATION_MODEL_CASE(DomePort, FOVCameraModel)                 \
  CAMERA_COMBINATION_MODEL_CASE(DomePort, ThinPrismFisheyeCameraModel)    \
  CAMERA_COMBINATION_MODEL_CASE(DomePort, MetashapeFisheyeCameraModel)
#endif

#ifndef CAMERA_REFRAC_MODEL_SWITCH_CASES
#define CAMERA_REFRAC_MODEL_SWITCH_CASES         \
  CAMERA_REFRAC_MODEL_CASES                      \
  default:                                       \
    CAMERA_REFRAC_MODEL_DOES_NOT_EXIST_EXCEPTION \
    break;
#endif

#define CAMERA_REFRAC_MODEL_DOES_NOT_EXIST_EXCEPTION \
  throw std::domain_error("Refractive camera model does not exist");

#define CAMERA_COMBINATION_MODEL_DOES_NOT_EXIST_EXCEPTION                     \
  throw std::domain_error(                                                    \
      "The combination of a perspective camera model and a refractive model " \
      "does not exist");

#ifndef CAMERA_COMBINATION_MODEL_IF_ELSE_CASES
#define CAMERA_COMBINATION_MODEL_IF_ELSE_CASES        \
  CAMERA_COMBINATION_MODEL_CASES                      \
                                                      \
  {                                                   \
    CAMERA_COMBINATION_MODEL_DOES_NOT_EXIST_EXCEPTION \
  }
#endif

// The "Curiously Recurring Template Pattern" (CRTP) is
// used here, so that we can reuse some shared
// functionality between all camera models - defined in
// the BaseCameraModel.
template <typename CameraRefracModel>
struct BaseCameraRefracModel {
  template <typename CameraModel, typename T>
  static inline void IterativeProjection(
      const T* cam_params, const T* refrac_params, T u, T v, T w, T* x, T* y);
};

// DoubleLayerPlanarRefractiveInterface (thick flat
// port camera system).
//
// Parameter list is expected in the following order:
//
// Nx, Ny, Nz, int_dist, int_thick, na, ng, nw (Nx, Ny,
// Nz is the planar interface normal in the local
// camera coordinate system)
//
// Note that (Nx, Ny, Nz) must be unit vector, e.g. (0,
// 0, 1) means the plane normal is the same as the
// camera's local Z-axis
//
// see:
// https://link.springer.com/chapter/10.1007/978-3-642-33715-4_61
struct FlatPort : public BaseCameraRefracModel<FlatPort> {
  CAMERA_REFRAC_MODEL_DEFINITIONS(0, "FLATPORT", 8)
};

// DoubleLayerSphericalRefractiveInterface (thick dome
// port camera system).
//
// Parameter list is expected in the following order:
//
// Cx, Cy, Cz, int_radius, int_thick, na, ng, nw (Cx,
// Cy, Cz is the spherical center in the local camera
// coordinate system)
//
// see:
// https://link.springer.com/chapter/10.1007/978-3-030-33676-9_6
struct DomePort : public BaseCameraRefracModel<DomePort> {
  CAMERA_REFRAC_MODEL_DEFINITIONS(1, "FLATPORT", 8)
};

// Check whether refractive camera with given name or identifier
// exists
bool ExistsCameraRefracModelWithName(const std::string& refrac_model_name);
bool ExistsCameraRefracModelWithId(const int refrac_model_id);

// Covnert camera non-single-view-point model to unique model identifier
//
// @param name      Unique name of non-SVP camera model
//
// @return          Unique identifier of non-SVP camera model
int CameraRefracModelNameToId(const std::string& refrac_model_name);

// Convert camera non-single-view-point model unique identifier to name
//
// @param model_id  Unique identifier of non-SVP camera model
//
// @param           Unique name of non-SVP camera model
std::string CameraRefracModelIdToName(const int refrac_model_id);

// Get human-readable information about the non-SVP parameter vector order
//
// @param model_id   Unique identifier of non-SVP camera model
std::string CameraRefracModelParamsInfo(const int refrac_model_id);

// Get the total number of parameters of a non-SVP camera model
//
// @param       Unique identifier of non-SVP camera model
size_t CameraRefracModelNumParams(const int refrac_model_id);

// Check whether parameters are valid, i.e. the parameter vector has
// the correct dimensions that match the specified non-SVP camera model.
//
// @param model_id      Unique identifier of non-SVP camera model.
// @param params        Array of camera parameters.
bool CameraRefracModelVerifyParams(const int refrac_model_id,
                                   const std::vector<double>& params);

// Transform world coordinates in camera coordinate system to image coordinates.
//
// This is the inverse of `CameraNonSvpModelImageToWorld`.
//
// @param model_id              Unique model_id of camera model as defined in
//                              `CAMERA_MODEL_NAME_TO_CODE`.
// @param non_svp_model_id      Unique identifier of non-SVP camera model
// @param camera_params         Array of camera parameters.
// @param non_svp_params        Array of non-SVP model parameters.
// @param world_point           Coordinates in camera system as (x, y, z).
// @param x, y                  Output image coordinates in pixels.
inline Eigen::Vector2d CameraRefracModelImgFromCam(
    const int model_id,
    const int refrac_model_id,
    const std::vector<double>& cam_params,
    const std::vector<double>& refrac_params,
    const Eigen::Vector3d& uvw);

// Transform image to 3D ray in camera frame using refractive model.
//
// This is the inverse of `CameraRefracModelImgFromCam`.
//
// @param model_id              Unique identifier of camera model.
// @param non_svp_model_id      Unique identifier of non-SVP camera model
// @param camera_params         Array of camera parameters.
// @param non_svp_params        Array of non-SVP model parameters.
// @param x, y                  Image coordinates in pixels.
// @param ray_ori, ray_dir      Output ray in camera system as origin and unit
// direction.
inline Ray3D CameraRefracModelCamFromImg(
    const int model_id,
    const int refrac_model_id,
    const std::vector<double>& cam_params,
    const std::vector<double>& refrac_params,
    const Eigen::Vector2d& xy);

////////////////////////////////////////////////////////////////////////////////
// BaseCameraRefracModel

template <typename CameraRefracModel>
template <typename CameraModel, typename T>
void BaseCameraRefracModel<CameraRefracModel>::IterativeProjection(
    const T* cam_params, const T* refrac_params, T u, T v, T w, T* x, T* y) {
  // Parameters for Newton iteration using numerical differentiation with
  // central differences, 100 iterations should be enough even for complex
  // camera models with higher order terms.
  const Eigen::Matrix<T, 3, 1> uvw(u, v, w);
  const T d = uvw.norm();
  const size_t kNumIterations = 100;
  const T kMaxStepNorm = T(1e-10);
  const T kRelStepSize = T(1e-7);
  const T kAbsStepSize = T(1e-6);

  Eigen::Matrix<T, 3, 2> J;
  const Eigen::Matrix<T, 2, 1> X0(*x, *y);
  Eigen::Matrix<T, 2, 1> X(*x, *y);
  Eigen::Matrix<T, 3, 1> err;
  Eigen::Matrix<T, 3, 1> dx_0b;
  Eigen::Matrix<T, 3, 1> dx_0f;
  Eigen::Matrix<T, 3, 1> dx_1b;
  Eigen::Matrix<T, 3, 1> dx_1f;

  for (size_t i = 0; i < kNumIterations; ++i) {
    const T step0 = std::max(kAbsStepSize, ceres::abs(kRelStepSize * X(0)));
    const T step1 = std::max(kAbsStepSize, ceres::abs(kRelStepSize * X(1)));
    CameraRefracModel::template CamFromImgPoint<CameraModel, T>(
        cam_params, refrac_params, X(0), X(1), d, &err);
    err = uvw - err;
    CameraRefracModel::template CamFromImgPoint<CameraModel, T>(
        cam_params, refrac_params, X(0) - step0, X(1), d, &dx_0b);

    CameraRefracModel::template CamFromImgPoint<CameraModel, T>(
        cam_params, refrac_params, X(0) + step0, X(1), d, &dx_0f);

    CameraRefracModel::template CamFromImgPoint<CameraModel, T>(
        cam_params, refrac_params, X(0), X(1) - step1, d, &dx_1b);

    CameraRefracModel::template CamFromImgPoint<CameraModel, T>(
        cam_params, refrac_params, X(0), X(1) + step1, d, &dx_1f);

    J.col(0) = (dx_0b - dx_0f) / (T(2.0) * step0);
    J.col(1) = (dx_1b - dx_1f) / (T(2.0) * step1);
    Eigen::Matrix<T, 2, 2> H = J.transpose() * J;
    Eigen::Matrix<T, 2, 1> b = T(-1.0) * J.transpose() * err;
    const Eigen::Matrix<T, 2, 1> step_x = H.ldlt().solve(b);
    X += step_x;
    if (step_x.squaredNorm() < kMaxStepNorm) break;
  }
  *x = X(0);
  *y = X(1);
}

}  // namespace colmap