#!/bin/bash

WORKSPACE_PATH=$HOME/colmap_active_workspace
BACKUP_PATH=$HOME/colmap_workspace

CAM_PARAM_F=340.5
CAM_PARAM_CX=556.5
CAM_PARAM_CY=417.5
CAM_PARAM_N="0.1658429,0.14809199,0.97496918"
CAM_PARAM_INT_DIST=5
CAM_PARAM_INT_THICK=2
CAM_PARAM_NA=1
CAM_PARAM_NG=1.52
CAM_PARAM_NW=1.334

rm -rf $WORKSPACE_PATH
cp -r $BACKUP_PATH $WORKSPACE_PATH

colmap_underwater feature_extractor \
   --database_path $WORKSPACE_PATH/database.db \
   --image_path $WORKSPACE_PATH/images \
   --ImageReader.camera_refrac_model FLATPORT \
   --ImageReader.camera_params "$CAM_PARAM_F,$CAM_PARAM_CX,$CAM_PARAM_CY" \
   --ImageReader.camera_refrac_params "$CAM_PARAM_N,$CAM_PARAM_INT_DIST,$CAM_PARAM_INT_THICK,$CAM_PARAM_NA,$CAM_PARAM_NG,$CAM_PARAM_NW" \
   --ImageReader.single_camera 1 \
   --ImageReader.camera_model SIMPLE_PINHOLE 


colmap_underwater exhaustive_matcher \
   --database_path $WORKSPACE_PATH/database.db

mkdir $WORKSPACE_PATH/sparse

colmap_underwater mapper \
    --database_path $WORKSPACE_PATH/database.db \
    --image_path $WORKSPACE_PATH/images \
    --output_path $WORKSPACE_PATH/sparse

mkdir $WORKSPACE_PATH/dense

colmap_underwater image_undistorter \
    --image_path $WORKSPACE_PATH/images \
    --input_path $WORKSPACE_PATH/sparse/0 \
    --output_path $WORKSPACE_PATH/dense \
    --output_type COLMAP \
    --max_image_size 2000

if [ -x "$(command -v nvidia-smi)" ] ; then
    colmap patch_match_stereo \
        --workspace_path $WORKSPACE_PATH/dense \
        --workspace_format COLMAP \
        --PatchMatchStereo.geom_consistency true

    colmap stereo_fusion \
        --workspace_path $WORKSPACE_PATH/dense \
        --workspace_format COLMAP \
        --input_type geometric \
        --output_path $WORKSPACE_PATH/dense/fused.ply

    colmap poisson_mesher \
        --input_path $WORKSPACE_PATH/dense/fused.ply \
        --output_path $WORKSPACE_PATH/dense/meshed-poisson.ply

    colmap delaunay_mesher \
        --input_path $WORKSPACE_PATH/dense \
        --output_path $WORKSPACE_PATH/dense/meshed-delaunay.ply
fi
