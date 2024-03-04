COLMAP Underwater
======

About
-----

This repo is a fork of the original Structure-from-Motion (SfM) framework COLMAP (https://github.com/colmap/colmap) which mainly focuses on 3D reconstruction in the underwater domain.
It supports Refractive Structure-from-Motion (RSfM) when using cameras underwater with waterproof housings (flat-ports / dome-ports), and it supports using navigation data in the reconstruction to help reducing drift in the large-scale robotic mapping scenarios.

Features include:
- Implementation of the commonly used refractive camera models (flat-port / dome-port camera models). 
- Refractive Structure-from-Motion.
- Bundle adjustment with pose priors (e.g. navigation data) in SfM to reduce potential drift in the large-scale robotic visual mapping scenario.

Since the main framework is created by the original COLMAP, please also cite the following papers if you use this project for your research:

    @inproceedings{schoenberger2016sfm,
        author={Sch\"{o}nberger, Johannes Lutz and Frahm, Jan-Michael},
        title={Structure-from-Motion Revisited},
        booktitle={Conference on Computer Vision and Pattern Recognition (CVPR)},
        year={2016},
    }

    @inproceedings{schoenberger2016mvs,
        author={Sch\"{o}nberger, Johannes Lutz and Zheng, Enliang and Pollefeys, Marc and Frahm, Jan-Michael},
        title={Pixelwise View Selection for Unstructured Multi-View Stereo},
        booktitle={European Conference on Computer Vision (ECCV)},
        year={2016},
    }

If you use the image retrieval / vocabulary tree engine, please also cite:

    @inproceedings{schoenberger2016vote,
        author={Sch\"{o}nberger, Johannes Lutz and Price, True and Sattler, Torsten and Frahm, Jan-Michael and Pollefeys, Marc},
        title={A Vote-and-Verify Strategy for Fast Spatial Verification in Image Retrieval},
        booktitle={Asian Conference on Computer Vision (ACCV)},
        year={2016},
    }



Compilation
--------

Building COLMAP Underwater from source is exactly the same as building the original COLMAP as there is no extra dependency. Therefore, please see  https://colmap.github.io/install.html

Why COLMAP Underwater?
---------------

To protect cameras from water and pressure, they are enclosed in waterproof housings, and observe the environment through a transparent window, typically with a *planar* or *spherical* shape. Light rays from the underwater scene change direction when they travel through these interfaces in a non-orthogonal manner, leading to distortion in the acquired images.

Example underwater cameras are:

<p align="center">
  <img src="doc/images/vive_dome.jpg" height="320" />
  <img src="doc/images/flatport_with_lights.jpg" height="320" /> 
  <img src="doc/images/flatport_deepsea.jpg" height="320" />
</p>

The planar shaped housing interfaces are referred to as **flat-ports**, and  the spherical shaped ones are referred to as **dome-ports**.

The housing interfaces are modeled with additional parameters in the image formation process to account for the refraction effects.

- The Flat-port parameters:
  -  the unit-length interface normal vector $\mathbf{n}_{\mathrm{int}} = (n_x, n_y, n_z)^T$. The normal points towards the positive $Z$-axis, with $\mathbf{n}_{\mathrm{int}} = (0, 0, 1)^T$ coinciding with the optical axis of the camera.
  - the camera-to-interface distance $d_{\mathrm{int}}$ (the orthogonal distance from the camera projection center to the interface plane). (unit: [$m$])
  - the thickness of the interface (unit: [$m$])
  - refraction indices of air, glass and water: $n_a, n_g, n_w$. For example: $n_a = 1.0$, $n_g = 1.49$, $n_w = 1.334$ 
  
- The Dome-port parameters:
  -  the dome center in the local camera coordinate frame $\mathbf{C}_d = (c_x, c_y, c_z)^T$ (unit: [$m$]). If $\mathbf{C}_d = (0, 0, 0)^T$, then the dome-port is perfectly centered with the camera, refraction will not occur at the interface.
  - the dome-port radius and thickness. (unit: [$m$])
  - refraction indices of air, glass and water: $n_a, n_g, n_w$.
  
Many traditional multi-view geometry techniques, e.g. two-view geometry, pose estimation, bundle adjustment are developed based on the pinhole camera model. We therefore need to adapt the framework such that these refractive camera models can be utilized for Structure-from-Motion.

What's Different?
---------------

Comming soon...

Usage tips
---------------

Comming soon...

License
-------

The COLMAP library is licensed under the new BSD license. Note that this text
refers only to the license for COLMAP itself, independent of its thirdparty
dependencies, which are separately licensed. Building COLMAP with these
dependencies may affect the resulting COLMAP license.

    Copyright (c) 2023, ETH Zurich and UNC Chapel Hill.
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

        * Redistributions of source code must retain the above copyright
          notice, this list of conditions and the following disclaimer.

        * Redistributions in binary form must reproduce the above copyright
          notice, this list of conditions and the following disclaimer in the
          documentation and/or other materials provided with the distribution.

        * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
          its contributors may be used to endorse or promote products derived
          from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.

