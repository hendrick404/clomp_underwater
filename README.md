COLMAP Underwater
======

About
-----

COLMAP underwater is a fork of the original Structure-from-Motion (SfM) framework COLMAP which mainly focuses on 3D reconstruction in the underwater domain.
It supports Refractive Structure-from-Motion (RSfM) when using cameras underwater with waterproof housings (flat-ports / dome-ports).

Features include:
- Implementation of the commonly used refractive camera models (flat-port / dome-port camera models). 
- Refractive Structure-from-Motion pipeline.
- Add pose priors (e.g. navigation data) as a soft constraint to the reconstruction, the resulting model is in the same scale as the pose priors.

Since the main incremental SfM algorithm remains unchanged as the original COLMAP, please also cite the following papers if you use this project for your research:

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

The latest source code is available at https://github.com/colmap/colmap. COLMAP
builds on top of existing works and when using specific algorithms within
COLMAP, please also cite the original authors, as specified in the source code.


Download
--------

Executables for Windows and Mac and other resources can be downloaded from
https://demuc.de/colmap/. Executables for Linux/Unix/BSD are available at
https://repology.org/metapackage/colmap/versions. To build COLMAP from source,
please see https://colmap.github.io/install.html.

Getting Started
---------------

1. Download the pre-built binaries from https://demuc.de/colmap/ or build the
   library manually as described in the documentation.
2. Download one of the provided datasets at https://demuc.de/colmap/datasets/
   or use your own images.
3. Use the **automatic reconstruction** to easily build models
   with a single click or command.


Documentation
-------------

The documentation is available at https://colmap.github.io/.


Support
-------

Please, use GitHub Discussions at https://github.com/colmap/colmap/discussions
for questions and the GitHub issue tracker at https://github.com/colmap/colmap
for bug reports, feature requests/additions, etc.


Acknowledgments
---------------

The library was originally written by Johannes L. Schönberger
(https://demuc.de/) with funding provided by his PhD advisors Jan-Michael Frahm
and Marc Pollefeys. Since then the project has benefited from countless
community contributions, including bug fixes, improvements, new features,
third-party tooling, and community support.


Contribution
------------

Contributions (bug reports, bug fixes, improvements, etc.) are very welcome and
should be submitted in the form of new issues and/or pull requests on GitHub.


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

