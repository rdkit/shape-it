/*******************************************************************************
gaussianVolume.h - Shape-it

Copyright 2012-2021 by Silicos-it, a division of Imacosi BVBA, Hans De Winter,
and the Shape-it contributors

This file is part of Shape-it.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Shape-it can be linked against either OpenBabel version 3 or the RDKit.

        OpenBabel is free software; you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation version 2 of the License.

***********************************************************************/

#ifndef __SILICOSIT_SHAPEIT_GAUSSIANVOLUME_H__
#define __SILICOSIT_SHAPEIT_GAUSSIANVOLUME_H__

// General
#include <algorithm>
#include <cmath>
#include <queue>
#include <set>
#include <vector>

// Shape-it
#include <alignmentInfo.h>
#include <atomGaussian.h>
#include <coordinate.h>
#include <options.h>
#include <siMath.h>

const double GCI = 2.828427125;
const double GLI = 1.480960979;
const double VCUTOFF = 2.0;
const unsigned int LEVEL = 6;
const double EPS = 0.03;
const double GRADSCALE = 0.9;
const double PENALTY = 5.00;

#ifndef USE_RDKIT
namespace OpenBabel {
class OBMol;
}
using Molecule = OpenBabel::OBMol;
#else
namespace RDKit {
class ROMol;
}
using Molecule = RDKit::ROMol;
#endif

class GaussianVolume {
public:
  double volume;       ///< Molecular volume
  double overlap;      ///< Self-overlap of the molecule
  Coordinate centroid; ///< center of the gaussian volume
  SiMath::Matrix
      rotation; ///< rotation matrix to align molecule to principal axes
  std::vector<AtomGaussian>
      gaussians; ///< vector of all atom gaussians and their overlaps
  std::vector<std::vector<unsigned int> *>
      childOverlaps; ///< vector to keep track of which overlaps are formed with
                     ///< one gaussian
  std::vector<unsigned int>
      levels; ///< indicates where in the vector the level of overlaps changes

  GaussianVolume(void);
  ~GaussianVolume(void);
};

void listAtomVolumes(const Molecule &mol, GaussianVolume &gv);
void initOrientation(GaussianVolume &);
double atomOverlap(const GaussianVolume &, const GaussianVolume &);
double GAlpha(unsigned int);
double getScore(const std::string &, double, double, double);
void checkVolumes(const GaussianVolume &, const GaussianVolume &,
                  AlignmentInfo &);

#endif
