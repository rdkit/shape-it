/*******************************************************************************
gaussianVolume.h - Shape-it

Copyright 2012 by Silicos-it, a division of Imacosi BVBA

This file is part of Shape-it.

        Shape-it is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published
        by the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        Shape-it is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with Shape-it.  If not, see <http://www.gnu.org/licenses/>.

Shape-it is linked against OpenBabel version 2.

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

void listAtomVolumes(Molecule &mol, GaussianVolume &gv);
void initOrientation(GaussianVolume &);
double atomOverlap(const GaussianVolume &, const GaussianVolume &);
double GAlpha(unsigned int);
double getScore(const std::string &, double, double, double);
void checkVolumes(const GaussianVolume &, const GaussianVolume &,
                  AlignmentInfo &);

#endif
