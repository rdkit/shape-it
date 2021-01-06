/*******************************************************************************
solutionInfo.h - Shape-it

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

#ifndef __SILICOSIT_SHAPEIT_SOLUTIONINFO_H__
#define __SILICOSIT_SHAPEIT_SOLUTIONINFO_H__

// General
#include <iomanip>
#include <string>

#ifndef USE_RDKIT
// OpenBabel
#include <openbabel/mol.h>
#else
#include <GraphMol/RWMol.h>
#endif

// Shape-it
#include <alignmentInfo.h>
#include <coordinate.h>
#include <gaussianVolume.h>
#include <options.h>
#include <siMath.h>

class SolutionInfo {
public:
  std::string refName;
  double refAtomVolume;
  Coordinate refCenter;
  SiMath::Matrix refRotation;

#ifndef USE_RDKIT
  OpenBabel::OBMol dbMol;
#else
  RDKit::RWMol dbMol;
#endif
  std::string dbName;
  double dbAtomVolume;
  Coordinate dbCenter;
  SiMath::Matrix dbRotation;

  double atomOverlap;
  double score;
  SiMath::Vector rotor;

  SolutionInfo(void);
  ~SolutionInfo(void);

  void printScores(Options &);
};

void setAllScores(SolutionInfo &);
void updateSolutionInfo(SolutionInfo &, AlignmentInfo &, double,
                        GaussianVolume &);

#endif
