/*******************************************************************************
alignLib.h - Shape-it

Copyright 2021 by Greg Landrum and the Shape-it contributors

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
#ifndef SHAPEIT_ALIGNLIB_H
#define SHAPEIT_ALIGNLIB_H

#include "bestResults.h"
#include "gaussianVolume.h"
#include "options.h"
#include "solutionInfo.h"
#include <string>

#ifndef USE_RDKIT
#include <openbabel/mol.h>
using Molecule = OpenBabel::OBMol;
#else
#include <GraphMol/ROMol.h>
using Molecule = RDKit::ROMol;
#endif

namespace shapeit {
SolutionInfo alignMols(const Molecule &refMol, const Molecule &dbMol,
                       const std::string &whichScore = tanimoto,
                       double maxIter = 0, double cutoff = 0.0,
                       BestResults *bestHits = nullptr);

SolutionInfo alignMolToVolume(const GaussianVolume &refVolume,
                              const Molecule &dbMol,
                              const std::string &whichScore = tanimoto,
                              double maxIter = 0, double cutoff = 0.0,
                              BestResults *bestHits = nullptr);

SolutionInfo alignVolumes(const GaussianVolume &refVolume,
                          const GaussianVolume &dbVolume,
                          const std::string &whichScore, double maxIter);
} // namespace shapeit
#endif