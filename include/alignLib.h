/*******************************************************************************
alignLib.h - Shape-it

Copyright 2021 by Greg Landrum

This file is part of Shape-it.

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

SolutionInfo alignMols(const Molecule &refMol, const Molecule &dbMol,
                       const std::string &whichScore = tanimoto,
                       double maxIter = 0, double cutoff = 0.0,
                       BestResults *bestHits = nullptr);

SolutionInfo alignVolumes(const GaussianVolume &refVolume,
                          const GaussianVolume &dbVolume,
                          const std::string &whichScore, double maxIter);
#endif