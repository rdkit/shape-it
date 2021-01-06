/*******************************************************************************
alignLib.cpp - Shape-it

Copyright 2021 by Greg Landrum

This file is part of Shape-it.

***********************************************************************/
#include <alignLib.h>

#include <alignmentInfo.h>
#include <bestResults.h>
#include <gaussianVolume.h>
#include <moleculeRotation.h>
#include <shapeAlignment.h>

SolutionInfo alignMols(const Molecule &refMol, const Molecule &dbMol,
                       const std::string &whichScore, double maxIter,
                       double cutoff, BestResults *bestHits) {
  // Create the refence set of Gaussians
  GaussianVolume refVolume;

  // List all Gaussians and their respective intersections
  listAtomVolumes(refMol, refVolume);

  // Move the Gaussian towards its center of geometry and align with principal
  // axes
  initOrientation(refVolume);

  // Create the set of Gaussians of database molecule
  GaussianVolume dbVolume;
  listAtomVolumes(dbMol, dbVolume);
  initOrientation(dbVolume);

  // Overlap with reference
  AlignmentInfo res;
  double bestScore(0.0);

  SolutionInfo bestSolution =
      std::move(alignVolumes(refVolume, dbVolume, whichScore, maxIter));
#ifndef USE_RDKIT
  bestSolution.dbMol = dbMol;
  bestSolution.refName = refMol.GetTitle();
  bestSolution.dbName = dbMol.GetTitle();
#else
  bestSolution.dbMol = dbMol;
  refMol.getProp("_Name", bestSolution.refName);
  dbMol.getProp("_Name", bestSolution.dbName);
#endif

  // Cleanup local pointers to atom-gaussians
  dbVolume.gaussians.clear();
  dbVolume.levels.clear();
  for (std::vector<std::vector<unsigned int> *>::iterator si =
           dbVolume.childOverlaps.begin();
       si != dbVolume.childOverlaps.end(); ++si) {
    if (*si != NULL) {
      delete *si;
      *si = NULL;
    }
  }
  dbVolume.childOverlaps.clear();
  refVolume.gaussians.clear();
  for (std::vector<std::vector<unsigned int> *>::iterator si =
           refVolume.childOverlaps.begin();
       si != refVolume.childOverlaps.end(); ++si) {
    if (*si != NULL) {
      delete *si;
      *si = NULL;
    }
  }
  refVolume.childOverlaps.clear();

  if (bestSolution.score < cutoff) {
    return bestSolution;
  }

  // Add the score properties
  setAllScores(bestSolution);

  // Translate and rotate the molecule towards its centroid and inertia
  // axes
  positionMolecule(bestSolution.dbMol, bestSolution.dbCenter,
                   bestSolution.dbRotation);

  // Rotate molecule with the optimal
  rotateMolecule(bestSolution.dbMol, bestSolution.rotor);

  // Rotate and translate the molecule with the inverse rotation and
  // translation of the reference molecule
  repositionMolecule(bestSolution.dbMol, refVolume.rotation,
                     refVolume.centroid);

  if (bestHits) {
    bestHits->add(bestSolution);
  }
  return bestSolution;
}

SolutionInfo alignVolumes(const GaussianVolume &refVolume,
                          const GaussianVolume &dbVolume,
                          const std::string &whichScore, double maxIter) {
  SolutionInfo res;
  res.refAtomVolume = refVolume.overlap;
  res.refCenter = refVolume.centroid;
  res.refRotation = refVolume.rotation;

  ShapeAlignment aligner(refVolume, dbVolume);
  aligner.setMaxIterations(maxIter);

  AlignmentInfo alignment;
  double bestScore = 0;
  for (unsigned int l = 0; l < 4; ++l) {
    SiMath::Vector quat(4, 0.0);
    quat[l] = 1.0;
    AlignmentInfo nextAlignment = aligner.gradientAscent(quat);
    checkVolumes(refVolume, dbVolume, nextAlignment);
    double ss = getScore(whichScore, nextAlignment.overlap, refVolume.overlap,
                         dbVolume.overlap);
    if (ss > bestScore) {
      alignment = nextAlignment;
      bestScore = ss;
    }

    if (bestScore > 0.98) {
      break;
    }
  }

  // Check if additional simulated annealing steps are requested and start
  // from the current best solution
  if (maxIter > 0) {
    AlignmentInfo nextRes = aligner.simulatedAnnealing(alignment.rotor);
    checkVolumes(refVolume, dbVolume, nextRes);
    double ss = getScore(whichScore, nextRes.overlap, refVolume.overlap,
                         dbVolume.overlap);
    if (ss > bestScore) {
      bestScore = ss;
      alignment = nextRes;
    }
  }
  // Optimal alignment information is stored in res and bestScore
  // => result reporting and post-processing
  updateSolutionInfo(res, alignment, bestScore, dbVolume);

  return std::move(res);
}
