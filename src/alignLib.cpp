/*******************************************************************************
alignLib.cpp - Shape-it

Copyright 2021 by Greg Landrum

This file is part of Shape-it.

***********************************************************************/

#include <alignmentInfo.h>
#include <bestResults.h>
#include <gaussianVolume.h>
#include <shapeAlignment.h>

SolutionInfo alignVolumes(GaussianVolume &refVolume, GaussianVolume &dbVolume,
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
