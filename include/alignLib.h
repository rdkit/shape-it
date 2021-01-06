/*******************************************************************************
alignLib.h - Shape-it

Copyright 2021 by Greg Landrum

This file is part of Shape-it.

***********************************************************************/
#ifndef SHAPEIT_ALIGNLIB_H
#define SHAPEIT_ALIGNLIB_H

#include <string>

class GaussianVolume;
class SolutionInfo;
SolutionInfo alignVolumes(GaussianVolume &refVolume, GaussianVolume &dbVolume,
                          const std::string &whichScore, double maxIter);
#endif