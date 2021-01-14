/*******************************************************************************
solutionInfo.cpp - Shape-it

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

#include <solutionInfo.h>

#ifndef USE_RDKIT
// OpenBabel
#include <openbabel/generic.h>
#endif

SolutionInfo::SolutionInfo(void)
    : refName(""), refAtomVolume(0.0), refCenter(0, 0, 0), refRotation(3, 3, 0),
      dbName(""), dbAtomVolume(0.0), dbMol(), dbCenter(0, 0, 0),
      dbRotation(3, 3, 0), atomOverlap(0.0), score(0.0), rotor(4, 0.0) {
  rotor[0] = 1.0;
}

SolutionInfo::~SolutionInfo(void) {}

void SolutionInfo::printScores(Options &uo) {
  *(uo.scoreOutStream)
      << dbName << "\t" << refName << "\t" << std::setprecision(3)
      << atomOverlap / (refAtomVolume + dbAtomVolume - atomOverlap) << "\t"
      << std::setprecision(3)
      << atomOverlap / (0.95 * refAtomVolume + 0.05 * dbAtomVolume) << "\t"
      << std::setprecision(3)
      << atomOverlap / (0.05 * refAtomVolume + 0.95 * dbAtomVolume) << "\t"
      << std::setprecision(5) << atomOverlap << "\t" << std::setprecision(5)
      << refAtomVolume << "\t" << std::setprecision(5) << dbAtomVolume
      << std::endl;
  return;
}

void updateSolutionInfo(SolutionInfo &s, const AlignmentInfo &res, double score,
                        const GaussianVolume &gv) {
  s.dbAtomVolume = gv.overlap;
  s.dbCenter = gv.centroid;
  s.dbRotation = gv.rotation;
  s.atomOverlap = res.overlap;
  s.score = score;
  s.rotor = res.rotor;
  return;
}

void setAllScores(SolutionInfo &res) {
  std::ostringstream ss;

  ss.str("");
  ss << res.atomOverlap /
            (res.refAtomVolume + res.dbAtomVolume - res.atomOverlap);
#ifndef USE_RDKIT
  OpenBabel::OBPairData *label1 = new OpenBabel::OBPairData();
  label1->SetAttribute(tanimoto);
  label1->SetValue(ss.str());
  res.dbMol.SetData(label1);
#else
  res.dbMol.setProp(tanimoto, ss.str());
#endif
  ss.str("");
  ss << res.atomOverlap / (0.95 * res.refAtomVolume + 0.05 * res.dbAtomVolume);
#ifndef USE_RDKIT
  OpenBabel::OBPairData *label2 = new OpenBabel::OBPairData();
  label2->SetAttribute(tversky_ref);
  label2->SetValue(ss.str());
  res.dbMol.SetData(label2);
#else
  res.dbMol.setProp(tversky_ref, ss.str());
#endif

  ss.str("");
  ss << res.atomOverlap / (0.05 * res.refAtomVolume + 0.95 * res.dbAtomVolume);
#ifndef USE_RDKIT
  OpenBabel::OBPairData *label3 = new OpenBabel::OBPairData();
  label3->SetAttribute(tversky_db);
  label3->SetValue(ss.str());
  res.dbMol.SetData(label3);
#else
  res.dbMol.setProp(tversky_db, ss.str());
#endif

  return;
}
