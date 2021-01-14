/*******************************************************************************
main.cpp - Shape-it

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

// General
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <vector>

#ifndef USE_RDKIT
// OpenBabel
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#else
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/ROMol.h>
#endif

// Pharao
#include <alignLib.h>
#include <alignmentInfo.h>
#include <atomGaussian.h>
#include <bestResults.h>
#include <gaussianVolume.h>
#include <mainErr.h>
#include <moleculeRotation.h>
#include <options.h>
#include <parseCommandLine.h>
#include <printHeader.h>
#include <printUsage.h>
#include <shapeAlignment.h>

//*--------------------------------------------------------------------------*//
//* MAIN                                                                MAIN
//*//
//*--------------------------------------------------------------------------*//
int main(int argc, char *argv[]) {
  // Initialise random number generator
  srandom(time(NULL));
  clock_t t0 = clock();

  // Print header
  printHeader();

  // Read options
  Options uo = parseCommandLine(argc, argv);
  if (uo.version) {
    printHeader();
    exit(0);
  }

  if (uo.help) {
    printUsage();
    exit(0);
  }
  std::cerr << uo.print();

  // Files
  if (uo.dbInpFile.empty()) {
    mainErr("Missing database file. This is a required option (-d).");
  }
  if (uo.refInpFile.empty()) {
    mainErr("Missing ref file. This is a required option (-r).");
  }

  if (uo.molOutFile.empty() && uo.scoreOutFile.empty()) {
    mainErr("At least one of the -o or -s option should be used.");
  }

  // Create a list to store the best results
  BestResults *bestHits = NULL;
  if (uo.bestHits != 0) {
    bestHits = new BestResults(uo.bestHits);
  }

  // Print header line to score output file
  if (!uo.scoreOutFile.empty()) {
    *(uo.scoreOutStream) << "dbName"
                         << "\t"
                         << "refName"
                         << "\t" << tanimoto << "\t" << tversky_ref << "\t"
                         << tversky_db << "\t"
                         << "overlap"
                         << "\t"
                         << "refVolume"
                         << "\t"
                         << "dbVolume" << std::endl;
  }

  // Create reference molecule
  std::string refName;
#ifndef USE_RDKIT
  OpenBabel::OBMol refMol;
  OpenBabel::OBConversion refInpReader;
  if (uo.format.empty()) {
    refInpReader.SetInFormat(refInpReader.FormatFromExt(uo.refInpFile));
  } else {
    refInpReader.SetInFormat(refInpReader.FindFormat(uo.format));
  }
  refInpReader.Read(&refMol, uo.refInpStream);
  refName = refMol.GetTitle();
#else
  bool takeOwnership = false;
  bool sanitize = true;

  bool removeHs = false;
  RDKit::ForwardSDMolSupplier refsuppl(uo.refInpStream, takeOwnership, sanitize,
                                       removeHs);
  std::unique_ptr<RDKit::ROMol> refmptr(refsuppl.next());
  if (!refmptr) {
    mainErr("Could not parse reference molecule");
  }
  RDKit::ROMol &refMol = *refmptr;
  refMol.getPropIfPresent("_Name", refName);
#endif
  if (refName == "") {
    refName = "Unnamed_ref";
  }

  // Create the refence set of Gaussians
  GaussianVolume refVolume;

  // List all Gaussians and their respective intersections
  listAtomVolumes(refMol, refVolume);

  // Move the Gaussian towards its center of geometry and align with principal
  // axes
  if (!uo.scoreOnly) {
    initOrientation(refVolume);
  }

  // Write reference molecule to output
#ifndef USE_RDKIT
  OpenBabel::OBConversion dbOutWriter;
  if (uo.format.empty()) {
    dbOutWriter.SetOutFormat(dbOutWriter.FormatFromExt(uo.molOutFile));
  } else {
    dbOutWriter.SetOutFormat(dbOutWriter.FindFormat(uo.format));
  }
  if (uo.showRef && !uo.molOutFile.empty()) {
    dbOutWriter.Write(&refMol, uo.molOutStream);
  }
#else
  RDKit::SDWriter dbOutWriter(uo.molOutStream, false);
  dbOutWriter.write(refMol);
#endif

  // Open database stream
  unsigned molCount(0);
  std::ostringstream ss;
#ifndef USE_RDKIT
  OpenBabel::OBMol dbMol;
  OpenBabel::OBConversion dbInpReader;
  if (uo.format.empty()) {
    dbInpReader.SetInFormat(dbInpReader.FormatFromExt(uo.dbInpFile));
  } else {
    dbInpReader.SetInFormat(dbInpReader.FindFormat(uo.format));
  }
  dbInpReader.SetInStream(uo.dbInpStream);

  while (dbInpReader.Read(&dbMol)) {
#else
  RDKit::ForwardSDMolSupplier dbsuppl(uo.dbInpStream, takeOwnership, sanitize,
                                      removeHs);
  while (!dbsuppl.atEnd()) {
    std::unique_ptr<RDKit::ROMol> dbmptr(dbsuppl.next());
    if (!dbmptr) {
      continue;
    }
    RDKit::ROMol &dbMol = *dbmptr;
#endif
    ++molCount;

    // Keep track of the number of molecules processed so far
    if ((molCount % 10) == 0) {
      std::cerr << ".";
      if ((molCount % 500) == 0) {
        std::cerr << " " << molCount << " molecules" << std::endl;
      }
    }

    std::string dbName;
#ifndef USE_RDKIT
    dbName = dbMol.GetTitle();
#else
    dbMol.getPropIfPresent("_Name", dbName);
#endif
    if (dbName == "") {
      ss.str("");
      ss << "MOL_" << molCount;
      dbName = ss.str();
#ifndef USE_RDKIT
      dbMol.SetTitle(dbName);
#else
      dbMol.setProp("_Name", dbName);
#endif
    }

    // Create the set of Gaussians of database molecule
    GaussianVolume dbVolume;
    listAtomVolumes(dbMol, dbVolume);

    // Overlap with reference
    AlignmentInfo res;
    double bestScore(0.0);

    SolutionInfo bestSolution;
    if (uo.scoreOnly) {
      res.overlap = atomOverlap(refVolume, dbVolume);
      res.rotor[0] = 1.0;
      bestScore = getScore(uo.whichScore, res.overlap, refVolume.overlap,
                           dbVolume.overlap);
    } else {
      initOrientation(dbVolume);

      bestSolution = std::move(shapeit::alignVolumes(
          refVolume, dbVolume, uo.whichScore, uo.maxIter));
    }

#ifndef USE_RDKIT
    bestSolution.dbMol = dbMol;
#else
    bestSolution.dbMol = dbMol;
#endif
    bestSolution.refName = refName;
    bestSolution.dbName = dbName;

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

    // At this point the information of the solution is stored in bestSolution
    // Check if the result is better than the cutoff
    if (bestSolution.score < uo.cutOff) {
      continue;
    }

    // Post-process molecules
    if (uo.bestHits || !uo.molOutFile.empty()) {
      // Add the score properties
      setAllScores(bestSolution);

      if (!uo.scoreOnly) {
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
      }

      if (uo.bestHits) {
        bestHits->add(bestSolution);
      } else if (!uo.molOutFile.empty()) {
#ifndef USE_RDKIT
        dbOutWriter.Write(&(bestSolution.dbMol), uo.molOutStream);
#else
        dbOutWriter.write(bestSolution.dbMol);
#endif
      }
    }

    if ((uo.bestHits == 0) && !uo.scoreOutFile.empty()) {
      bestSolution.printScores(uo);
    }

#ifndef USE_RDKIT
    // Clear current molecule
    dbMol.Clear();
#endif
  }

  if (uo.bestHits) {
    if (!uo.molOutFile.empty()) {
      for (const auto molptr : bestHits->getBestList()) {
        if (molptr != nullptr) {
#ifndef USE_RDKIT
          dbOutWriter.Write(&(molptr->dbMol), uo.molOutStream);
#else
          dbOutWriter.write(molptr->dbMol);
#endif
        }
      }
      delete uo.molOutStream;
      uo.molOutStream = NULL;
    }
    if (!uo.scoreOutFile.empty()) {
      bestHits->writeScores(&uo);
      delete uo.scoreOutStream;
      uo.scoreOutStream = NULL;
    }
  }

  // Clear current streams
  if (uo.dbInpStream != NULL) {
    delete uo.dbInpStream;
    uo.dbInpStream = NULL;
  }
  if (uo.refInpStream != NULL) {
    delete uo.refInpStream;
    uo.refInpStream = NULL;
  }

  // Done processing database
  std::cerr << std::endl;
  std::cerr << "Processed " << molCount << " molecules" << std::endl;
  double tt = (double)(clock() - t0) / CLOCKS_PER_SEC;
  std::cerr << molCount << " molecules in " << tt << " seconds (";
  std::cerr << molCount / tt << " molecules per second)" << std::endl;

  // Cleanup local db volume
  refVolume.gaussians.clear();
  for (std::vector<std::vector<unsigned int> *>::iterator si =
           refVolume.childOverlaps.begin();
       si != refVolume.childOverlaps.end(); ++si) {
    if (*si != NULL) {
      delete *si;
      *si = NULL;
    }
  }

  exit(0);
}
