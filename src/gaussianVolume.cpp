/*******************************************************************************
gaussianVolume.cpp - Shape-it

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

#include <gaussianVolume.h>

#ifndef USE_RDKIT
// OpenBabel
#include <openbabel/atom.h>
#include <openbabel/data.h>
#include <openbabel/elements.h>
#include <openbabel/mol.h>
#else
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/ROMol.h>
#endif

// OpenBabel

GaussianVolume::GaussianVolume(void)
    : volume(0.0), overlap(0.0), centroid(0.0, 0.0, 0.0), rotation(3, 3, 0.0),
      gaussians(), childOverlaps(), levels() {}

GaussianVolume::~GaussianVolume(void) {}

double GAlpha(unsigned int an) {
  switch (an) {
  case 1: ///< H
    return 1.679158285;
    break;
  case 3: ///< Li
    return 0.729980658;
    break;
  case 5: ///< B
    return 0.604496983;
    break;
  case 6: ///< C
    return 0.836674025;
    break;
  case 7: ///< N
    return 1.006446589;
    break;
  case 8: ///< O
    return 1.046566798;
    break;
  case 9: ///< F
    return 1.118972618;
    break;
  case 11: ///< Na
    return 0.469247983;
    break;
  case 12: ///< Mg
    return 0.807908026;
    break;
  case 14: ///< Si
    return 0.548296583;
    break;
  case 15: ///< P
    return 0.746292571;
    break;
  case 16: ///< S
    return 0.746292571;
    break;
  case 17: ///< Cl
    return 0.789547080;
    break;
  case 19: ///< K
    return 0.319733941;
    break;
  case 20: ///< Ca
    return 0.604496983;
    break;
  case 26: ///< Fe
    return 1.998337133;
    break;
  case 29: ///< Cu
    return 1.233667312;
    break;
  case 30: ///< Zn
    return 1.251481772;
    break;
  case 35: ///< Br
    return 0.706497569;
    break;
  case 53: ///< I
    return 0.616770720;
    break;
  default: ///< *
    return 1.074661303;
  }
  return 1.074661303;
};

namespace {
#ifndef USE_RDKIT

unsigned int initFromMol(const Molecule &mol, GaussianVolume &gv) {
  // Prepare the vector to store the atom and overlap volumes;
  unsigned int N = 0;
  for (unsigned int i = 1; i <= mol.NumAtoms(); ++i) {
    const OpenBabel::OBAtom *a = mol.GetAtom(i);
    if (a->GetAtomicNum() == 1) {
      continue;
    } else {
      ++N;
    }
  }
  gv.gaussians.resize(N);
  gv.childOverlaps.resize(N);
  gv.levels.push_back(N); // first level
  gv.volume = 0.0;
  gv.centroid.x = 0;
  gv.centroid.y = 0;
  gv.centroid.z = 0;
  int atomIndex = 0; // keeps track of the atoms processed so far
  int vecIndex = N;  // keeps track of the last element added to the vectors
  for (unsigned int i = 1; i <= mol.NumAtoms(); ++i) {
    const OpenBabel::OBAtom *a = mol.GetAtom(i);
    // Skip hydrogens
    if (a->GetAtomicNum() == 1) {
      continue;
    }

    // First atom append self to the list
    // Store it at [index]
    gv.gaussians[atomIndex].center.x = a->GetX();
    gv.gaussians[atomIndex].center.y = a->GetY();
    gv.gaussians[atomIndex].center.z = a->GetZ();
    gv.gaussians[atomIndex].alpha = GAlpha(a->GetAtomicNum());
    gv.gaussians[atomIndex].C = GCI;
    double radius = OpenBabel::OBElements::GetVdwRad(a->GetAtomicNum());
    gv.gaussians[atomIndex].volume =
        (4.0 * PI / 3.0) * radius * radius * radius;
    ++atomIndex;
  }
  return N;
}
#else
unsigned int initFromMol(const Molecule &mol, GaussianVolume &gv) {
  // Prepare the vector to store the atom and overlap volumes;
  unsigned int N = 0;
  for (const auto a : mol.atoms()) {
    if (a->getAtomicNum() != 1) {
      ++N;
    }
  }
  gv.gaussians.resize(N);
  gv.childOverlaps.resize(N);
  gv.levels.push_back(N); // first level
  gv.volume = 0.0;
  gv.centroid.x = 0;
  gv.centroid.y = 0;
  gv.centroid.z = 0;
  int atomIndex = 0; // keeps track of the atoms processed so far
  int vecIndex = N;  // keeps track of the last element added to the vectors
  const auto &conf = mol.getConformer();
  for (const auto a : mol.atoms()) {
    // Skip hydrogens
    if (a->getAtomicNum() == 1) {
      continue;
    }

    // First atom append self to the list
    // Store it at [index]
    const auto &p = conf.getAtomPos(a->getIdx());
    gv.gaussians[atomIndex].center.x = p.x;
    gv.gaussians[atomIndex].center.y = p.y;
    gv.gaussians[atomIndex].center.z = p.z;
    gv.gaussians[atomIndex].alpha = GAlpha(a->getAtomicNum());
    gv.gaussians[atomIndex].C = GCI;
    // double radius = et.GetVdwRad(a->GetAtomicNum());
    double radius =
        RDKit::PeriodicTable::getTable()->getRvdw(a->getAtomicNum());
    gv.gaussians[atomIndex].volume =
        (4.0 * PI / 3.0) * radius * radius * radius;
    ++atomIndex;
  }
  return N;
}

#endif

} // namespace

void listAtomVolumes(const Molecule &mol, GaussianVolume &gv) {
  // Prepare the vector to store the atom and overlap volumes;
  unsigned int N = initFromMol(mol, gv);

  // Vector to keep track of parents of an overlap
  std::vector<std::pair<unsigned int, unsigned int>> parents(N);

  // Create a vector to keep track of the overlaps
  // Overlaps are stored as sets
  std::vector<std::set<unsigned int> *> overlaps(N);
  std::set<unsigned int>::iterator setIter;

  // Start by iterating over the single atoms and build map of overlaps
  int vecIndex = N; // keeps track of the last element added to the vectors
  for (unsigned int atomIndex = 0; atomIndex < N; ++atomIndex) {
    // Add empty child overlaps
    std::vector<unsigned int> *vec = new std::vector<unsigned int>();
    gv.childOverlaps[atomIndex] = vec;

    // Update volume and centroid
    gv.volume += gv.gaussians[atomIndex].volume;
    gv.centroid.x +=
        gv.gaussians[atomIndex].volume * gv.gaussians[atomIndex].center.x;
    gv.centroid.y +=
        gv.gaussians[atomIndex].volume * gv.gaussians[atomIndex].center.y;
    gv.centroid.z +=
        gv.gaussians[atomIndex].volume * gv.gaussians[atomIndex].center.z;

    // Add new empty set of possible overlaps
    std::set<unsigned int> *tmp = new std::set<unsigned int>();
    overlaps[atomIndex] = tmp;

    // Loop over the current list of processed atoms and add overlaps
    for (int i = 0; i < atomIndex; ++i) {
      // Create overlap gaussian
      AtomGaussian ga =
          atomIntersection(gv.gaussians[i], gv.gaussians[atomIndex]);

      // Check if the atom-atom overlap volume is large enough
      if (ga.volume / (gv.gaussians[i].volume + gv.gaussians[atomIndex].volume -
                       ga.volume) <
          EPS) {
        continue;
      }

      // Add gaussian volume, and empty overlap set
      gv.gaussians.push_back(ga);
      std::vector<unsigned int> *vec = new std::vector<unsigned int>();
      gv.childOverlaps.push_back(vec);

      // Update local variables of parents and possible overlaps
      parents.push_back(std::make_pair(i, atomIndex));
      std::set<unsigned int> *dummy = new std::set<unsigned int>();
      overlaps.push_back(dummy);

      // Update volume and centroid (negative contribution of atom-atom overlap)
      gv.volume -= ga.volume;
      gv.centroid.x -= ga.volume * ga.center.x;
      gv.centroid.y -= ga.volume * ga.center.y;
      gv.centroid.z -= ga.volume * ga.center.z;

      // Update overlap information of the parent atom
      overlaps[i]->insert(atomIndex);
      gv.childOverlaps[i]->push_back(vecIndex);

      // Move to next index in vector
      ++vecIndex;
    }
  }

  // Position in list of gaussians where atom gaussians end
  unsigned int startLevel = N;
  unsigned int nextLevel = gv.gaussians.size();

  // Update level information
  gv.levels.push_back(nextLevel);

  // Loop overall possible levels of overlaps from 2 to 6
  for (unsigned int l = 2; l < LEVEL; ++l) {
    // List of atom-atom overlaps is made => gv.gaussians[startLevel ..
    // nextLevel-1]; Now update the overlap lists for each overlap in this level
    // Create the next overlap Gaussian
    // And add it to the vector of overlaps
    for (unsigned int i = startLevel; i < nextLevel; ++i) {
      // Parent indices
      unsigned int a1 = parents[i].first;
      unsigned int a2 = parents[i].second;

      // Append volume to end of overlap vector
      // Add new empty set
      std::set<unsigned int> *tmp = overlaps[i];
      std::set_intersection(
          overlaps[a1]->begin(), overlaps[a1]->end(), overlaps[a2]->begin(),
          overlaps[a2]->end(),
          std::insert_iterator<std::set<unsigned int>>(*tmp, tmp->begin()));

      // Check if the overlap list is empty
      if (overlaps[i]->empty()) {
        continue;
      }

      // Get the possible overlaps from the parent gaussians
      // and create the new overlap volume
      for (setIter = overlaps[i]->begin(); setIter != overlaps[i]->end();
           ++setIter) {
        if (*setIter <= a2) {
          continue;
        }

        // Create a new overlap gaussian
        AtomGaussian ga =
            atomIntersection(gv.gaussians[i], gv.gaussians[*setIter]);

        // Check if the volume is large enough
        if (ga.volume / (gv.gaussians[i].volume +
                         gv.gaussians[*setIter].volume - ga.volume) <
            EPS) {
          continue;
        }

        gv.gaussians.push_back(ga);
        std::vector<unsigned int> *vec = new std::vector<unsigned int>();
        gv.childOverlaps.push_back(vec);

        // Update local variables
        parents.push_back(std::make_pair(i, *setIter));
        std::set<unsigned int> *tmp = new std::set<unsigned int>();
        overlaps.push_back(tmp);

        // Update volume, centroid and moments
        // Overlaps consisting of an even number of atoms have a negative
        // contribution
        if ((ga.nbr % 2) == 0) {
          // Update volume and centroid
          gv.volume -= ga.volume;
          gv.centroid.x -= ga.volume * ga.center.x;
          gv.centroid.y -= ga.volume * ga.center.y;
          gv.centroid.z -= ga.volume * ga.center.z;
        } else {
          // Update volume and centroid
          gv.volume += ga.volume;
          gv.centroid.x += ga.volume * ga.center.x;
          gv.centroid.y += ga.volume * ga.center.y;
          gv.centroid.z += ga.volume * ga.center.z;
        }

        // Update child list of the first
        gv.childOverlaps[i]->push_back(vecIndex);

        // Move to next index in vector
        ++vecIndex;
      }
    }

    // Update levels
    startLevel = nextLevel;
    nextLevel = gv.gaussians.size();

    // Update level information
    gv.levels.push_back(nextLevel);
  }

  // cleanup current set of computed overlaps
  for (std::vector<std::set<unsigned int> *>::iterator si = overlaps.begin();
       si != overlaps.end(); ++si) {
    if (*si != NULL) {
      delete *si;
      *si = NULL;
    }
  }

  parents.clear();

  // Update self-overlap
  gv.overlap = atomOverlap(gv, gv);

  return;
}

void initOrientation(GaussianVolume &gv) {
  double x(0.0), y(0.0), z(0.0);

  // Scale centroid and moments with self volume
  gv.centroid.x /= gv.volume;
  gv.centroid.y /= gv.volume;
  gv.centroid.z /= gv.volume;

  // Compute moments of inertia from mass matrix
  SiMath::Matrix mass(3, 3, 0.0);

  // Loop over all gaussians
  for (std::vector<AtomGaussian>::iterator i = gv.gaussians.begin();
       i != gv.gaussians.end(); ++i) {
    // Translate to center
    i->center.x -= gv.centroid.x;
    i->center.y -= gv.centroid.y;
    i->center.z -= gv.centroid.z;

    x = i->center.x;
    y = i->center.y;
    z = i->center.z;

    if ((i->nbr % 2) == 0) {
      // Update upper triangle
      mass[0][0] -= i->volume * x * x;
      mass[0][1] -= i->volume * x * y;
      mass[0][2] -= i->volume * x * z;
      mass[1][1] -= i->volume * y * y;
      mass[1][2] -= i->volume * y * z;
      mass[2][2] -= i->volume * z * z;
    } else {
      // Update upper triangle
      mass[0][0] += i->volume * x * x;
      mass[0][1] += i->volume * x * y;
      mass[0][2] += i->volume * x * z;
      mass[1][1] += i->volume * y * y;
      mass[1][2] += i->volume * y * z;
      mass[2][2] += i->volume * z * z;
    }
  }

  // Set lower triangle
  mass[1][0] = mass[0][1];
  mass[2][0] = mass[0][2];
  mass[2][1] = mass[1][2];

  // Normalize mass matrix
  mass /= gv.volume;

  // Compute SVD of the mass matrix
  SiMath::SVD svd(mass, true, true);
  gv.rotation = svd.getU();

  double det = gv.rotation[0][0] * gv.rotation[1][1] * gv.rotation[2][2] +
               gv.rotation[2][1] * gv.rotation[1][0] * gv.rotation[0][2] +
               gv.rotation[0][1] * gv.rotation[1][2] * gv.rotation[2][0] -
               gv.rotation[0][0] * gv.rotation[1][2] * gv.rotation[2][1] -
               gv.rotation[1][1] * gv.rotation[2][0] * gv.rotation[0][2] -
               gv.rotation[2][2] * gv.rotation[0][1] * gv.rotation[1][0];

  // Check if it is a rotation matrix and not a mirroring
  if (det < 0) {
    // Switch sign of third column
    gv.rotation[0][2] = -gv.rotation[0][2];
    gv.rotation[1][2] = -gv.rotation[1][2];
    gv.rotation[2][2] = -gv.rotation[2][2];
  }

  // Rotate all gaussians
  for (std::vector<AtomGaussian>::iterator i = gv.gaussians.begin();
       i != gv.gaussians.end(); ++i) {
    x = i->center.x;
    y = i->center.y;
    z = i->center.z;
    i->center.x =
        gv.rotation[0][0] * x + gv.rotation[1][0] * y + gv.rotation[2][0] * z;
    i->center.y =
        gv.rotation[0][1] * x + gv.rotation[1][1] * y + gv.rotation[2][1] * z;
    i->center.z =
        gv.rotation[0][2] * x + gv.rotation[1][2] * y + gv.rotation[2][2] * z;
  }

  return;
}

double atomOverlap(const GaussianVolume &gRef, const GaussianVolume &gDb) {
  // Create a queue to hold the pairs to process
  std::queue<std::pair<unsigned int, unsigned int>> processQueue;

  // loop over the single atom volumes of both molecules and make the
  // combinations
  unsigned int N1(gRef.levels[0]);
  unsigned int N2(gDb.levels[0]);

  double Cij(0.0), Vij(0.0);

  double dx(0.0), dy(0.0), dz(0.0);

  std::vector<unsigned int> *d1(NULL), *d2(NULL);
  std::vector<unsigned int>::iterator it1;

  // Overlap volume
  double overlapVol(0.0);

  // First compute atom-atom overlaps
  for (unsigned int i(0); i < N1; ++i) {
    for (unsigned int j(0); j < N2; ++j) {
      // Scaling constant
      Cij = gRef.gaussians[i].alpha * gDb.gaussians[j].alpha /
            (gRef.gaussians[i].alpha + gDb.gaussians[j].alpha);

      // Variables to store sum and difference of components
      dx = (gRef.gaussians[i].center.x - gDb.gaussians[j].center.x);
      dx *= dx;
      dy = (gRef.gaussians[i].center.y - gDb.gaussians[j].center.y);
      dy *= dy;
      dz = (gRef.gaussians[i].center.z - gDb.gaussians[j].center.z);
      dz *= dz;

      // Compute overlap volume
      Vij = gRef.gaussians[i].C * gDb.gaussians[j].C *
            pow(PI / (gRef.gaussians[i].alpha + gDb.gaussians[j].alpha), 1.5) *
            exp(-Cij * (dx + dy + dz));

      // Check if overlap is sufficient enough
      if (Vij / (gRef.gaussians[i].volume + gDb.gaussians[j].volume - Vij) <
          EPS) {
        continue;
      }

      // Add to overlap volume
      overlapVol += Vij;

      // Loop over child nodes and add to queue
      d1 = gRef.childOverlaps[i];
      d2 = gDb.childOverlaps[j];

      // First add (i,child(j))
      if (d2 != NULL) {
        for (it1 = d2->begin(); it1 != d2->end(); ++it1) {
          processQueue.push(std::make_pair(i, *it1));
        }
      }
      // Second add (child(i,j))
      if (d1 != NULL) {
        for (it1 = d1->begin(); it1 != d1->end(); ++it1) {
          // add (child(i),j)
          processQueue.push(std::make_pair(*it1, j));
        }
      }
    }
  }

  while (!processQueue.empty()) {
    // Get next element from queue
    std::pair<unsigned int, unsigned int> nextPair = processQueue.front();
    processQueue.pop();

    unsigned int i = nextPair.first;
    unsigned int j = nextPair.second;

    // Scaling constant
    Cij = gRef.gaussians[i].alpha * gDb.gaussians[j].alpha /
          (gRef.gaussians[i].alpha + gDb.gaussians[j].alpha);

    // Variables to store sum and difference of components
    dx = (gRef.gaussians[i].center.x - gDb.gaussians[j].center.x);
    dx *= dx;
    dy = (gRef.gaussians[i].center.y - gDb.gaussians[j].center.y);
    dy *= dy;
    dz = (gRef.gaussians[i].center.z - gDb.gaussians[j].center.z);
    dz *= dz;

    // Compute overlap volume
    Vij = gRef.gaussians[i].C * gDb.gaussians[j].C *
          pow(PI / (gRef.gaussians[i].alpha + gDb.gaussians[j].alpha), 1.5) *
          exp(-Cij * (dx + dy + dz));

    // Check if overlap is sufficient enough
    if (Vij / (gRef.gaussians[i].volume + gDb.gaussians[j].volume - Vij) <
        EPS) {
      continue;
    }

    // Even number of overlap atoms => addition to volume
    // Odd number => substraction
    if ((gRef.gaussians[i].nbr + gDb.gaussians[j].nbr) % 2 == 0) {
      overlapVol += Vij;
    } else {
      overlapVol -= Vij;
    }

    // Loop over child nodes and add to queue
    d1 = gRef.childOverlaps[i];
    d2 = gDb.childOverlaps[j];
    if (d1 != NULL && gRef.gaussians[i].nbr > gDb.gaussians[j].nbr) {
      for (it1 = d1->begin(); it1 != d1->end(); ++it1) {
        // Add (child(i),j)
        processQueue.push(std::make_pair(*it1, j));
      }
    } else {
      // First add (i,child(j))
      if (d2 != NULL) {
        for (it1 = d2->begin(); it1 != d2->end(); ++it1) {
          processQueue.push(std::make_pair(i, *it1));
        }
      }
      if (d1 != NULL && gDb.gaussians[j].nbr - gRef.gaussians[i].nbr < 2) {
        for (it1 = d1->begin(); it1 != d1->end(); ++it1) {
          // add (child(i),j)
          processQueue.push(std::make_pair(*it1, j));
        }
      }
    }
  }

  return overlapVol;
}

double getScore(const std::string &id, double Voa, double Vra, double Vda) {
  // set the score by which molecules are being compared
  if (id == tanimoto) {
    return Voa / (Vra + Vda - Voa);
  } else if (id == tversky_ref) {
    return Voa / (0.95 * Vra + 0.05 * Vda);
  } else if (id == tversky_db) {
    return Voa / (0.05 * Vra + 0.95 * Vda);
  }

  return 0.0;
}

void checkVolumes(const GaussianVolume &gRef, const GaussianVolume &gDb,
                  AlignmentInfo &res) {
  if (res.overlap > gRef.overlap) {
    res.overlap = gRef.overlap;
  }
  if (res.overlap > gDb.overlap) {
    res.overlap = gDb.overlap;
  }
  return;
}
