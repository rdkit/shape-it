/*******************************************************************************
moleculeRotation.cpp - Shape-it

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

#include <moleculeRotation.h>

#ifndef USE_RDKIT
// OpenBabel
#include <openbabel/atom.h>
#include <openbabel/mol.h>

void positionMolecule(OpenBabel::OBMol &m, const Coordinate &centroid,
                      const SiMath::Matrix &rotation) {
  std::vector<OpenBabel::OBAtom *>::iterator i;
  for (OpenBabel::OBAtom *a = m.BeginAtom(i); a; a = m.NextAtom(i)) {
    // Translate point
    double x = a->GetX() - centroid.x;
    double y = a->GetY() - centroid.y;
    double z = a->GetZ() - centroid.z;

    // Rotate according to eigenvectors SVD
    a->SetVector(rotation[0][0] * x + rotation[1][0] * y + rotation[2][0] * z,
                 rotation[0][1] * x + rotation[1][1] * y + rotation[2][1] * z,
                 rotation[0][2] * x + rotation[1][2] * y + rotation[2][2] * z);
  }
  return;
}

void repositionMolecule(OpenBabel::OBMol &m, const SiMath::Matrix &rotation,
                        const Coordinate &centroid) {
  std::vector<OpenBabel::OBAtom *>::iterator i;
  for (OpenBabel::OBAtom *a = m.BeginAtom(i); a; a = m.NextAtom(i)) {
    // Get coordinates
    double x = a->GetX();
    double y = a->GetY();
    double z = a->GetZ();

    // Rotate according to eigenvectors SVD
    double xx = rotation[0][0] * x + rotation[0][1] * y + rotation[0][2] * z;
    double yy = rotation[1][0] * x + rotation[1][1] * y + rotation[1][2] * z;
    double zz = rotation[2][0] * x + rotation[2][1] * y + rotation[2][2] * z;

    a->SetVector(xx + centroid.x, yy + centroid.y, zz + centroid.z);
  }
  return;
}

void rotateMolecule(OpenBabel::OBMol &m, const SiMath::Vector &rotor) {
  // Build rotation matrix
  SiMath::Matrix rot(3, 3, 0.0);
  double r1 = rotor[1] * rotor[1];
  double r2 = rotor[2] * rotor[2];
  double r3 = rotor[3] * rotor[3];

  rot[0][0] = 1.0 - 2.0 * r2 - 2.0 * r3;
  rot[0][1] = 2.0 * (rotor[1] * rotor[2] - rotor[0] * rotor[3]);
  rot[0][2] = 2.0 * (rotor[1] * rotor[3] + rotor[0] * rotor[2]);
  rot[1][0] = 2.0 * (rotor[1] * rotor[2] + rotor[0] * rotor[3]);
  rot[1][1] = 1.0 - 2 * r3 - 2 * r1;
  rot[1][2] = 2.0 * (rotor[2] * rotor[3] - rotor[0] * rotor[1]);
  rot[2][0] = 2.0 * (rotor[1] * rotor[3] - rotor[0] * rotor[2]);
  rot[2][1] = 2.0 * (rotor[2] * rotor[3] + rotor[0] * rotor[1]);
  rot[2][2] = 1.0 - 2 * r2 - 2 * r1;

  std::vector<OpenBabel::OBAtom *>::iterator i;
  for (OpenBabel::OBAtom *a = m.BeginAtom(i); a; a = m.NextAtom(i)) {
    // Translate point
    double x = a->GetX();
    double y = a->GetY();
    double z = a->GetZ();

    // rotate according to eigenvectors SVD
    a->SetVector(rot[0][0] * x + rot[0][1] * y + rot[0][2] * z,
                 rot[1][0] * x + rot[1][1] * y + rot[1][2] * z,
                 rot[2][0] * x + rot[2][1] * y + rot[2][2] * z);
  }
  return;
}
#else
#include <Geometry/point.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/RDKitBase.h>

void positionMolecule(RDKit::ROMol &m, const Coordinate &centroid,
                      const SiMath::Matrix &rotation) {
  RDKit::Conformer &conf = m.getConformer();
  RDGeom::Point3D rdcentroid(centroid.x, centroid.y, centroid.z);
  for (unsigned int i = 0; i < m.getNumAtoms(); ++i) {
    RDGeom::Point3D tp = conf.getAtomPos(i) - rdcentroid;
    conf.setAtomPos(
        i, RDGeom::Point3D(rotation[0][0] * tp.x + rotation[1][0] * tp.y +
                               rotation[2][0] * tp.z,
                           rotation[0][1] * tp.x + rotation[1][1] * tp.y +
                               rotation[2][1] * tp.z,
                           rotation[0][2] * tp.x + rotation[1][2] * tp.y +
                               rotation[2][2] * tp.z));
  }
}

void repositionMolecule(RDKit::ROMol &m, const SiMath::Matrix &rotation,
                        const Coordinate &centroid) {
  RDKit::Conformer &conf = m.getConformer();
  RDGeom::Point3D rdcentroid(centroid.x, centroid.y, centroid.z);
  for (unsigned int i = 0; i < m.getNumAtoms(); ++i) {
    RDGeom::Point3D tp = conf.getAtomPos(i);
    conf.setAtomPos(
        i, RDGeom::Point3D(rotation[0][0] * tp.x + rotation[0][1] * tp.y +
                               rotation[0][2] * tp.z,
                           rotation[1][0] * tp.x + rotation[1][1] * tp.y +
                               rotation[1][2] * tp.z,
                           rotation[2][0] * tp.x + rotation[2][1] * tp.y +
                               rotation[2][2] * tp.z));
    conf.getAtomPos(i) += rdcentroid;
  }
  return;
}

void rotateMolecule(RDKit::ROMol &m, const SiMath::Vector &rotor) {
  RDKit::Conformer &conf = m.getConformer();
  // Build rotation matrix
  SiMath::Matrix rot(3, 3, 0.0);
  double r1 = rotor[1] * rotor[1];
  double r2 = rotor[2] * rotor[2];
  double r3 = rotor[3] * rotor[3];

  rot[0][0] = 1.0 - 2.0 * r2 - 2.0 * r3;
  rot[0][1] = 2.0 * (rotor[1] * rotor[2] - rotor[0] * rotor[3]);
  rot[0][2] = 2.0 * (rotor[1] * rotor[3] + rotor[0] * rotor[2]);
  rot[1][0] = 2.0 * (rotor[1] * rotor[2] + rotor[0] * rotor[3]);
  rot[1][1] = 1.0 - 2 * r3 - 2 * r1;
  rot[1][2] = 2.0 * (rotor[2] * rotor[3] - rotor[0] * rotor[1]);
  rot[2][0] = 2.0 * (rotor[1] * rotor[3] - rotor[0] * rotor[2]);
  rot[2][1] = 2.0 * (rotor[2] * rotor[3] + rotor[0] * rotor[1]);
  rot[2][2] = 1.0 - 2 * r2 - 2 * r1;

  for (unsigned int i = 0; i < m.getNumAtoms(); ++i) {
    RDGeom::Point3D tp = conf.getAtomPos(i);
    conf.setAtomPos(
        i, RDGeom::Point3D(
               rot[0][0] * tp.x + rot[0][1] * tp.y + rot[0][2] * tp.z,
               rot[1][0] * tp.x + rot[1][1] * tp.y + rot[1][2] * tp.z,
               rot[2][0] * tp.x + rot[2][1] * tp.y + rot[2][2] * tp.z));
  }
}

#endif