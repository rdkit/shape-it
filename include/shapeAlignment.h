/*******************************************************************************
shapeAlignment.h - Shape-it

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

#ifndef __SILICOSIT_SHAPEIT_SHAPEALIGNMENT_H__
#define __SILICOSIT_SHAPEIT_SHAPEALIGNMENT_H__

// General
#include <map>
#include <queue>
#include <vector>

// OpenBabel

// Shape-it
#include <alignmentInfo.h>
#include <atomGaussian.h>
#include <gaussianVolume.h>
#include <siMath.h>

typedef std::map<unsigned int, double *> MatrixMap;
typedef std::map<unsigned int, double *>::iterator MatIter;

class ShapeAlignment {
private:
  const GaussianVolume *_gRef;
  const GaussianVolume *_gDb;

  unsigned int _rAtoms;
  unsigned int _rGauss;
  unsigned int _dAtoms;
  unsigned int _dGauss;
  unsigned int _maxSize;
  unsigned int _maxIter;

  MatrixMap _matrixMap;

  double *_updateMatrixMap(const AtomGaussian &, const AtomGaussian &);

public:
  ShapeAlignment(const GaussianVolume &, const GaussianVolume &);
  ~ShapeAlignment(void);

  AlignmentInfo gradientAscent(SiMath::Vector);
  AlignmentInfo simulatedAnnealing(SiMath::Vector);

  void setMaxIterations(unsigned int);
};

#endif