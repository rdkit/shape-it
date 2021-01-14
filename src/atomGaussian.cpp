/*******************************************************************************
atomGaussian.cpp - Shape-it

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

#include <atomGaussian.h>

AtomGaussian::AtomGaussian(void)
    : center(0.0, 0.0, 0.0), alpha(0.0), volume(0.0), C(0.0), nbr(0) {}

AtomGaussian::~AtomGaussian(void) {}

AtomGaussian atomIntersection(AtomGaussian &a, AtomGaussian &b) {
  AtomGaussian c;

  // new alpha
  c.alpha = a.alpha + b.alpha;

  // new center
  c.center.x = (a.alpha * a.center.x + b.alpha * b.center.x) / c.alpha;
  c.center.y = (a.alpha * a.center.y + b.alpha * b.center.y) / c.alpha;
  c.center.z = (a.alpha * a.center.z + b.alpha * b.center.z) / c.alpha;

  // self-volume
  double d = (a.center.x - b.center.x) * (a.center.x - b.center.x) +
             (a.center.y - b.center.y) * (a.center.y - b.center.y) +
             (a.center.z - b.center.z) * (a.center.z - b.center.z);

  c.C = a.C * b.C * exp(-a.alpha * b.alpha / c.alpha * d);

  double scale = PI / (c.alpha);

  c.volume = c.C * scale * sqrt(scale);

  // set the number of gaussians
  c.nbr = a.nbr + b.nbr;

  return c;
}
