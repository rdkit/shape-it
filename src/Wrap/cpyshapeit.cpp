/*******************************************************************************

Copyright 2021 by Greg Landrum and the Shape-it contributors

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

***********************************************************************/

#include <boost/python.hpp>

#include "alignLib.h"
#include "options.h"

namespace python = boost::python;
using namespace RDKit;

namespace {
double alignMol(const ROMol &ref, ROMol &probe, const std::string &whichScore,
                double maxIter, double cutoff) {
  auto sinfo = shapeit::alignMols(ref, probe, whichScore, maxIter, cutoff);
  const Conformer &conf = sinfo.dbMol.getConformer();
  probe.clearConformers();
  probe.addConformer(new Conformer(conf));
  return sinfo.score;
}
} // namespace

void wrap_pyshapeit() {
  python::def("AlignMol", &alignMol,
              (python::arg("ref"), python::arg("probe"),
               python::arg("whichScore") = tanimoto,
               python::arg("maxIter") = 0., python::arg("cutoff") = 0.),
              "aligns probe to ref, probe is modified");
}

BOOST_PYTHON_MODULE(cpyshapeit) { wrap_pyshapeit(); }
