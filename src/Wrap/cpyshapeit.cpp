//
//  Copyright (C) 2021 Greg Landrum
//
#include <boost/python.hpp>

#include "alignLib.h"
#include "options.h"

namespace python = boost::python;
using namespace RDKit;

namespace {
std::string hello() { return "hello world"; }
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
  python::def("hello", &hello, "hello world");
  python::def("AlignMol", &alignMol,
              (python::arg("ref"), python::arg("probe"),
               python::arg("whichScore") = tanimoto,
               python::arg("maxIter") = 0., python::arg("cutoff") = 0.),
              "aligns probe to ref, probe is modified");
}

BOOST_PYTHON_MODULE(cpyshapeit) { wrap_pyshapeit(); }
