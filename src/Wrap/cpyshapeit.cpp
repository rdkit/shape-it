//
//  Copyright (C) 2021 Greg Landrum
//
#include "alignLib.h"
#include <boost/python.hpp>

namespace python = boost::python;
namespace {
std::string hello() { return "hello world"; }
} // namespace

void wrap_pyshapeit() { python::def("hello", &hello, "hello world"); }

BOOST_PYTHON_MODULE(cpyshapeit) { wrap_pyshapeit(); }
