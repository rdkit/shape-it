/*******************************************************************************
options.cpp - Shape-it

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

#include <options.h>
#include <sstream>

Options::Options(void) {
  refInpFile = "";
  refInpStream = NULL;

  dbInpFile = "";
  dbInpStream = NULL;

  molOutFile = "";
  molOutStream = NULL;

  scoreOutFile = "";
  scoreOutStream = NULL;

  bestHits = 0;
  cutOff = 0.0;
  maxIter = 0;

  whichScore = tanimoto;

  scoreOnly = false;
  showRef = true;

  version = false;
  help = false;
}

Options::~Options(void) {
  // reference input
  if (!refInpFile.empty()) {
    refInpFile = "";
  };
  if (refInpStream) {
    delete refInpStream;
    refInpStream = NULL;
  };

  // Database input
  if (!dbInpFile.empty()) {
    dbInpFile = "";
  };
  if (dbInpStream) {
    delete dbInpStream;
    dbInpStream = NULL;
  };

  // Molecule output
  if (!molOutFile.empty()) {
    molOutFile = "";
  };
  if (molOutStream) {
    delete molOutStream;
    molOutStream = NULL;
  };

  // Score output
  if (!scoreOutFile.empty()) {
    scoreOutFile = "";
  };
  if (scoreOutStream) {
    delete scoreOutStream;
    scoreOutStream = NULL;
  };
}

std::string Options::print(void) const {
  std::ostringstream os;
  os << std::endl;
  os << "COMMAND_LINE OPTIONS:" << std::endl;
  os << std::endl;
  os << "  -> Reference file:    " << refInpFile << std::endl;
  os << "  -> Database file:     " << dbInpFile << std::endl;
  os << "  -> Output file:       " << (molOutFile.empty() ? "no" : molOutFile)
     << std::endl;
  os << "  -> Scores file:       "
     << (scoreOutFile.empty() ? "no" : scoreOutFile) << std::endl;
  os << "  -> Best hits:         ";
  if (bestHits) {
    os << bestHits << std::endl;
  } else {
    os << "no" << std::endl;
  }
  os << "  -> Scoring only:      " << (scoreOnly ? "yes" : "no") << std::endl;
  os << "  -> Extra iterations:  ";
  if (maxIter) {
    os << maxIter << std::endl;
  } else {
    os << "no" << std::endl;
  }
  os << "  -> Rank by:           " << whichScore << std::endl;
  os << "  -> Cutoff:            ";
  if (cutOff) {
    os << cutOff << std::endl;
  } else {
    os << "no" << std::endl;
  }
  os << "  -> Output reference   " << (showRef ? "yes" : "no") << std::endl;

  os << std::endl;
  std::string r = os.str();
  return r;
}
