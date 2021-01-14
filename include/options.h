/*******************************************************************************
options.h - Shape-it

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

#ifndef __SILICOSIT_SHAPEIT_OPTIONS_H__
#define __SILICOSIT_SHAPEIT_OPTIONS_H__

// General
#include <fstream>
#include <iostream>
#include <string>

// Shape-it
const std::string tanimoto = "Shape-it::Tanimoto";
const std::string tversky_ref = "Shape-it::Tversky_Ref";
const std::string tversky_db = "Shape-it::Tversky_Db";

class Options {
public:
  std::string format;     //  -f
  std::string refInpFile; //  -r  --reference
  std::ifstream *refInpStream;

  std::string dbInpFile; //  -d  --dbase
  std::ifstream *dbInpStream;

  std::string molOutFile; //  -o  --out
  std::ofstream *molOutStream;

  std::string scoreOutFile; //  -s  --score
  std::ofstream *scoreOutStream;

  unsigned int bestHits; //      --best
  double cutOff;         //      --cutOff
  double maxIter;        //      --maxIterations

  std::string whichScore; //      --rankBy

  bool scoreOnly; //      --scoreOnly
  bool showRef;   //      --noRef

  bool version; //  -v  --version
  bool help;    //  -h  --help

  Options(void);
  ~Options(void);

  std::string print(void) const;
};

#endif
