/*******************************************************************************
printHeader.cpp - Shape-it

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

#include <printHeader.h>
#ifndef USE_RDKIT
// OpenBabel
#include <openbabel/mol.h>
#else
#include <RDGeneral/versions.h>
#endif

void printHeader() {
  std::cerr << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
               "++++++++++++"
            << std::endl;
  std::cerr << "  Shape-it v" << SHAPEIT_VERSION << "." << SHAPEIT_RELEASE
            << "." << SHAPEIT_SUBRELEASE << " | ";
  std::cerr << __DATE__ " " << __TIME__ << std::endl;
  std::cerr << std::endl;
  std::cerr << "  -> GCC:       " << __VERSION__ << std::endl;
#ifndef USE_RDKIT
  std::cerr << "  -> OpenBabel: " << BABEL_VERSION << std::endl;
#else
  std::cerr << "  -> RDKit: " << RDKit::rdkitVersion << std::endl;
#endif
  std::cerr << std::endl;
  std::cerr
      << R"DOC(  Copyright 2012-2021 by Silicos-it, a division of Imacosi BVBA, Hans De Winter,
     and the Shape-it contributors

  Shape-it is free software: you can redistribute it and/or modify
  it under the terms of the MIT license.

  Shape-it is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  MIT License for more details.
)DOC" << std::endl;
#ifndef USE_RDKIT
  std::cerr << "  Shape-it is linked against OpenBabel version 2." << std::endl;
  std::cerr
      << "  OpenBabel is free software; you can redistribute it and/or modify"
      << std::endl;
  std::cerr << "  it under the terms of the GNU General Public License as "
               "published by"
            << std::endl;
  std::cerr << "  the Free Software Foundation version 2 of the License."
            << std::endl;
#else
  std::cerr << "  Shape-it is linked against the RDKit (https://www.rdkit.org)."
            << std::endl;
#endif
  std::cerr << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
               "++++++++++++"
            << std::endl;
  std::cerr << std::endl;
  return;
}
