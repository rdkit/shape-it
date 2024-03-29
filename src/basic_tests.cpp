//
//  Copyright (C) 2021 Greg Landrum and the Shape-it contributors
/*
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
*/

#include "catch2/catch_all.hpp"

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/RDKitBase.h>
#include <alignLib.h>

using namespace RDKit;

TEST_CASE("align mols", "[library][molecules]") {
  std::string ref = R"CTAB(3l5u_lig_ZEC
                    3D
 Structure written by MMmdl.
 20 21  0  0  1  0            999 V2000
   15.0500  -34.9220  -18.1430 O   0  0  0  0  0  0
   14.9110  -34.7040  -19.4790 C   0  0  0  0  0  0
   14.7350  -35.7750  -20.3500 C   0  0  0  0  0  0
   14.6060  -35.5430  -21.7160 C   0  0  0  0  0  0
   14.3620  -36.6080  -23.0370 S   0  0  0  0  0  0
   14.3210  -35.3400  -24.1850 C   0  0  0  0  0  0
   14.0030  -35.5290  -25.8940 S   0  0  0  0  0  0
   15.1750  -34.7990  -26.7570 N   0  0  0  0  0  0
   12.6760  -34.9070  -26.2170 O   0  0  0  0  0  0
   13.9630  -36.9930  -26.2180 O   0  0  0  0  0  0
   14.4970  -34.1790  -23.5590 N   0  0  0  0  0  0
   14.6510  -34.2470  -22.2350 C   0  0  0  0  0  0
   14.8270  -33.1870  -21.3480 C   0  0  0  0  0  0
   14.9560  -33.4070  -19.9800 C   0  0  0  0  0  0
   15.1610  -34.0820  -17.6920 H   0  0  0  0  0  0
   14.6990  -36.7830  -19.9650 H   0  0  0  0  0  0
   15.1440  -34.8130  -27.7660 H   0  0  0  0  0  0
   15.9340  -34.3310  -26.2820 H   0  0  0  0  0  0
   14.8640  -32.1730  -21.7190 H   0  0  0  0  0  0
   15.0910  -32.5720  -19.3090 H   0  0  0  0  0  0
  1  2  1  0  0  0
  1 15  1  0  0  0
  2  3  2  0  0  0
  2 14  1  0  0  0
  3  4  1  0  0  0
  3 16  1  0  0  0
  4  5  1  0  0  0
  4 12  2  0  0  0
  5  6  1  0  0  0
  6  7  1  0  0  0
  6 11  2  0  0  0
  7  8  1  0  0  0
  7  9  2  0  0  0
  7 10  2  0  0  0
  8 17  1  0  0  0
  8 18  1  0  0  0
 11 12  1  0  0  0
 12 13  1  0  0  0
 13 14  2  0  0  0
 13 19  1  0  0  0
 14 20  1  0  0  0
M  END)CTAB";
  std::string probe = R"CTAB(3hof_lig_DHC
                    3D
 Structure written by MMmdl.
 20 20  0  0  1  0            999 V2000
   14.6290  -34.5170  -18.4190 C   0  0  0  0  0  0
   15.6070  -34.6620  -17.5400 O   0  0  0  0  0  0
   14.9220  -34.5200  -19.8370 C   0  0  0  0  0  0
   14.7370  -35.7220  -20.3520 C   0  0  0  0  0  0
   14.9680  -35.9740  -21.7740 C   0  0  0  0  0  0
   14.8780  -34.9380  -22.6930 C   0  0  0  0  0  0
   15.1020  -35.2380  -24.0360 C   0  0  0  0  0  0
   15.4390  -36.6310  -24.4550 C   0  0  0  0  0  0
   15.5160  -37.6070  -23.4830 C   0  0  0  0  0  0
   15.2760  -37.2740  -22.1560 C   0  0  0  0  0  0
   15.6830  -36.9520  -25.7670 O   0  0  0  0  0  0
   15.0160  -34.2570  -24.9550 O   0  0  0  0  0  0
   13.4860  -34.4200  -18.0300 O   0  5  0  0  0  0
   15.2430  -33.6100  -20.3240 H   0  0  0  0  0  0
   14.4110  -36.5360  -19.7210 H   0  0  0  0  0  0
   14.6410  -33.9380  -22.3590 H   0  0  0  0  0  0
   15.7620  -38.6220  -23.7550 H   0  0  0  0  0  0
   15.3330  -38.0570  -21.4140 H   0  0  0  0  0  0
   15.1950  -34.6170  -25.8260 H   0  0  0  0  0  0
   15.8806  -37.8889  -25.8363 H   0  0  0  0  0  0
  1  2  2  0  0  0
  1  3  1  0  0  0
  1 13  1  0  0  0
  3  4  2  0  0  0
  3 14  1  0  0  0
  4  5  1  0  0  0
  4 15  1  0  0  0
  5  6  2  0  0  0
  5 10  1  0  0  0
  6  7  1  0  0  0
  6 16  1  0  0  0
  7  8  2  0  0  0
  7 12  1  0  0  0
  8  9  1  0  0  0
  8 11  1  0  0  0
  9 10  2  0  0  0
  9 17  1  0  0  0
 10 18  1  0  0  0
 11 20  1  0  0  0
 12 19  1  0  0  0
M  CHG  1  13  -1
M  END)CTAB";

  bool sanitize = true;
  bool removeHs = false;
  std::unique_ptr<ROMol> refMol{MolBlockToMol(ref, sanitize, removeHs)};
  REQUIRE(refMol);
  std::unique_ptr<ROMol> prbMol{MolBlockToMol(probe, sanitize, removeHs)};
  SECTION("basics") {
    auto solution = shapeit::alignMols(*refMol, *prbMol);
    CHECK(solution.score == Catch::Approx(0.647).margin(0.01));
  }
  SECTION("self") {
    ROMol m2(*refMol);
    auto solution = shapeit::alignMols(*refMol, m2);
    CHECK(solution.score == Catch::Approx(1.00).margin(0.01));
  }
  SECTION("align to volume") {
    GaussianVolume refVolume;
    listAtomVolumes(*refMol, refVolume);
    initOrientation(refVolume);
    auto solution = shapeit::alignMolToVolume(refVolume, *prbMol);

    refVolume.gaussians.clear();
    for (auto &childOverlap : refVolume.childOverlaps) {
      delete childOverlap;
      childOverlap = nullptr;
    }
    refVolume.childOverlaps.clear();
  }
}