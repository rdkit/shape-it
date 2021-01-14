/*******************************************************************************
parseCommandLine.cpp - Shape-it

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

#include <parseCommandLine.h>

Options parseCommandLine(int argc, char *argv[]) {
  static struct option Arguments[] = {
      {"version", no_argument, NULL, 'v'},
      {"reference", required_argument, NULL, 'r'},
      {"dbase", required_argument, NULL, 'd'},
      {"scores", required_argument, NULL, 's'},
      {"out", required_argument, NULL, 'o'},
      {"format", required_argument, NULL, 'f'},
      {"scoreOnly", no_argument, NULL, 1},
      {"rankBy", required_argument, NULL, 2},
      {"best", required_argument, NULL, 4},
      {"addIterations", required_argument, NULL, 5},
      {"cutoff", required_argument, NULL, 6},
      {"noRef", no_argument, NULL, 11},
      {"help", no_argument, NULL, 'h'}};

  Options o;

  int choice;
  opterr = 0;
  int optionIndex = 0;
  std::string s;

  while ((choice = getopt_long(argc, argv, "vhpr:d:s:o:f:", Arguments,
                               &optionIndex)) != -1) {
    switch (choice) {
    case 'v': //....................................................version
      o.version = true;
      break;

    case 'r': //..................................................reference
      o.refInpFile = optarg;
      o.refInpStream = new std::ifstream(optarg);
      if (!o.refInpStream->good()) {
        mainErr("Error opening input file for reference (-r)");
      }
      break;

    case 'd': //......................................................dbase
      o.dbInpFile = optarg;
      o.dbInpStream = new std::ifstream(optarg);
      if (!o.dbInpStream->good()) {
        mainErr("Error opening input file for database (-d)");
      }
      break;

    case 's': //.....................................................scores
      o.scoreOutFile = optarg;
      o.scoreOutStream = new std::ofstream(optarg);
      if (!o.scoreOutStream->good()) {
        mainErr("Error opening output file for scores (-s)");
      }
      break;

    case 'o': //........................................................out
      o.molOutFile = optarg;
      o.molOutStream = new std::ofstream(optarg);
      if (!o.molOutStream->good()) {
        mainErr("Error opening output file for molecules (-o)");
      }
      break;

    case 'f': //.....................................................format
      o.format = optarg;
#ifdef USE_RDKIT
      if (!o.format.empty() && o.format != "SDF") {
        mainErr("RDKit implementation currently only supports SDF (-f)");
      }
#endif
      break;

    case 1: //....................................................scoreOnly
      o.scoreOnly = true;
      break;

    case 2: //.......................................................rankBy
      s = optarg;
      transform(s.begin(), s.end(), s.begin(), toupper);
      if (s == "TANIMOTO") {
        o.whichScore = tanimoto;
      } else if (s == "TVERSKY_DB") {
        o.whichScore = tversky_db;
      } else if (s == "TVERSKY_REF") {
        o.whichScore = tversky_ref;
      }
      break;

    case 4: //.........................................................best
      o.bestHits = strtol(optarg, NULL, 10);
      break;

    case 5: //................................................addIterations
      o.maxIter = strtol(optarg, NULL, 10);
      break;

    case 6: //.......................................................cutoff
      o.cutOff = strtod(optarg, NULL);
      if (o.cutOff > 1) {
        o.cutOff = 1.0;
      } else if (o.cutOff < 0) {
        o.cutOff = 0.0;
      }
      break;

    case 11: //.......................................................noRef
      o.showRef = false;
      break;

    case 'h': //.......................................................help
      o.help = true;
      break;

    default:
      mainErr("Unknown command line option");
    }
  }

  // If no options are given print the help
  if (optind == 1) {
    o.help = true;
  }

  argc -= optind;
  argv += optind;
  return o;
}
