/*******************************************************************************
bestResults.cpp - Shape-it

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

#include <bestResults.h>

BestResults::BestResults(unsigned int n) {
  _bestList.clear();
  _size = n;
  _lowest = 0.0;
  _filled = 0;
}

BestResults::~BestResults(void) {
  std::vector<SolutionInfo *>::iterator it;
  for (it = _bestList.begin(); it != _bestList.end(); ++it) {
    if (*it != NULL) {
      delete *it;
      *it = NULL;
    }
  }
}

bool BestResults::add(SolutionInfo &res) {
  std::vector<SolutionInfo *>::reverse_iterator it;
  if (_filled < _size) {
    SolutionInfo *i = new SolutionInfo(res);
    _bestList.push_back(i);
    ++_filled;
  } else if (res.score < _lowest) {
    return false;
  } else {
    // delete last element
    it = _bestList.rbegin();
    if (*it != NULL) {
      delete *it;
      *it = NULL;
    }

    // make new info element in the list
    *it = new SolutionInfo(res);
  }

  std::sort(_bestList.begin(), _bestList.end(), BestResults::_compInfo());
  it = _bestList.rbegin();
  _lowest = (*it)->score;

  return true;
}

#if 0
void
BestResults::writeMolecules(Options* uo)
{
	if (uo->molOutWriter == NULL)
   {
		return;
   }
	
	std::vector<SolutionInfo* >::iterator it;
	for (it = _bestList.begin(); it != _bestList.end(); ++it)
	{
		if (*it != NULL)
      {
			uo->molOutWriter->Write(&((*it)->dbMol), uo->molOutStream);
      }
	}
	return;
}
#endif

void BestResults::writeScores(Options *uo) {
  if (uo->scoreOutStream == NULL) {
    return;
  }

  std::vector<std::string> score(3);
  std::vector<SolutionInfo *>::iterator it;
  for (it = _bestList.begin(); it != _bestList.end(); ++it) {
    if (*it != NULL) {
      (*it)->printScores(*uo);
    }
  }
  return;
}
