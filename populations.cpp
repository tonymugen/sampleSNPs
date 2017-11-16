/*
 * Copyright (c) 2017 Anthony J. Greenberg
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 * IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/// Connect lines with populations
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2017 Anthony J. Greenberg
 * \version 1.0
 *
 * Implementation of the class that relates individual lines to populations they belong to.
 *
 */

#include <iostream>
#include <fstream>

#include "populations.hpp"

using std::cerr;
using std::endl;
using std::ifstream;

using namespace sampFiles;

/*
 * PopIndex methods
 */

PopIndex::PopIndex(const int *arr, const size_t &N) : _N(N) {
	for (size_t elInd = 0; elInd < N; elInd++) {
		if (arr[elInd] < 0) {
			cerr << "ERROR: population index #" << elInd << " (base-0) is less than zero (" << arr[elInd] << ") in PopIndex array-based initialization." << endl;
			exit(1);
		}
		
		if (arr[elInd] > _popInd.size()) {
			_popInd.resize(arr[elInd]); // maybe this popID is a few steps ahead of the current
			_popInd[arr[elInd] - 1].push_back(elInd);
		}
		else {
			_popInd[arr[elInd] - 1].push_back(elInd);
		}
		
	}
}
PopIndex::PopIndex(const string &inFileName){
	ifstream idxFl(inFileName.c_str());
	
	if (!idxFl) {
		cerr << "ERROR: cannot open file " << inFileName << endl;
		exit(1);
	}
	int tmpIn;
	_N = 0;
	while (idxFl >> tmpIn) {
		if (tmpIn < 0) {
			cerr << "ERROR: negative population IDs not allowed when initializing PopIndex" << endl;
			exit(1);
		}
		else if (tmpIn == 0) {
			cerr << "ERROR: population indexes should be base-1 when initializing PopIndex" << endl;
			exit(2);
		}
		if (tmpIn > _popInd.size()) {
			_popInd.resize(tmpIn); // maybe this popID is a few steps ahead of the current
			_popInd[tmpIn - 1].push_back(_N);
			_N++;
		}
		else {
			_popInd[tmpIn - 1].push_back(_N);
			_N++;
		}
	}
}

