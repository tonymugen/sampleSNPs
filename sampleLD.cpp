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

/// Sample-based linkage disequilibrium
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2017 Anthony J. Greenberg
 * \version 0.5
 *
 * Using the _varfiles_ library to sampple SNPs calculate pairwise LD. 
 *
 */

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>

#include "varfiles.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::unordered_map;

using namespace sampFiles;

void parseCL(int&, char**, unordered_map<char, string> &);

/** \brief Command line parser
 *
 * Maps flags to values. Flags assumed to be of the form -x.
 *
 * \param[in] argc size of the `argv` array
 * \param[in] argv command line input array
 * \param[out] cli map of tags to values
 */
void parseCL(int &argc, char **argv, unordered_map<char, string> &cli){
	// set to true after encountering a flag token (the character after the dash)
	bool val = false;
	// store the token value here
	char curFlag;
	
	for (int iArg = 1; iArg < argc; iArg++) {
		const char *pchar = argv[iArg];
		
		if (pchar[0] == '-') { // encountered the dash, look for the token after it
			if (!pchar[1]) {
				cerr << "WARNING: forgot character after dash. Ignoring." << endl;
				continue;
			}
			// what follows the dash?
			val     = true;
			curFlag = pchar[1];
			
		} else {
			if (val) {
				val = false;
				cli[curFlag] = pchar;
			} else {
				cerr << "WARNING: command line value " << pchar << " ignored because it is not preceded by a flag" << endl;
			}
			
		}
		
	}
}

int main(int argc, char *argv[]){
	string inFileStub;  // input file name stub (without extension)
	string popFileName; // population index file name
	string nSampTxt;    // string version of the number of samples
	uint64_t n;         // number of samples
	bool pop = false;   // are we stratifying by population?
	
	// set usage message
	string cliHelp = "Command line flags required (in any order):\n  -s NNN (number of samples)\n  -i fileStub (input file name, without extension)\n  -p popIndexFile (optional file with the population index)";
	unordered_map <char, string> clInfo;
	parseCL(argc, argv, clInfo);
	auto clIter = clInfo.begin(); // iterator of the command line flags
	
	// processing CL input
	clIter = clInfo.find('s');
	if (clIter == clInfo.end()) { // if not there, complain
		cerr << "ERROR: specification of the number of samples is required" << endl;
		cerr << cliHelp << endl;
		exit(1);
	} else {
		n        = stoi(clIter->second);
		nSampTxt = clIter->second;
	}
	clIter = clInfo.find('i');
	if (clIter == clInfo.end()) { // if not there, complain
		cerr << "ERROR: specification of the input file name is required" << endl;
		cerr << cliHelp << endl;
		exit(1);
	} else {
		inFileStub = clIter->second;
	}
	
	clIter = clInfo.find('p');
	if (clIter != clInfo.end()) {
		pop         = true;
		popFileName = clIter->second;
	}
	if (pop) { // stratify by population
		PopIndex popIdx(popFileName);
		BedFileI inBed(inFileStub);
		inBed.sampleLD(popIdx, n);
	} else {
		BedFileI inBed(inFileStub);
		inBed.sampleLD(n);
	}
}


