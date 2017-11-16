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

/// Sample SNPs
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2017 Anthony J. Greenberg
 * \version 0.5
 *
 * Using the _varfiles_ library to sampple SNPs from a variety of file formats.
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
	string fileType;   // type variant file; Supporting what varfiles supports
	string inFileName; // input file name; will strip off the extension if necessary
	string nSampTxt;   // string version of the number of samples
	uint64_t n;        // number of samples
	
	// set usage message
	string cliHelp = "Command line flags required (in any order):\n  -s NNN (number of samples)\n  -i fileName (input file name, with or without extension)\n  -t fileType (type of SNP file; required if no extension in the input file name)";
	unordered_map <char, string> clInfo;
	parseCL(argc, argv, clInfo);
	auto clIter = clInfo.begin(); // iterator of the command line flags
	
	// Process CLI flags
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
		exit(2);
	} else {
		inFileName = clIter->second;
	}

	clIter = clInfo.find('t');
	if (clIter == clInfo.end()) { // if not there, set to nothing
		fileType = "";
	} else {
		fileType = clIter->second;
	}
	
	// Processing file types first
	if (fileType == "BED") {
		// Doing BED. First figure out if there is an extension on the file name
		string ext;
		bool haveExt = false;
		auto flNit   = inFileName.end() - 1;
		for (; flNit != inFileName.begin(); --flNit) {
			if ((*flNit) == '.') {
				haveExt = true;
				break;
			}
			ext = (*flNit) + ext;
		}
		// get rid of the extension
		if (haveExt) {
			if (ext != "bed") {
				cerr << "WARNING: extension " << ext << " non-standard for BED files" << endl;
			}
			inFileName.erase(flNit, inFileName.end());
		}
		string outFileStub = inFileName + "_s" + nSampTxt;
		BedFileI bedIn(inFileName);
		BedFileO bedOut(outFileStub);
		bedIn.sample(bedOut, n);
		
	} else if (fileType == "TPED"){
		// Doing TPED. First figure out if there is an extension on the file name
		string ext;
		bool haveExt = false;
		auto flNit   = inFileName.end();
		for (; flNit != inFileName.begin(); --flNit) {
			if ((*flNit) == '.') {
				haveExt = true;
				break;
			}
			ext = (*flNit) + ext;
		}
		// get rid of the extension
		if (haveExt) {
			if (ext != "tped") {
				cerr << "WARNING: extension " << ext << " non-standard for TPED files" << endl;
			}
			inFileName.erase(flNit, inFileName.end());
		}
		string outFileStub = inFileName + "_s" + nSampTxt;
		TpedFileI tpedIn(inFileName);
		TpedFileO tpedOut(outFileStub);
		tpedIn.sample(tpedOut, n);

	} else if (fileType == "VCF"){
		// Doing VCF. First figure out if there is an extension on the file name
		string ext;
		bool haveExt = false;
		auto flNit   = inFileName.end() - 1;
		for (; flNit != inFileName.begin(); --flNit) {
			if ((*flNit) == '.') {
				haveExt = true;
				break;
			}
			ext = (*flNit) + ext;
		}
		// get rid of the extension
		if (haveExt) {
			if (ext != "vcf") {
				cerr << "WARNING: extension " << ext << " non-standard for VCF files" << endl;
			}
			inFileName.erase(flNit, inFileName.end());
		}
		string outFileName = inFileName + "_s" + nSampTxt + ".vcf";
		inFileName         = inFileName + ".vcf";
		
		VcfFileI vcfIn(inFileName);
		VcfFileO vcfOut(outFileName);
		vcfIn.sample(vcfOut, n);
		
	} else if (fileType == "HMP"){
		// Doing hapmap. First figure out if there is an extension on the file name
		string ext;
		bool haveExt = false;
		auto flNit   = inFileName.end() - 1;
		for (; flNit != inFileName.begin(); --flNit) {
			if ((*flNit) == '.') {
				if (ext == "txt") { // because typically it is .hmp.txt
					ext = (*flNit) + ext;
					continue;
				}
				haveExt = true;
				break;
			}
			ext = (*flNit) + ext;
		}
		// get rid of the extension
		if (haveExt) {
			if (ext != "hmp.txt") {
				cerr << "WARNING: extension " << ext << " non-standard for HapMap files" << endl;
			}
			inFileName.erase(flNit, inFileName.end());
		}
		string outFileName = inFileName + "_s" + nSampTxt + ".hmp.txt";
		inFileName         = inFileName + ".hmp.txt";
		HmpFileI hmpIn(inFileName);
		HmpFileO hmpOut(outFileName);
		hmpIn.sample(hmpOut, n);
		
	} else { // no (or an unrecognized) file type specified on CL
		// First figure out if there is an extension on the file name
		string ext;
		bool haveExt = false;
		auto flNit   = inFileName.end() - 1;
		for (; flNit != inFileName.begin(); --flNit) {
			if ((*flNit) == '.') {
				if (ext == "txt") { // because typically it is .hmp.txt
					ext = (*flNit) + ext;
					continue;
				}
				haveExt = true;
				break;
			}
			ext = (*flNit) + ext;
		}
		if (haveExt) {
			if (ext == "bed") {
				if (!fileType.empty()) {
					cerr << "WARNING: unrecognized file type " << fileType << " specified. Assuming BED file based on extension" << endl;
				}
				inFileName.erase(flNit, inFileName.end());
				string outFileStub = inFileName + "_s" + nSampTxt;
				BedFileI bedIn(inFileName);
				BedFileO bedOut(outFileStub);
				bedIn.sample(bedOut, n);

			} else if (ext == "tped") {
				if (!fileType.empty()) {
					cerr << "WARNING: unrecognized file type " << fileType << " specified. Assuming TPED file based on extension" << endl;
				}
				inFileName.erase(flNit, inFileName.end());
				string outFileStub = inFileName + "_s" + nSampTxt;
				TpedFileI tpedIn(inFileName);
				TpedFileO tpedOut(outFileStub);
				tpedIn.sample(tpedOut, n);

			} else if (ext == "vcf") {
				if (!fileType.empty()) {
					cerr << "WARNING: unrecognized file type " << fileType << " specified. Assuming VCF file based on extension" << endl;
				}
				inFileName.erase(flNit, inFileName.end());
				string outFileName = inFileName + "_s" + nSampTxt + ".vcf";
				inFileName         = inFileName + "." + ext;
				VcfFileI vcfIn(inFileName);
				VcfFileO vcfOut(outFileName);
				vcfIn.sample(vcfOut, n);

			} else if (ext == "hmp.txt") {
				if (!fileType.empty()) {
					cerr << "WARNING: unrecognized file type " << fileType << " specified. Assuming HapMap file based on extension" << endl;
				}
				inFileName.erase(flNit, inFileName.end());
				string outFileName = inFileName + "_s" + nSampTxt + ".hmp.txt";
				inFileName         = inFileName + "." + ext;
				HmpFileI hmpIn(inFileName);
				HmpFileO hmpOut(outFileName);
				hmpIn.sample(hmpOut, n);

			} else {
				if (fileType.empty()) {
					cerr << "ERROR: urecognized extension (" << ext << ") in input file name " << inFileName << " and no file type specified" << endl;
					exit(3);
				} else {
					cerr << "ERROR: unrecognized extension (" << ext << ") in input file name " << inFileName << " and unrecognized file type (" << fileType  << ") specified" << endl;
					exit(4);
				}

			}
		} else {
			if (fileType.empty()) {
				cerr << "ERROR: no extension in input file name " << inFileName << " and no file type specified" << endl;
				exit(3);
			} else {
				cerr << "ERROR: no extension in input file name " << inFileName << " and unrecognized file type (" << fileType  << ") specified" << endl;
				exit(4);
			}
		}
	}
	
}




