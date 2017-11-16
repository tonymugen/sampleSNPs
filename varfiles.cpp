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

/// Read and write genetic variant files
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2017 Anthony J. Greenberg
 * \version 0.1
 *
 * Implementation of classes that read and write various genetic variant file formats.
 *
 *
 */

#include "varfiles.hpp"
#include "random.hpp"
#include "populations.hpp"

#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <sstream>
#include <cstdint>
#include <cmath>
#include <limits>

using std::fstream;
using std::ofstream;
using std::stringstream;
using std::string;
using std::vector;
using std::unordered_map;
using std::cerr;
using std::endl;
using std::flush;
using std::system_error;
using std::ios;
using std::bad_alloc;
using std::streamsize;

using namespace sampFiles;

// GbinFile methods

uint64_t GbinFileI::_numLines(){
    
    if (_varFile.is_open()) {
        _varFile.close();
    }
    try {
        _varFile.open(_fileName.c_str(), ios::in | ios::binary | ios::ate);
    } catch (system_error &error) {
        cerr << "ERROR: cannot open binary file " << _fileName << " for input: " << error.code().message() << flush;
        perror(" ");
        exit(1);
        
    }
    size_t Ntot = _varFile.tellg();
    _varFile.close();
    if (Ntot % (_elemSize*_nCols)) {
        throw string("Number of elements not divisible by row size");
    }
    return Ntot/(_elemSize*_nCols);
}

void GbinFile::close(){
    if (_varFile.is_open()) {
        _varFile.close();
    }
}


void GbinFileI::open(){
    if (_varFile.is_open()) {
        _varFile.close();
    }
    try {
        _varFile.open(_fileName.c_str(), ios::in | ios::binary);
    } catch (system_error &error) {
        cerr << "ERROR: cannot open binary file " << _fileName << " for input: " << error.code().message() << flush;
        perror(" ");
        exit(1);
    }
}

void GbinFileI::sample(GbinFileO &out, const uint64_t &n){
    
    if (n == 0) {
        cerr << "WARNING: zero SNPs requested. Nothing to be done." << endl;
        return;
    }
    // start by figuring out the number of rows
    
    
    // calculate number of rows in the input binary file; file is closed before the function exits.
    uint64_t N = _numLines();
    
    // Test for potential problems
    if (N < n) {
        cerr << "ERROR: requested a sample of " << n << " rows that is greater than the number of rows (" << N << ") in the input file." << endl;
        return;
    } else if (N == n) {
        cerr << "WARNING: sample size (" << n << ") the same as the number of rows (" << N << ") in the input file. Simply copying the file." << endl;
        char *buf;
        const size_t bufSize = BUF_SIZE;
        try {
            buf = new char[bufSize];
        } catch (bad_alloc&) {
            cerr << "ERROR: failed to allocate buffer" << endl;
            exit(4);
        }
        
        try {
            _varFile.open(_fileName.c_str(), ios::in | ios::binary);
        } catch (system_error &error) {
            cerr << "ERROR: cannot open binary file " << _fileName << " for input: " << error.code().message() << flush;
            perror(" ");
            exit(1);
        }
        try {
            out._varFile.open(out._fileName.c_str(), ios::out | ios::binary | ios::trunc);
        } catch (system_error &error) {
            cerr << "ERROR: cannot open binary file " << out._fileName << " for output: " << error.code().message() << flush;
            perror(" ");
            exit(1);
        }
        
        while (_varFile) {
            _varFile.read(buf, bufSize);
            out._varFile.write(buf, _varFile.gcount()); // must be gcount() in case we got to the end of the file
        }
        _varFile.close();
        out._varFile.close();
        delete [] buf;
        return;
    }
    
    // Passed all the tests, proceed to sampling
    char *ROWbuf;
    size_t rowSize = _nCols*_elemSize;
    try {
        ROWbuf = new char[rowSize]; // .read() does not append a null-character
    } catch (bad_alloc &) {
        cerr << "ERROR: failed to allocate the row buffer" << endl;
        exit(4);
    }
    
    try {
        _varFile.open(_fileName.c_str(), ios::in | ios::binary);
    } catch (system_error &error) {
        cerr << "ERROR: cannot open binary file " << _fileName << " for input: " << error.code().message() << flush;
        perror(" ");
        exit(1);
    }
    
    if (out._varFile.is_open()) {
        out.close();
    }
    try {
        out._varFile.open(out._fileName.c_str(), ios::out | ios::binary | ios::trunc);
    } catch (system_error &error) {
        cerr << "ERROR: cannot open binary file " << out._fileName << " for output: " << error.code().message() << flush;
        perror(" ");
        exit(1);
    }
    
    uint64_t nloc = n; // local copy of n
    uint64_t S;
    uint64_t cumS = 0; // cumulative position in the file
    RanDraw *rowSamp;  // so that I can catch RanDraw constructor exceptions
    try {
        rowSamp = new RanDraw();
    } catch (string error) {
        cerr << "ERROR: " << error << endl;
        exit(5);
    }
    while (nloc) {
        try {
            S = rowSamp->vitter(nloc, N); // sample the number of rows to skip; keep track of the running total
        } catch (string error) {
            cerr << "ERROR: " << error << endl;
            exit(5);
        }
        cumS += S;
        N    -= S + 1;
        _varFile.seekg(cumS*rowSize); // seekg() index is base-0
        _varFile.read(ROWbuf, rowSize);
        cumS++;                        // step up cumS because seekg() will not change it (unlike the getline() that automatially advances)
        out._varFile.write(ROWbuf, rowSize);
        nloc--;
        
    }
    _varFile.close();
    out._varFile.close();
    delete [] ROWbuf;
    delete rowSamp;
    
}

void GbinFileO::open(){
    if (_varFile.is_open()) {
        _varFile.close();
    }
    try {
        _varFile.open(_fileName.c_str(), ios::out | ios::binary);
    } catch (system_error &error) {
        cerr << "ERROR: cannot open binary file " << _fileName << " for output: " << error.code().message() << flush;
        perror(" ");
        exit(1);
        
    }
}

// BedFile methods
BedFile::BedFile() : GbinFile() {
    _bimFile.exceptions(fstream::badbit); // for some reason, setting failbit | badbit throws errors no matter what
	_famFile.exceptions(fstream::badbit);
    // _elemSize is set to sizeof(char) in the GbinFile default constructor
}
BedFile::BedFile(const string &stubName) : GbinFile(), _fileStub{stubName} {
	_fileName = _fileStub + ".bed";
    _bimFile.exceptions(fstream::badbit);
	_famFile.exceptions(fstream::badbit);
    // _elemSize is set to sizeof(char) in the GbinFile default constructor
}

const vector<char> BedFile::_masks = {static_cast<char>(0x03), static_cast<char>(0x0C), static_cast<char>(0x30), static_cast<char>(0xC0)};
// M is missing, H is heterozygous, homozygous derived is all 0x00, homozygous ancestral is never tested
const unordered_map<char, string> BedFile::_tests = {
    {'M', {static_cast<char>(0x01), static_cast<char>(0x04), static_cast<char>(0x10), static_cast<char>(0x40)}},
    {'H', {static_cast<char>(0x02), static_cast<char>(0x08), static_cast<char>(0x20), static_cast<char>(0x80)}}
};

BedFile::~BedFile(){
    if (_bimFile.is_open()) {
        _bimFile.close();
    }
    if (_famFile.is_open()) {
        _famFile.close();
    }
}

void BedFile::close(){
    if (_varFile.is_open()) {
        _varFile.close();
    }
    if (_bimFile.is_open()) {
        _bimFile.close();
    }
    if (_famFile.is_open()) {
        _famFile.close();
    }
}

uint64_t BedFileI::_numLines(){
    if (!_nCols) {
        _nCols = _famLines();
    }
    if (_varFile.is_open()) {
        _varFile.close();
    }
    try {
        _varFile.open(_fileName.c_str(), ios::in | ios::binary | ios::ate);
    } catch (system_error &error) {
        cerr << "ERROR: cannot open binary file " << _fileName << " for output: " << error.code().message() << flush;
        perror(" ");
        exit(1);
        
    }
    uint64_t Nbed = (_nCols/4UL) + static_cast<uint64_t>( (_nCols%4UL) > 0UL );
    uint64_t N = static_cast<uint64_t>(_varFile.tellg()) - 3UL; // 3 is the number of magic bytes
    if (N % Nbed) {
        cerr << "ERROR: total file size " << N << " not divisible by the number of compressed colums " << Nbed << " in BedFile _numLines()" << endl;
        exit(7);
    }
    _varFile.close();
    
    return N/Nbed;
}

uint64_t BedFileI::_famLines(){
	uint64_t N = 0;
	
	/*
	 * I am using the line-end counting method. It is > 2-fold faster than reading lines with getline().
	 */
	if (_famFile.is_open()) {
        _famFile.close();
    }
    string famName = _fileStub + ".fam";
    try {
        _famFile.open(famName.c_str(), ios::in);
    } catch (system_error &error) {
        cerr << "ERROR: cannot open .fam file " << famName << " for input: " << error.code().message() << flush;
        perror(" ");
        exit(1);
        
    }
    
	const size_t bufSize = BUF_SIZE;
	char *buf;
	try {
		buf = new char[bufSize];
	} catch (bad_alloc& error) {
		cerr << "ERROR: failed to allocate buffer in BedFile::_famLines(): " << error.what() << endl;
		exit(4);
	}
	
	while (_famFile) {
		_famFile.read(buf, bufSize);
		for (size_t i = 0; i < _famFile.gcount(); i++) {
			if (buf[i] == '\n') {
				N++;
			}
		}
	}
	delete [] buf;
	_famFile.close();
	
	return N;
}

uint64_t BedFileI::_famLines(fstream &fam){
	uint64_t N = 0;
	
	/*
	 * I am using the line-end counting method. It is > 2-fold faster than reading lines with getline().
	 */
    if (_famFile.is_open()) {
        _famFile.close();
    }
    string famName = _fileStub + ".fam";
    try {
        _famFile.open(famName.c_str(), ios::in);
    } catch (system_error &error) {
        cerr << "ERROR: cannot open .fam file " << famName << " for input: " << error.code().message() << flush;
        perror(" ");
        exit(1);
        
    }
	const size_t bufSize = BUF_SIZE;
	char *buf;
	try {
		buf = new char[bufSize];
	} catch (bad_alloc& error) {
		cerr << "ERROR: failed to allocate buffer in BedFile::_famLines(): " << error.what() << endl;
		exit(4);
	}
    
    if (!fam.is_open()) {
        throw string("Output .fam filestream not open");
    }
	while (_famFile) {
		_famFile.read(buf, bufSize);
		for (size_t i = 0; i < _famFile.gcount(); i++) {
			if (buf[i] == '\n') {
				N++;
			}
		}
		fam.write(buf, _famFile.gcount());
	}
	delete [] buf;
	_famFile.close();
	fam.close();
	
	return N;
}

void BedFileI::_ld(const char *snp1, const char *snp2, const size_t &N, const unsigned short &pad, double &rSq, double &Dprime, double &dcnt1, double &dcnt2){
	// Using the Gaunt et al. (2007) formulation
	// The data in .bed files are strictly biallelic by definition of the format, so that makes things easier.
	// Using derived alleles because they are mostly at low frequency so the expected number of calculations is low
	// However, then convert everything to the most common allele frequencies to make things compatible with Gaunt et al.
	// Using the notation in that paper, too, to make things more clear
	
	// constants
	const char allMiss = static_cast<char>(0x55); // all genotypes missing in a byte
	const char homDerv = static_cast<char>(0x00); // homozygous derived for all bit-pair positions
	
	// modifiable parameters
	double D     = 0.0;
	double p1    = 0.0; // frequency of the derived allele at snp1
	double q1    = 0.0; // frequency of the derived allele at snp2
	double f11   = 0.0; // frequency of the (hom. maj. allele @1, hom. maj.allel @2) haplotype; to be estimated; D = f11 - p1*q1
	double n11   = 0.0; // # of 1,1 homozygote haplotypes
	double n12   = 0.0; // # of hom. 1, heterozygous 2 haplotypes
	double n21   = 0.0; // # of het. 1, hom. 2 haplotypes
	double n22   = 0.0; // # of het. 1, het. 2 haplotypes
	double Npres = 0.0; // # of present genotypes (present at both loci)
	double p1q1;        // p1*q1
	double p2q2;        // p2*q2
	// f11 min and max
	double f11Min;
	double f11Max;
	// re-used intermediate variables
	double opq;         // 1-2p-2q
	double smn;         // 2n11+n12+n21
	double ppq;         // p1+q1
	// qubic equation coefficients
	double a;
	double b;
	double c;
	double d;
	// decision variables
	double bda;    // b/a
	double xN;
	double deltaSq;
	double hSq;
	double gammaN;
	double Delta3;
	
	char curBitPair1;
	char curBitPair2;
	
	// Data missing at either locus are ignored (i.e. I am using pairwise-complete observations)
	for (size_t ind = 0; ind < N-1; ind++) { // leave the last byte out for now since it has padding bits (in general)
		if ( (snp1[ind] == allMiss) || (snp2[ind] == allMiss) ) { // a byte with all data missing; ignore (NOTE the OR)
			continue;
		}
		
		// go through the byte from the end
		for (size_t inByte = 0; inByte < 4; inByte++) {
			curBitPair1 = snp1[ind] & _masks[inByte];
			curBitPair2 = snp2[ind] & _masks[inByte];
			
			// Test missingness
			if ( (curBitPair1 == _tests.at('M')[inByte]) || (curBitPair2 == _tests.at('M')[inByte]) ) { // must use .at() because _tests is const
				continue;
			}
			Npres += 1.0;
			// Go through all the pairwise possibilities; adding up gametes (one per diploid)
			if (curBitPair1 == homDerv) {
				p1 += 1.0;
				if (curBitPair2 == homDerv) {
					q1  += 1.0;
					n11 += 1.0;
				} else if (curBitPair2 == _tests.at('H')[inByte]){
					q1  += 0.5;
					n12 += 1.0;
				}
			} else if (curBitPair1 == _tests.at('H')[inByte]){
				p1 += 0.5;
				if (curBitPair2 == homDerv) {
					q1  += 1.0;
					n21 += 1.0;
				} else if (curBitPair2 == _tests.at('H')[inByte]){
					q1  += 0.5;
					n22 += 1.0;
				}
			} else {
				if (curBitPair2 == homDerv) {
					q1 += 1.0;
				} else if (curBitPair2 == _tests.at('H')[inByte]){
					q1 += 0.5;
				}
			}
		}
	}
	// finally, deal with the padded byte
	for (size_t inByte = 0; inByte < (4 - pad); inByte++) { // not testing for sane values of pad, since this is a protected function and I can expect sane values in correct class implementation
		curBitPair1 = snp1[N-1] & _masks[inByte];
		curBitPair2 = snp2[N-1] & _masks[inByte];
		
		// Test missingness
		if ( (curBitPair1 == _tests.at('M')[inByte]) || (curBitPair2 == _tests.at('M')[inByte]) ) { // must use .at() because _tests is const
			continue;
		}
		Npres += 1.0;
		// Go through all the pairwise possibilities; adding up gametes (one per diploid)
		if (curBitPair1 == homDerv) {
			p1 += 1.0;
			if (curBitPair2 == homDerv) {
				q1  += 1.0;
				n11 += 1.0;
			} else if (curBitPair2 == _tests.at('H')[inByte]){
				q1  += 0.5;
				n12 += 1.0;
			}
		} else if (curBitPair1 == _tests.at('H')[inByte]){
			p1 += 0.5;
			if (curBitPair2 == homDerv) {
				q1  += 1.0;
				n21 += 1.0;
			} else if (curBitPair2 == _tests.at('H')[inByte]){
				q1  += 0.5;
				n22 += 1.0;
			}
		} else {
			if (curBitPair2 == homDerv) {
				q1 += 1.0;
			} else if (curBitPair2 == _tests.at('H')[inByte]){
				q1 += 0.5;
			}
		}
	}
	
	// minor allele counts
	dcnt1 = (p1 <= Npres - p1 ? 2.0*p1 : 2.0*(Npres - p1));
	dcnt2 = (q1 <= Npres - q1 ? 2.0*q1 : 2.0*(Npres - q1));
	
	if (Npres <= 1.0) {
		rSq    = -9.0;
		Dprime = -9.0;
		return;
	}
	if ( (dcnt1 == 0.0) || (dcnt2 == 0.0) ) { // even if the SNPs are pre-screened for polymorphism, the alternatives may correspond to missing at the other locus
		rSq    = -9.0;
		Dprime = -9.0;
		return;
	}
	if (n22 == Npres) { // if everyone is het, cannot determine LD because we don't know the phase
		rSq    = -9.0;
		Dprime = -9.0;
		return;
	}
	// Finish calculating statistics if everything is sane
	p1   = p1/Npres;
	q1   = q1/Npres;
	// intermediate values
	p1q1 = p1*q1;

	if (n22 == 0.0) { // no het/het match-ups; means I can directly enumerate haplotypes
		f11  = (n11 + 0.5*(n12 + n21))/Npres;
		D    = f11 - p1q1;
		p2q2 = (1.0-p1)*(1.0-q1);
		rSq  = (D*D)/(p1q1*p2q2);
		
		if (D < -10.0*EPS) {
			double Dmax = ( p1q1 <= p2q2 ? -p1q1 : -p2q2 ); // (9) of Gaunt et al.
			Dprime      = D/Dmax;
		} else if (D > 10.0*EPS) {
			p1q1 = p1*(1.0 - q1);
			p2q2 = (1.0 - p1)*q1;
			
			double Dmax = ( p1q1 <= p2q2 ? p1q1 : p2q2 ); // (9) of Gaunt et al.
			Dprime      = D/Dmax;
		} else {
			Dprime = 0.0;
		}

		return;
	}
	
	// There is n22, so need to estimate f11 via ML (using Gaunt et al. cubic equation); they say f11 is between major alleles, but that is not necessary. Results will be the same for minor or any other combination
	// Calculate the bracketing values for f11:
	f11Min = n11 + 0.5*(n12 + n21); // i.e., all the het/het haplotypes are crossovers
	f11Max = f11Min + 0.5*n22;      // i.e., all het/het haplotypes are in phase
	f11Min = f11Min/Npres - 100.0*EPS; // some padding to take into account round-off errors
	f11Max = f11Max/Npres + 100.0*EPS;
	// re-used auxiliary variables
	ppq = p1 + q1;
	opq = 1.0 - 2.0*ppq;
	smn = 2.0*n11 + n12 + n21;
	// coefficients
	a   = 2.0*Npres; // will need to multiply again later to make 4N; doing this to re-use the 2N
	b   = a*opq - 2.0*smn - n22;
	c   = a*p1q1 - smn*opq - n22*(1.0 - ppq);
	d   = -smn*p1q1;
	a  *= 2.0;
	// decision variables
	bda     = b/a;
	xN      = -bda/3.0;
	deltaSq = bda*bda/9.0 - c/(3.0*a);
	hSq     = 4.0*a*a*deltaSq*deltaSq*deltaSq;
	gammaN  = xN*( xN*( a*xN + b ) + c ) + d;
	Delta3  = gammaN*gammaN - hSq;
	
	// now decide how many roots we have
	if (Delta3 > 100.0*EPS) { // everything is cool, only one distinct root
		f11 = xN + cbrt((sqrt(Delta3) - gammaN)/(2.0*a)) + cbrt(-(sqrt(Delta3) + gammaN)/(2.0*a));
	} else if (Delta3 < -100.0*EPS){ // Three roots; worst case scenario
		double theta  = acos(-gammaN/sqrt(hSq))/3.0;
		double ddelta = 2.0*sqrt(deltaSq);
		double alpha  = xN + ddelta*cos(theta);
		double beta   = xN + ddelta*cos(2.0*PI/3.0 + theta);
		double gamma  = xN + ddelta*cos(4.0*PI/3.0 + theta);
		
		if ( (alpha < f11Min) || (alpha > f11Max) ) {
			if ( (beta < f11Min) || (beta > f11Max) ) {
				if ( (gamma < f11Min) || (gamma > f11Max) ) {
					rSq    = -9.0;
					Dprime = -9.0;
					
					return;

				} else {
					f11 = gamma;
				}
			} else if ( (gamma < f11Min) || (gamma > f11Max) ) {
				f11 = beta;
			} else {
				// two plausible roots; choose the one with smallest |D|
				f11 = (fabs(beta - p1q1) <= fabs(gamma - p1q1) ? beta : gamma);
			}
		} else {
			if ( (beta < f11Min) || (beta > f11Max) ) {
				if ( (gamma < f11Min) || (gamma > f11Max) ) {
					f11 = alpha;
				} else {
					// two plausible roots; choose the one with smallest |D|
					f11 = (fabs(alpha - p1q1) <= fabs(gamma - p1q1) ? alpha : gamma);
				}
			} else if ( (gamma < f11Min) || (gamma > f11Max) ) {
				// two plausible roots; choose the one with smallest |D|
				f11 = (fabs(beta - p1q1) <= fabs(alpha - p1q1) ? beta : alpha);
			} else {
				// three plausible roots; choose the one with smallest |D|
				double minAB = (fabs(beta - p1q1) <= fabs(alpha - p1q1) ? beta : alpha);
				f11 = (fabs(minAB - p1q1) <= fabs(gamma - p1q1) ? minAB : gamma);
			}

		}
	} else { // Delta3 == 0; two different roots
		double mu    = cbrt(gammaN/(2.0*a));
		double alpha = xN + mu;
		double gamma = xN - 2.0*mu;
		if ( (alpha < f11Min) || (alpha > f11Max) ) { // f11 cannot be greater than min(p1, q1)
			if ( (gamma < f11Min) || (gamma > f11Max) ) {
				rSq    = -9.0;
				Dprime = -9.0;
				
				return;
			} else {
				f11 = gamma;
			}
		} else {
			if ( (gamma < f11Min) || (gamma > f11Max) ) { // impossible gamma; going with alpha
				f11 = alpha;
			} else { // both doable; conservatively choose the one minimzing |D|
				f11 = (fabs(alpha - p1q1) <= fabs(gamma - p1q1) ? alpha : gamma);
			}
		}
		
	}
	
	D    = f11 - p1q1;
	p2q2 = (1.0-p1)*(1.0-q1);
	rSq  = (D*D)/(p1q1*p2q2);
	
	if (D < -10.0*EPS) {
		double Dmax = ( p1q1 <= p2q2 ? -p1q1 : -p2q2 ); // (9) of Gaunt et al.
		Dprime      = D/Dmax;
	} else if (D > 10.0*EPS) {
		p1q1 = p1*(1.0 - q1);
		p2q2 = (1.0 - p1)*q1;
		
		double Dmax = ( p1q1 <= p2q2 ? p1q1 : p2q2 ); // (9) of Gaunt et al.
		Dprime      = D/Dmax;
	} else {
		Dprime = 0.0;
	}
	
}

void BedFileI::_ld(const char *snp1, const char *snp2, const PopIndex &popID, vector<double> &rSq, vector<double> &Dprime, vector<double> &dcnt1, vector<double> &dcnt2){
	
	// The data in .bed files are strictly biallelic by definition of the format, so that makes things easier.
	// Using derived alleles because they are mostly at low frequency so the expected number of calculations is low
	
	// constants
	const char homDerv = static_cast<char>(0x00); // homozygous derived for all bit-pair positions
	
	// modifiable parameters
	double D;
	double p1;     // frequency of the derived allele at snp1
	double q1;     // frequency of the derived allele at snp2
	double f11;    // frequency of the (hom. maj. allele @1, hom. maj.allel @2) haplotype; to be estimated; D = f11 - p1*q1
	double n11;    // # of 1,1 homozygote haplotypes
	double n12;    // # of hom. 1, heterozygous 2 haplotypes
	double n21;    // # of het. 1, hom. 2 haplotypes
	double n22;    // # of het. 1, het. 2 haplotypes
	double Npres;  // # of present genotypes (present at both loci)
	double p1q1;   // p1*q1
	double p2q2;   // p2*q2
	// f11 min and max
	double f11Min;
	double f11Max;
	// re-used intermediate variables
	double opq;       // 1-2p-2q
	double smn;       // 2n11+n12+n21
	double ppq;       // p1+q1
	// qubic equation coefficients
	double a;
	double b;
	double c;
	double d;
	// decision variables
	double bda;    // b/a
	double xN;
	double deltaSq;
	double hSq;
	double gammaN;
	double Delta3;
	
	char curBitPair1;
	char curBitPair2;

	for (size_t iPop = 0; iPop < popID.popNumber(); iPop++) {
		// reset everything
		p1    = 0.0;
		q1    = 0.0;
		f11   = 0.0;
		n11   = 0.0;
		n12   = 0.0;
		n21   = 0.0;
		n22   = 0.0;
		Npres = 0.0;

		for (auto popIt = popID[iPop].begin(); popIt != popID[iPop].end(); ++popIt) {
			size_t byteInd    = (*popIt)/4;
			size_t bitPairInd = (*popIt) % 4; // will automatically index from the end of the byte
			curBitPair1 = snp1[byteInd] & _masks[bitPairInd];
			curBitPair2 = snp2[byteInd] & _masks[bitPairInd];
			
			// Test miss
			if ( (curBitPair1 == _tests.at('M')[bitPairInd]) || (curBitPair2 == _tests.at('M')[bitPairInd]) ) { // must use .at() because _tests is const
				continue;
			}
			Npres += 1.0;
			// Go through all the pairwise possibilities; adding up gametes (one per diploid)
			if (curBitPair1 == homDerv) {
				p1 += 1.0;
				if (curBitPair2 == homDerv) {
					q1  += 1.0;
					n11 += 1.0;
				} else if (curBitPair2 == _tests.at('H')[bitPairInd]){
					q1  += 0.5;
					n12 += 1.0;
				}
			} else if (curBitPair1 == _tests.at('H')[bitPairInd]){
				p1 += 0.5;
				if (curBitPair2 == homDerv) {
					q1  += 1.0;
					n21 += 1.0;
				} else if (curBitPair2 == _tests.at('H')[bitPairInd]){
					q1  += 0.5;
					n22 += 1.0;
				}
			} else {
				if (curBitPair2 == homDerv) {
					q1 += 1.0;
				} else if (curBitPair2 == _tests.at('H')[bitPairInd]){
					q1 += 0.5;
				}
			}
		}
		dcnt1[iPop]  = (p1 <= Npres - p1 ? 2.0*p1 : 2.0*(Npres - p1));
		dcnt2[iPop]  = (q1 <= Npres - q1 ? 2.0*q1 : 2.0*(Npres - q1));
		
		if (Npres <= 1.0) {
			rSq[iPop]    = -9.0;
			Dprime[iPop] = -9.0;
			continue;
		}
		if ( (dcnt1[iPop] == 0.0) || (dcnt2[iPop] == 0.0) ) { // even if the SNPs are pre-screened for polymorphism, the alternatives may correspond to missing at the other locus
			rSq[iPop]    = -9.0;
			Dprime[iPop] = -9.0;
			continue;
		}
		if (n22 == Npres) { // if everyone is het, cannot determine LD because we don't know the phase
			rSq[iPop]    = -9.0;
			Dprime[iPop] = -9.0;
			continue;
		}
		// Finish calculating statistics if everything is sane
		p1   = p1/Npres;
		q1   = q1/Npres;
		// intermediate values
		p1q1 = p1*q1;
		
		if (n22 == 0.0) { // no het/het match-ups; means I can directly enumerate haplotypes
			f11        = (n11 + 0.5*(n12 + n21))/Npres;
			D          = f11 - p1q1;
			p2q2       = (1.0-p1)*(1.0-q1);
			rSq[iPop]  = (D*D)/(p1q1*p2q2);
			
			if (D < -10.0*EPS) {
				double Dmax  = ( p1q1 <= p2q2 ? -p1q1 : -p2q2 ); // (9) of Gaunt et al.
				Dprime[iPop] = D/Dmax;
			} else if (D > 10.0*EPS) {
				p1q1 = p1*(1.0 - q1);
				p2q2 = (1.0 - p1)*q1;
				
				double Dmax  = ( p1q1 <= p2q2 ? p1q1 : p2q2 ); // (9) of Gaunt et al.
				Dprime[iPop] = D/Dmax;
			} else {
				Dprime[iPop] = 0.0;
			}
			
		} else {
			// There is n22, so need to estimate f11 via ML (using Gaunt et al. cubic equation); they say f11 is between major alleles, but that is not necessary. Results will be the same for minor or any other combination
			// Calculate the bracketing values for f11:
			f11Min = n11 + 0.5*(n12 + n21); // i.e., all the het/het haplotypes are crossovers
			f11Max = f11Min + 0.5*n22;      // i.e., all het/het haplotypes are in phase
			f11Min = f11Min/Npres - 100.0*EPS; // some padding to take into account round-off errors
			f11Max = f11Max/Npres + 100.0*EPS;
			// re-used auxiliary variables
			ppq = p1 + q1;
			opq = 1.0 - 2.0*ppq;
			smn = 2.0*n11 + n12 + n21;
			// coefficients
			a   = 2.0*Npres; // will need to multiply again later to make 4N; doing this to re-use the 2N
			b   = a*opq - 2.0*smn - n22;
			c   = a*p1q1 - smn*opq - n22*(1.0 - ppq);
			d   = -smn*p1q1;
			a  *= 2.0;
			// decision variables
			bda     = b/a;
			xN      = -bda/3.0;
			deltaSq = bda*bda/9.0 - c/(3.0*a);
			hSq     = 4.0*a*a*deltaSq*deltaSq*deltaSq;
			gammaN  = xN*( xN*( a*xN + b ) + c ) + d;
			Delta3  = gammaN*gammaN - hSq;
			
			// now decide how many roots we have
			if (Delta3 > 100.0*EPS) { // everything is cool, only one distinct root
				f11 = xN + cbrt((sqrt(Delta3) - gammaN)/(2.0*a)) + cbrt(-(sqrt(Delta3) + gammaN)/(2.0*a));
			} else if (Delta3 < -100.0*EPS){ // Three roots; worst case scenario
				double theta  = acos(-gammaN/sqrt(hSq))/3.0;
				double ddelta = 2.0*sqrt(deltaSq);
				double alpha  = xN + ddelta*cos(theta);
				double beta   = xN + ddelta*cos(2.0*PI/3.0 + theta);
				double gamma  = xN + ddelta*cos(4.0*PI/3.0 + theta);
				
				if ( (alpha < f11Min) || (alpha > f11Max) ) {
					if ( (beta < f11Min) || (beta > f11Max) ) {
						if ( (gamma < f11Min) || (gamma > f11Max) ) {
							rSq[iPop]    = -9.0;
							Dprime[iPop] = -9.0;
							continue;
						} else {
							f11 = gamma;
						}
					} else if ( (gamma < f11Min) || (gamma > f11Max) ) {
						f11 = beta;
					} else {
						// two plausible roots; choose the one with smallest |D|
						f11 = (fabs(beta - p1q1) <= fabs(gamma - p1q1) ? beta : gamma);
					}
				} else {
					if ( (beta < f11Min) || (beta > f11Max) ) {
						if ( (gamma < f11Min) || (gamma > f11Max) ) {
							f11 = alpha;
						} else {
							// two plausible roots; choose the one with smallest |D|
							f11 = (fabs(alpha - p1q1) <= fabs(gamma - p1q1) ? alpha : gamma);
						}
					} else if ( (gamma < f11Min) || (gamma > f11Max) ) {
						// two plausible roots; choose the one with smallest |D|
						f11 = (fabs(beta - p1q1) <= fabs(alpha - p1q1) ? beta : alpha);
					} else {
						// three plausible roots; choose the one with smallest |D|
						double minAB = (fabs(beta - p1q1) <= fabs(alpha - p1q1) ? beta : alpha);
						f11 = (fabs(minAB - p1q1) <= fabs(gamma - p1q1) ? minAB : gamma);
					}
					
				}
			} else { // Delta3 == 0; two different roots
				double mu    = cbrt(gammaN/(2.0*a));
				double alpha = xN + mu;
				double gamma = xN - 2.0*mu;
				if ( (alpha < f11Min) || (alpha > f11Max) ) { // f11 cannot be greater than min(p1, q1)
					if ( (gamma < f11Min) || (gamma > f11Max) ) {
						rSq[iPop]    = -9.0;
						Dprime[iPop] = -9.0;
						continue;
					} else {
						f11 = gamma;
					}
				} else {
					if ( (gamma < f11Min) || (gamma > f11Max) ) { // impossible gamma; going with alpha
						f11 = alpha;
					} else { // both doable; conservatively choose the one minimzing |D|
						f11 = (fabs(alpha - p1q1) <= fabs(gamma - p1q1) ? alpha : gamma);
					}
				}
				
			}
			
			D         = f11 - p1q1;
			p2q2      = (1.0-p1)*(1.0-q1);
			rSq[iPop] = (D*D)/(p1q1*p2q2);
			
			if (D < -10.0*EPS) {
				double Dmax  = ( p1q1 <= p2q2 ? -p1q1 : -p2q2 ); // (9) of Gaunt et al.
				Dprime[iPop] = D/Dmax;
			} else if (D > 10.0*EPS) {
				p1q1 = p1*(1.0 - q1);
				p2q2 = (1.0 - p1)*q1;
				
				double Dmax  = ( p1q1 <= p2q2 ? p1q1 : p2q2 ); // (9) of Gaunt et al.
				Dprime[iPop] = D/Dmax;
			} else {
				Dprime[iPop] = 0.0;
			}

		}
	}
}

void BedFileI::open(){
	string bimName = _fileStub + ".bim";
	string famName = _fileStub + ".fam";
	
	
	try {
		_varFile.open(_fileName.c_str(), ios::in | ios::binary);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open BED file " << _fileName << " for input: " << error.code().message() << flush;
		perror(" ");
		exit(1);

	}
	// check the magic bytes and SNP-major status
	char magic[3];
	_varFile.read(magic, 3);
	if ( (magic[0] != static_cast<char>(0x6C)) || (magic[1] != static_cast<char>(0x1B)) ) {
		cerr << "ERROR: binary file " << _fileName << " not recognized as BED" << endl;
		exit(2);
	} else if (magic[2] != static_cast<char>(0x01)){
		cerr << "ERROR: binary file " << _fileName << " not SNP-major. Run a newer version of plink to fix it." << endl;
		exit(3);
	}
	
	try {
		_bimFile.open(bimName.c_str(), ios::in);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open .bim file " << bimName << " for input: " << error.code().message() << flush;
		perror(" ");
		exit(1);

	}
	
	try {
		_famFile.open(famName.c_str(), ios::in);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open .fam file " << famName << " for input: " << error.code().message() << flush;
		perror(" ");
		exit(1);

	}
	
}

void BedFileI::sample(BedFileO &out, const uint64_t &n){
	
	if (n == 0) {
		cerr << "WARNING: zero SNPs requested. Nothing to be done." << endl;
		return;
	}
	// start by figuring out the number of SNPs in the .bed file
	if (out._varFile.is_open()) {
		out.close();
	}
	// First have to know the number of lines
	string outFam = out._fileStub + ".fam";
	
	try {
		_famFile.open(_fileName.c_str(), ios::in); // _fileName has the .bed file name, with the extension
	} catch (system_error &error) {
		cerr << "ERROR: cannot open .fam file " << _fileName << " for input: " << error.code().message() << flush;
		perror(" ");
		exit(1);
	}
	
	try {
		out._famFile.open(outFam.c_str(), ios::out | ios::trunc);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open .fam file " << out._fileName << " for output: " << error.code().message() << flush;
		perror(" ");
		exit(1);
	}
    
    uint64_t N    = _numLines(); // number of SNPs in the .bed file
    uint64_t Nbed = (_nCols/4UL) + static_cast<uint64_t>( (_nCols%4UL) > 0UL );
	_varFile.close();
	
	// Test for potential problems
	string inBim  = _fileStub + ".bim";
	string outBim = out._fileStub + ".bim";
	if (N < n) {
		cerr << "ERROR: requested a sample of " << n << " SNPs that is greater than the number of SNPs (" << N << ") in the input file." << endl;
		return;
	} else if (N == n) {
		cerr << "WARNING: sample size (" << n << ") the same as the number of SNPs (" << N << ") in the input file. Simply copying the files." << endl;
		char *buf;
		const size_t bufSize = BUF_SIZE;
		try {
			buf = new char[bufSize];
		} catch (bad_alloc&) {
			cerr << "ERROR: failed to allocate buffer" << endl;
			exit(4);
		}
		
		// Copy .bed
		try {
			_varFile.open(_fileName.c_str(), ios::in | ios::binary);
			// check the magic bytes and SNP-major status
			char magic[3];
			_varFile.read(magic, 3);
			if ( (magic[0] != static_cast<char>(0x6C)) || (magic[1] != static_cast<char>(0x1B)) ) {
				cerr << "ERROR: binary file " << _fileName << " not recognized as BED" << endl;
				exit(2);
			} else if (magic[2] != static_cast<char>(0x01)){
				cerr << "ERROR: binary file " << _fileName << " not SNP-major. Run a newer version of plink to fix it." << endl;
				exit(3);
			}
			
		} catch (system_error &error) {
			cerr << "ERROR: cannot open .bed file " << _fileName << " for input: " << error.code().message() << flush;
			perror(" ");
			exit(1);
		}
		try {
			out._varFile.open(out._fileName.c_str(), ios::out);
		} catch (system_error &error) {
			cerr << "ERROR: cannot open .bed file " << out._fileName << " for output: " << error.code().message() << flush;
			perror(" ");
			exit(1);
		}
		while (_varFile) {
			_varFile.read(buf, bufSize);
			out._varFile.write(buf, _varFile.gcount()); // must be gcount() in case we got to the end of the file
		}
		_varFile.close();
		out._varFile.close();
		
		// Copy .bim
		try {
			_bimFile.open(inBim.c_str(), ios::in);
		} catch (system_error &error) {
			cerr << "ERROR: cannot open .bim file " << inBim << " for input: " << error.code().message() << flush;
			perror(" ");
			exit(1);
		}
		try {
			out._bimFile.open(outBim.c_str(), ios::out);
		} catch (system_error &error) {
			cerr << "ERROR: cannot open .bim file " << outBim << " for output: " << error.code().message() << flush;
			perror(" ");
			exit(1);
		}
		while (_bimFile) {
			_bimFile.read(buf, bufSize);
			out._bimFile.write(buf, _bimFile.gcount());
		}
		_bimFile.close();
		out._bimFile.close();
		delete [] buf;
		return;
	}
	
	// Passed all the tests, proceed to sampling
	char *SNPbuf;
	try {
		SNPbuf = new char[Nbed]; // .read() does not append a null-character
	} catch (bad_alloc &) {
		cerr << "ERROR: failed to allocate the SNP buffer" << endl;
		exit(4);

	}
	string bimLine;
	
	try {
		_bimFile.open(inBim.c_str(), ios::in);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open .bim file " << inBim << " for input: " << error.code().message() << flush;
		perror(" ");
		exit(1);
	}
	try {
		out._bimFile.open(outBim.c_str(), ios::out);
		
	} catch (system_error &error) {
		cerr << "ERROR: cannot open .bim file " << outBim << " for output: " << error.code().message() << flush;
		perror(" ");
		exit(1);
	}
	
	try {
		_varFile.open(_fileName.c_str(), ios::in | ios::binary);
		// check the magic bytes and SNP-major status
		char magic[3];
		_varFile.read(magic, 3);
		if ( (magic[0] != static_cast<char>(0x6C)) || (magic[1] != static_cast<char>(0x1B)) ) {
			cerr << "ERROR: binary file " << _fileName << " not recognized as BED" << endl;
			exit(2);
		} else if (magic[2] != static_cast<char>(0x01)){
			cerr << "ERROR: binary file " << _fileName << " not SNP-major. Run a newer version of plink to fix it." << endl;
			exit(3);
		}
	} catch (system_error &error) {
		cerr << "ERROR: cannot open .bed file " << _fileName << " for input: " << error.code().message() << flush;
		perror(" ");
		exit(1);
	}
	try {
		out._varFile.open(out._fileName.c_str(), ios::out | ios::binary | ios::trunc);
		// write the magic bytes and the SNP-major status byte
		char magic[3] = {static_cast<char>(0x6C), static_cast<char>(0x1B), static_cast<char>(0x01)};
		out._varFile.write(magic, 3);

	} catch (system_error &error) {
		cerr << "ERROR: cannot open .bed file " << out._fileName << " for output: " << error.code().message() << flush;
		perror(" ");
		exit(1);
	}
	
	uint64_t nloc = n; // local copy of n
	uint64_t S;
	uint64_t cumS = 0; // cumulative position in the file
	RanDraw *snpSamp;  // so that I can catch RanDraw constructor exceptions
	try {
		snpSamp = new RanDraw();
	} catch (string error) {
		cerr << "ERROR: " << error << endl;
		exit(5);
	}
	while (nloc) {
		try {
			S = snpSamp->vitter(nloc, N); // sample the number of SNPs to skip; keep track of the running total
		} catch (string error) {
			cerr << "ERROR: " << error << endl;
			exit(5);
		}
		cumS += S;
		N    -= S + 1;
		_varFile.seekg(cumS*Nbed + 3); // seekg() index is base-0
		_varFile.read(SNPbuf, Nbed);
		cumS++;                        // step up cumS because seekg() will not change it (unlike the getline() that automatially advances)
		out._varFile.write(SNPbuf, Nbed);
		while (S) { // skipping S lines in the bim file
			_bimFile.ignore(numeric_limits<streamsize>::max(), '\n');
			S--;
		}
		getline(_bimFile, bimLine);
		out._bimFile << bimLine << endl;
		nloc--;
		
	}
	_varFile.close();
	out._varFile.close();
	_bimFile.close();
	out._bimFile.close();
	delete [] SNPbuf;
	delete snpSamp;
	
}

void BedFileI::sampleLD(const uint64_t &n){
	if (n == 0) {
		cerr << "WARNING: zero SNPs requested. Nothing to be done." << endl;
		return;
	}
	
    if (_bimFile.is_open()) {
        _bimFile.close();
    }
    try {
        string bimName = _fileStub + ".bim";
        _bimFile.open(bimName.c_str(), ios::in);
    } catch (system_error &error) {
        cerr << "ERROR: cannot open BIM file for input: " << error.code().message() << flush;
        perror(" ");
        exit(1);
    }
	uint64_t nSamp = 2UL * n; // n is the number of pairs
	
    const uint16_t Npad  = _nCols % 4UL;
    const uint64_t Nbed  = (_nCols/4UL) + static_cast<uint64_t>(Npad > 0UL);
    uint64_t N           = _numLines();
	const string outName = _fileStub + "_LD.tsv";
	if (nSamp > N) {
		cerr << "ERROR: requested a sample of " << nSamp << " SNPs that is greater than the number of SNPs (" << N << ") in the input file." << endl;
		return;
	} else if (nSamp == N) {
		cerr << "WARNING: the number of sampled SNPs the same as the total number of SNPs in the file. LD will be between neighboring loci." << endl;
		// use neighboring loci
		char *locus1  = new char[Nbed];
		char *locus2  = new char[Nbed];
		uint64_t cumS = 0; // cumulative position in the file
		double rSq;
		double Dprime;
		double cnt1;
		double cnt2;
		string chrID;   // chromosome ID
		uint64_t pos1;  // position of locus1
		uint64_t pos2;  // position of locus2
		string bimLine; // line in the .bim file
		string field;   // single field in the .bim file
		stringstream lineStream;
		
		ofstream outFile(outName);
		outFile << "ChrID\tpos1\tpos2\tDistance\trSq\tDprime\tn1\tn2"  << endl;
		
        try {
            _varFile.open(_fileName.c_str(), ios::in | ios::binary);
            // check the magic bytes and SNP-major status
            char magic[3];
            _varFile.read(magic, 3);
            if ( (magic[0] != static_cast<char>(0x6C)) || (magic[1] != static_cast<char>(0x1B)) ) {
                cerr << "ERROR: binary file " << _fileName << " not recognized as BED" << endl;
                exit(2);
            } else if (magic[2] != static_cast<char>(0x01)){
                cerr << "ERROR: binary file " << _fileName << " not SNP-major. Run a newer version of plink to fix it." << endl;
                exit(3);
            }
        } catch (system_error &error) {
            cerr << "ERROR: cannot open .bed file " << _fileName << " for input: " << error.code().message() << flush;
            perror(" ");
            exit(1);
        }
		while (nSamp) {
			// first locus in the pair
			_varFile.seekg(cumS*Nbed + 3); // seekg() index is base-0
			_varFile.read(locus1, Nbed);
			getline(_bimFile, bimLine);
			cumS++;
			nSamp--;
			
			lineStream.str(bimLine); // replaces whatever was there from before
			lineStream >> field;
			chrID = field;
			lineStream >> field >> field >> field;
			pos1 = atoi(field.c_str());
			
			// second locus
			_varFile.seekg(cumS*Nbed + 3);
			_varFile.read(locus2, Nbed);
			getline(_bimFile, bimLine);
			cumS++;
			nSamp--;
			
			lineStream.str(bimLine);
			lineStream >> field;
			if (chrID != field) { // if the chromosome number changes, we discard this pair and continue
				continue;
			}
			lineStream >> field >> field >> field;
			pos2 = atoi(field.c_str());
			if (pos2 <= pos1) {
				cerr << "WARNING: second SNP position (" << pos2 << ") no larger than the first (" << pos1 << ") on chromosome " << chrID << ". Skipping this pair without re-trying." << endl;
				continue;
			}
			_ld(locus1, locus2, Nbed, Npad, rSq, Dprime, cnt1, cnt2);
			outFile << chrID << "\t" << pos1 << "\t" << pos2 << "\t" << pos2 - pos1 << "\t" << rSq << "\t" << Dprime << "\t" << cnt1 << "\t" << cnt2 << endl;
		}
		
		delete [] locus1;
		delete [] locus2;
		_varFile.close();
		_bimFile.close();
		outFile.close();

	}
	try {
		_varFile.open(_fileName.c_str(), ios::in | ios::binary);
		// check the magic bytes and SNP-major status
		char magic[3];
		_varFile.read(magic, 3);
		if ( (magic[0] != static_cast<char>(0x6C)) || (magic[1] != static_cast<char>(0x1B)) ) {
			cerr << "ERROR: binary file " << _fileName << " not recognized as BED" << endl;
			exit(2);
		} else if (magic[2] != static_cast<char>(0x01)){
			cerr << "ERROR: binary file " << _fileName << " not SNP-major. Run a newer version of plink to fix it." << endl;
			exit(3);
		}
	} catch (system_error &error) {
		cerr << "ERROR: cannot open BED file " << _fileName << ": " << error.code().message() << flush;
		perror(" ");
		exit(1);
	}
	// sample loci and calculate LD
	char *locus1  = new char[Nbed];
	char *locus2  = new char[Nbed];
	uint64_t cumS = 0; // cumulative position in the file
	uint64_t S;
	double rSq;
	double Dprime;
	double cnt1;
	double cnt2;
	string chrID;   // chromosome ID
	uint64_t pos1;  // position of locus1
	uint64_t pos2;  // position of locus2
	string bimLine; // line in the .bim file
	string field;   // single field in the .bim file
	stringstream lineStream;
	
	ofstream outFile(outName);
	outFile << "ChrID\tpos1\tpos2\tDistance\trSq\tDprime\tn1\tn2"  << endl;
	
	RanDraw *snpSamp;  // so that I can catch RanDraw constructor exceptions
	try {
		snpSamp = new RanDraw();
	} catch (string error) {
		cerr << "ERROR: " << error << endl;
		exit(5);
	}

	while (nSamp) {
		// first locus in the pair
		try {
			S = snpSamp->vitter(nSamp, N); // sample the number of SNPs to skip; keep track of the running total
		} catch (string error) {
			cerr << "ERROR: " << error << endl;
			exit(5);
		}
		cumS += S;
		N    -= S + 1;
		_varFile.seekg(cumS*Nbed + 3); // seekg() index is base-0
		_varFile.read(locus1, Nbed);
		cumS++;                        // step up cumS because seekg() will not change it (unlike the getline() that automatially advances)
		while (S) { // skipping S lines in the bim file
			_bimFile.ignore(numeric_limits<streamsize>::max(), '\n');
			S--;
		}
		getline(_bimFile, bimLine);
		lineStream.str(bimLine); // replaces whatever was there from before
		lineStream >> field;
		chrID = field;
		lineStream >> field >> field >> field;
		pos1 = atoi(field.c_str());
		nSamp--;
		
		// second locus
		try {
			S = snpSamp->vitter(nSamp, N);
		} catch (string error) {
			cerr << "ERROR: " << error << endl;
			exit(5);
		}
		cumS += S;
		N    -= S + 1;
		_varFile.seekg(cumS*Nbed + 3);
		_varFile.read(locus2, Nbed);
		cumS++;
		while (S) {
			_bimFile.ignore(numeric_limits<streamsize>::max(), '\n');
			S--;
		}
		getline(_bimFile, bimLine);
		
		lineStream.str(bimLine);
		lineStream >> field;
		if (chrID != field) { // if the chromosome number changes, we discard this pair and start over
			nSamp++;          // back up so that we end up with the correct number of pairs regardless of how many chromosome boundaries we cross
			continue;
		}
		lineStream >> field >> field >> field;
		pos2 = atoi(field.c_str());
		if (pos2 <= pos1) {
			cerr << "WARNING: second SNP position (" << pos2 << ") no larger than the first (" << pos1 << ") on chromosome " << chrID << ". Skipping this pair without re-trying." << endl;
			nSamp--; // will not re-try in this case
			continue;
		}
		_ld(locus1, locus2, Nbed, Npad, rSq, Dprime, cnt1, cnt2);
		outFile << chrID << "\t" << pos1 << "\t" << pos2 << "\t" << pos2 - pos1 << "\t" << rSq << "\t" << Dprime << "\t" << cnt1 << "\t" << cnt2 << endl;
		nSamp--;
	}
	
	delete [] locus1;
	delete [] locus2;
	_varFile.close();
	_bimFile.close();
	outFile.close();
}

void BedFileI::sampleLD(const PopIndex &popID, const uint64_t &n){
	if (n == 0) {
		cerr << "WARNING: zero SNPs requested. Nothing to be done." << endl;
		return;
	}
    if (_varFile.is_open()) {
        _varFile.close();
    }
	if (!_bimFile.is_open()) {
		try {
			string bimName = _fileStub + ".bim";
			_bimFile.open(bimName.c_str(), ios::in);
		} catch (system_error &error) {
			cerr << "ERROR: cannot open BIM file for input: " << error.code().message() << flush;
			perror(" ");
			exit(1);
		}
	}
	if (_nCols != popID.size()) {
		cerr << "Sample size the population index (" << popID.size() << ") not equal to sample size (" << _nCols << ") in the .fam file" << endl;
		exit(5);
	}
    
	const uint64_t Nbed = (_nCols/4UL) + static_cast<uint64_t>( (_nCols%4UL) > 0UL ); // 4 SNPs/byte plus padding in BED
    uint64_t N          = _numLines();
    uint64_t nSamp      = 2 * n; // n is the number of pairs
    
    try { // has to be here, because _numLines() closes the file
        _varFile.open(_fileName.c_str(), ios::in | ios::binary);
        // check the magic bytes and SNP-major status
        char magic[3];
        _varFile.read(magic, 3);
        if ( (magic[0] != static_cast<char>(0x6C)) || (magic[1] != static_cast<char>(0x1B)) ) {
            cerr << "ERROR: binary file " << _fileName << " not recognized as BED" << endl;
            exit(2);
        } else if (magic[2] != static_cast<char>(0x01)){
            cerr << "ERROR: binary file " << _fileName << " not SNP-major. Run a newer version of plink to fix it." << endl;
            exit(3);
        }
    } catch (system_error &error) {
        cerr << "ERROR: cannot open BED file " << _fileName << ": " << error.code().message() << flush;
        perror(" ");
        exit(1);
    }

	const string outName = _fileStub + "_LD.tsv";
	if (nSamp > N) {
		cerr << "ERROR: requested a sample of " << nSamp << " SNPs that is greater than the number of SNPs (" << N << ") in the input file." << endl;
		return;
	} else if (nSamp == N) {
		cerr << "WARNING: the number of sampled SNPs the same as the total number of SNPs in the file. LD will be between neighboring loci." << endl;
		// use neighboring loci
		char *locus1  = new char[Nbed];
		char *locus2  = new char[Nbed];
		uint64_t cumS = 0; // cumulative position in the file
		vector<double> rSq(popID.popNumber());
		vector<double> Dprime(popID.popNumber());
		vector<double> cnt1(popID.popNumber());
		vector<double> cnt2(popID.popNumber());
		string chrID;   // chromosome ID
		uint64_t pos1;  // position of locus1
		uint64_t pos2;  // position of locus2
		string bimLine; // line in the .bim file
		string field;   // single field in the .bim file
		stringstream lineStream;
		
		ofstream outFile(outName);
		outFile << "ChrID\tpos1\tpos2\tDistance"  << flush;
		for (unsigned short iPop = 1; iPop <= popID.popNumber(); iPop++) {
			outFile << "\trSq_" << iPop << "\tDprime_" << iPop << "\tn1_" << iPop << "\tn2_" << iPop << flush;
		}
		outFile << endl;
		
		while (nSamp) {
			// first locus in the pair
			_varFile.seekg(cumS*Nbed + 3); // seekg() index is base-0
			_varFile.read(locus1, Nbed);
			getline(_bimFile, bimLine);
			cumS++;
			nSamp--;
			
			lineStream.str(bimLine); // replaces whatever was there from before
			lineStream >> field;
			chrID = field;
			lineStream >> field >> field >> field;
			pos1 = atoi(field.c_str());
			
			// second locus
			_varFile.seekg(cumS*Nbed + 3);
			_varFile.read(locus2, Nbed);
			getline(_bimFile, bimLine);
			cumS++;
			nSamp--;
			
			lineStream.str(bimLine);
			lineStream >> field;
			if (chrID != field) { // if the chromosome number changes, we discard this pair and continue
				continue;
			}
			lineStream >> field >> field >> field;
			pos2 = atoi(field.c_str());
			if (pos2 <= pos1) {
				cerr << "WARNING: second SNP position (" << pos2 << ") no larger than the first (" << pos1 << ") on chromosome " << chrID << ". Skipping this pair without re-trying." << endl;
				continue;
			}
			_ld(locus1, locus2, popID, rSq, Dprime, cnt1, cnt2);
			outFile << chrID << "\t" << pos1 << "\t" << pos2 << "\t" << pos2 - pos1 << flush;
			for (size_t iPop = 0; iPop < popID.popNumber(); iPop++) {
				outFile <<  "\t" << rSq[iPop] << "\t" << Dprime[iPop] << "\t" << cnt1[iPop] << "\t" << cnt2[iPop] << flush;
			}
			outFile << endl;
		}
		
		delete [] locus1;
		delete [] locus2;
		_varFile.close();
		_bimFile.close();
		outFile.close();
        return;
	}
	// sample loci and calculate LD
	char *locus1  = new char[Nbed];
	char *locus2  = new char[Nbed];
	uint64_t cumS = 0; // cumulative position in the file
	uint64_t S;
	vector<double> rSq(popID.popNumber());
	vector<double> Dprime(popID.popNumber());
	vector<double> cnt1(popID.popNumber());
	vector<double> cnt2(popID.popNumber());
	string chrID;   // chromosome ID
	uint64_t pos1;  // position of locus1
	uint64_t pos2;  // position of locus2
	string bimLine; // line in the .bim file
	string field;   // single field in the .bim file
	stringstream lineStream;
	
	ofstream outFile(outName);
	outFile << "ChrID\tpos1\tpos2\tDistance"  << flush;
	for (unsigned short iPop = 1; iPop <= popID.popNumber(); iPop++) {
		outFile << "\trSq_" << iPop << "\tDprime_" << iPop << "\tn1_" << iPop << "\tn2_" << iPop << flush;
	}
	outFile << endl;
	
	RanDraw *snpSamp;  // so that I can catch RanDraw constructor exceptions
	try {
		snpSamp = new RanDraw();
	} catch (string error) {
		cerr << "ERROR: " << error << endl;
		exit(5);
	}
	
	while (nSamp) {
		// first locus in the pair
		try {
			S = snpSamp->vitter(nSamp, N); // sample the number of SNPs to skip; keep track of the running total
		} catch (string error) {
			cerr << "ERROR: " << error << endl;
			exit(5);
		}
		cumS += S;
		N    -= S + 1;
		_varFile.seekg(cumS*Nbed + 3); // seekg() index is base-0
		_varFile.read(locus1, Nbed);
		cumS++;                        // step up cumS because seekg() will not change it (unlike the getline() that automatially advances)
		while (S) { // skipping S lines in the bim file
			_bimFile.ignore(numeric_limits<streamsize>::max(), '\n');
 			S--;
		}
		getline(_bimFile, bimLine);
		lineStream.str(bimLine); // replaces whatever was there from before
		lineStream >> field;
		chrID = field;
		lineStream >> field >> field >> field;
		pos1 = atoi(field.c_str());
		nSamp--;
		
		// second locus
		try {
			S = snpSamp->vitter(nSamp, N);
		} catch (string error) {
			cerr << "ERROR: " << error << endl;
			exit(5);
		}
		cumS += S;
		N    -= S + 1;
		_varFile.seekg(cumS*Nbed + 3);
		_varFile.read(locus2, Nbed);
		cumS++;
		while (S) {
			_bimFile.ignore(numeric_limits<streamsize>::max(), '\n');
			S--;
		}
		getline(_bimFile, bimLine);
		
		lineStream.str(bimLine);
		lineStream >> field;
		if (chrID != field) { // if the chromosome number changes, we discard this pair and start over
			nSamp++;          // back up so that we end up with the correct number of pairs regardless of how many chromosome boundaries we cross
			continue;
		}
		lineStream >> field >> field >> field;
		pos2 = atoi(field.c_str());
		if (pos2 <= pos1) {
			cerr << "WARNING: second SNP position (" << pos2 << ") no larger than the first (" << pos1 << ") on chromosome " << chrID << ". Skipping this pair without re-trying." << endl;
			nSamp--; // will not re-try in this case
			continue;
		}
		_ld(locus1, locus2, popID, rSq, Dprime, cnt1, cnt2);
		outFile << chrID << "\t" << pos1 << "\t" << pos2 << "\t" << pos2 - pos1 << flush;
		for (size_t iPop = 0; iPop < popID.popNumber(); iPop++) {
			outFile <<  "\t" << rSq[iPop] << "\t" << Dprime[iPop] << "\t" << cnt1[iPop] << "\t" << cnt2[iPop] << flush;
		}
		outFile << endl;
		nSamp--;
	}
	
	delete [] locus1;
	delete [] locus2;
	_varFile.close();
	_bimFile.close();
	outFile.close();

}

void BedFileO::open(){
	string bimName = _fileStub + ".bim";
	string famName = _fileStub + ".fam";
    if (_varFile.is_open()) {
        _varFile.close();
    }
	try {
		_varFile.open(_fileName.c_str(), ios::out | ios::binary);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open BED file " << _fileName << " for output: " << error.code().message() << flush;
		perror(" ");
		exit(1);
		
	}
	// write the magic bytes and the SNP-major status byte
	char magic[3] = {static_cast<char>(0x6C), static_cast<char>(0x1B), static_cast<char>(0x01)};
	_varFile.write(magic, 3);
	
	try {
		_bimFile.open(bimName.c_str(), ios::out);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open .bim file " << bimName << " for output: " << error.code().message() << flush;
		perror(" ");
		exit(1);
		
	}
	
	try {
		_famFile.open(famName.c_str(), ios::out);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open .fam file " << famName << " for output: " << error.code().message() << flush;
		perror(" ");
		exit(1);
		
	}
	
}

// GtxtFile methods

void GtxtFile::close(){
    if (_varFile.is_open()) {
        _varFile.close();
    }
}

void GtxtFileI::open(){
    if (_varFile.is_open()) {
        _varFile.close();
    }
    try {
        _varFile.open(_fileName.c_str(), ios::in);
    } catch (system_error &error) {
        cerr << "ERROR: cannot open text file " << _fileName << " for input: " << error.code().message() << flush;
        perror(" ");
        exit(1);
        
    }
    
}

uint64_t GtxtFileI::_numLines(){
    uint64_t N = 0;
    
    /*
     * I am using the line-end counting method. It is > 2-fold faster than reading lines with getline().
     */
    if (_varFile.is_open()) {
        _varFile.close();
    }
    try {
        _varFile.open(_fileName.c_str(), ios::in);
    } catch (system_error &error) {
        cerr << "ERROR: cannot open text file " << _fileName << " for input: " << error.code().message() << flush;
        perror(" ");
        exit(1);
        
    }
    
    const size_t bufSize = BUF_SIZE;
    char *buf;
    try {
        buf = new char[bufSize];
    } catch (bad_alloc& error) {
        cerr << "ERROR: failed to allocate buffer in GtxtFile::_numLines(): " << error.what() << endl;
        exit(4);
    }
    
    while (_varFile) {
        _varFile.read(buf, bufSize);
        for (size_t i = 0; i < _varFile.gcount(); i++) {
            if (buf[i] == '\n') {
                N++;
            }
        }
    }
    _varFile.close();
    delete [] buf;
    
    return N - _head; // subtracting 1 for the header if present
    
}

void GtxtFileI::sample(GtxtFileO &out, const uint64_t &n, const bool &headSkip){
    if (n == 0) {
        cerr << "WARNING: zero rows requested. Nothing to be done." << endl;
        return;
    }
    // start by figuring out the number of rows in the file
    if (out._varFile.is_open()) {
        out._varFile.close();
    }
    
    uint64_t N = _numLines(); // Number of rows
    
    try {
        _varFile.open(_fileName.c_str(), ios::in);
    } catch (system_error &error) {
        cerr << "ERROR: cannot open text file " << _fileName << " for input: " << error.code().message() << flush;
        perror(" ");
        exit(1);
    }
    
    try {
        out._varFile.open(out._fileName.c_str(), ios::out | ios::trunc);
    } catch (system_error &error) {
        cerr << "ERROR: cannot open text file " << out._fileName << " for output: " << error.code().message() << flush;
        perror(" ");
        exit(1);
    }
    
    // Test for potential problems
    if (N < n) {
        cerr << "ERROR: requested a sample of " << n << " rows that is greater than the number of rows (" << N << ") in the input file." << endl;
        return;
    } else if (N == n) {
        cerr << "WARNING: sample size (" << n << ") the same as the number of rows (" << N << ") in the input file. Simply copying the file." << endl;
        char *buf;
        const size_t bufSize = BUF_SIZE;
        try {
            buf = new char[bufSize];
        } catch (bad_alloc&) {
            cerr << "ERROR: failed to allocate buffer" << endl;
            exit(4);
        }
        if (_head) {
            if (headSkip) {
                _varFile.ignore(numeric_limits<streamsize>::max(), '\n');
            } else {
                string header;
                getline(_varFile, header);
                out._varFile << header << endl;
            }
        }
        
        // Copy the file
        while (_varFile) {
            _varFile.read(buf, bufSize);
            out._varFile.write(buf, _varFile.gcount());
        }
        _varFile.close();
        out._varFile.close();
        delete [] buf;
        return;
    }
    
    // Passed all the tests, proceed to sampling
    string curLine;
    if (_head) {
        if (headSkip) {
            _varFile.ignore(numeric_limits<streamsize>::max(), '\n');
        } else {
            getline(_varFile, curLine);
            out._varFile << curLine << endl;
            curLine.erase();
        }
    }
    
    uint64_t nloc = n; // local copy of n
    uint64_t S;
    RanDraw *snpSamp;  // so that I can catch RanDraw constructor exceptions
    try {
        snpSamp = new RanDraw();
    } catch (string error) {
        cerr << "ERROR: " << error << endl;
        exit(5);
    }
    while (nloc) {
        try {
            S = snpSamp->vitter(nloc, N); // sample the number of SNPs to skip
        } catch (string error) {
            cerr << "ERROR: " << error << endl;
            exit(5);
        }
        N -= S + 1;
        while (S) { // skipping S lines in the HMP file
            _varFile.ignore(numeric_limits<streamsize>::max(), '\n');
            S--;
        }
        getline(_varFile, curLine);
        out._varFile << curLine << endl;
        nloc--;
        
    }
    _varFile.close();
    out._varFile.close();
    delete snpSamp;
    
}

void GtxtFileI::sample(const uint64_t &n, const bool &headSkip, const char &delim, vector<string> &out){
    if (n == 0) {
        cerr << "WARNING: zero rows requested. Nothing to be done." << endl;
        return;
    }
    if (out.size()) {
        out.resize(0);
    }
    
    uint64_t N = _numLines(); // Number of rows
    
    try {
        _varFile.open(_fileName.c_str(), ios::in);
    } catch (system_error &error) {
        cerr << "ERROR: cannot open text file " << _fileName << " for input: " << error.code().message() << flush;
        perror(" ");
        exit(1);
    }
    
    // Test for potential problems
    if (N < n) {
        cerr << "ERROR: requested a sample of " << n << " rows that is greater than the number of rows (" << N << ") in the input file." << endl;
        return;
    } else if (N == n) {
        cerr << "WARNING: sample size (" << n << ") the same as the number of rows (" << N << ") in the input file. Simply copying the file." << endl;
        string curLine;
        if (_head) {
            if (headSkip) {
                _varFile.ignore(numeric_limits<streamsize>::max(), '\n');
            } else {
                getline(_varFile, curLine);
                string field;
                stringstream lineSS(curLine);
                while (lineSS) {
                    getline(lineSS, field, delim);
                    out.push_back(field);
                    
                }
            }
        }
        
        // Copy the file
        string field;
        while (_varFile) {
            getline(_varFile, curLine);
            stringstream lineSS(curLine);
            while (lineSS) {
                getline(lineSS, field, delim);
                out.push_back(field);
            }
        }
        _varFile.close();
        return;
    }
    
    // Passed all the tests, proceed to sampling
    string curLine;
    string field;
    if (_head) {
        if (headSkip) {
            _varFile.ignore(numeric_limits<streamsize>::max(), '\n');
        } else {
            getline(_varFile, curLine);
            stringstream lineSS(curLine);
            while (lineSS) {
                getline(lineSS, field, delim);
                out.push_back(field);
                
            }
        }
    }
    
    uint64_t nloc = n; // local copy of n
    uint64_t S;
    RanDraw *snpSamp;  // so that I can catch RanDraw constructor exceptions
    try {
        snpSamp = new RanDraw();
    } catch (string error) {
        cerr << "ERROR: " << error << endl;
        exit(5);
    }
    while (nloc) {
        try {
            S = snpSamp->vitter(nloc, N); // sample the number of SNPs to skip
        } catch (string error) {
            cerr << "ERROR: " << error << endl;
            exit(5);
        }
        N -= S + 1;
        while (S) { // skipping S lines in the HMP file
            _varFile.ignore(numeric_limits<streamsize>::max(), '\n');
            S--;
        }
        getline(_varFile, curLine);
        stringstream lineSS(curLine);
        while (lineSS) {
            getline(lineSS, field, delim);
            out.push_back(field);
            
        }
        nloc--;
        
    }
    _varFile.close();
    delete snpSamp;
    
}

void GtxtFileO::open(){
    try {
        _varFile.open(_fileName.c_str(), ios::out);
    } catch (system_error &error) {
        cerr << "ERROR: cannot open text file " << _fileName << " for output: " << error.code().message() << flush;
        perror(" ");
        exit(1);
        
    }
}

//TpedFile methods
TpedFile::~TpedFile(){
    if (_tfamFile.is_open()) {
        _tfamFile.close();
    }
}

void TpedFile::close(){
    if (_varFile.is_open()) {
        _varFile.close();
    }
    if (_tfamFile.is_open()) {
        _tfamFile.close();
    }
}

uint64_t TpedFileI::_famLines(){
	uint64_t N = 0;
	
	/*
	 * I am using the line-end counting method. It is > 2-fold faster than reading lines with getline().
	 */
	if (_tfamFile.is_open()) {
        _tfamFile.close();
    }
    string tfamName = _fileStub + ".tfam";
    try {
        _tfamFile.open(tfamName.c_str(), ios::in);
    } catch (system_error &error) {
        cerr << "ERROR: cannot open .tfam file " << tfamName << " for input: " << error.code().message() << flush;
        perror(" ");
        exit(1);
        
    }
    
	const size_t bufSize = BUF_SIZE; // 10M seems optimal based on experiments on my MacBook Pro (with SSD)
	char *buf;
	try {
		buf = new char[bufSize];
	} catch (bad_alloc& error) {
		cerr << "ERROR: failed to allocate buffer in TpedFile::_numLines(): " << error.what() << endl;
		exit(4);
	}
	
	while (_tfamFile) {
		_tfamFile.read(buf, bufSize);
		for (size_t i = 0; i < _tfamFile.gcount(); i++) {
			if (buf[i] == '\n') {
				N++;
			}
		}
	}
	delete [] buf;
	_tfamFile.close();
	
	return N;
}
uint64_t TpedFileI::_famLines(fstream &fam){
	uint64_t N = 0;
	
	/*
	 * I am using the line-end counting method. It is > 2-fold faster than reading lines with getline().
	 */
	if (_tfamFile.is_open()) {
        _tfamFile.close();
    }
    string tfamName = _fileStub + ".tfam";
    try {
        _tfamFile.open(tfamName.c_str(), ios::in);
    } catch (system_error &error) {
        cerr << "ERROR: cannot open .tfam file " << tfamName << " for input: " << error.code().message() << flush;
        perror(" ");
        exit(1);
    }
    if (!fam.is_open()) {
        throw string("Output .fam filestream not open");
    }
    const size_t bufSize = BUF_SIZE;
    char *buf;
    try {
        buf = new char[bufSize];
    } catch (bad_alloc& error) {
        cerr << "ERROR: failed to allocate buffer in TpedFile::_famLines(): " << error.what() << endl;
        exit(4);
    }
    
    while (_tfamFile) {
        _tfamFile.read(buf, bufSize);
        for (size_t i = 0; i < _tfamFile.gcount(); i++) {
            if (buf[i] == '\n') {
                N++;
            }
        }
        fam.write(buf, _tfamFile.gcount());
    }
    _tfamFile.close();
    fam.close();
    delete [] buf;

	return N;
}

void TpedFileI::_famCopy(fstream &fam){
	if (_tfamFile.is_open()) {
        _tfamFile.close();
    }
    string tfamName = _fileStub + ".tfam";
    try {
        _tfamFile.open(tfamName.c_str(), ios::in);
    } catch (system_error &error) {
        cerr << "ERROR: cannot open .tfam file " << tfamName << " for input: " << error.code().message() << flush;
        perror(" ");
        exit(1);
        
    }
    if (!fam.is_open()) {
        throw string("Output .fam filestream not open");
    }

    const size_t bufSize = BUF_SIZE;
    char *buf;
    try {
        buf = new char[bufSize];
    } catch (bad_alloc& error) {
        cerr << "ERROR: failed to allocate buffer in TpedFile::_famCopy(): " << error.what() << endl;
        exit(4);
    }
    
    while (_tfamFile) {
        _tfamFile.read(buf, bufSize);
        fam.write(buf, _tfamFile.gcount());
    }
    _tfamFile.close();
    delete [] buf;

}

uint64_t TpedFileI::_numLines(){
    uint64_t N = 0;
    
    /*
     * I am using the line-end counting method. It is > 2-fold faster than reading lines with getline().
     */
    if (_varFile.is_open()) {
        _varFile.close();
    }
    try {
        _varFile.open(_fileName.c_str(), ios::in);
    } catch (system_error &error) {
        cerr << "ERROR: cannot open .tped file " << _fileName << " for input: " << error.code().message() << flush;
        perror(" ");
        exit(1);
        
    }
    
    const size_t bufSize = BUF_SIZE;
    char *buf;
    try {
        buf = new char[bufSize];
    } catch (bad_alloc& error) {
        cerr << "ERROR: failed to allocate buffer in TpedFile::_numLines(): " << error.what() << endl;
        exit(4);
    }
    
    while (_varFile) {
        _varFile.read(buf, bufSize);
        for (size_t i = 0; i < _varFile.gcount(); i++) {
            if (buf[i] == '\n') {
                N++;
            }
        }
    }
    _varFile.close();
    delete [] buf;
    
    return N;
}

void TpedFileI::open(){
	string tfamName = _fileStub + ".tfam";
	
	try {
		_varFile.open(_fileName.c_str(), ios::in);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open TPED file " << _fileName << " for input: " << error.code().message() << flush;
		perror(" ");
		exit(1);
		
	}
	
	try {
		_tfamFile.open(tfamName.c_str(), ios::in);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open .tfam file " << tfamName << " for input: " << error.code().message() << flush;
		perror(" ");
		exit(1);
		
	}
	
}

void TpedFileI::sample(TpedFileO &out, const uint64_t &n){
	
	if (n == 0) {
		cerr << "WARNING: zero SNPs requested. Nothing to be done." << endl;
		return;
	}
	// start by figuring out the number of SNPs in the .bed file
	if (out._varFile.is_open()) {
		out.close();
	}
	
	string inFam  = _fileStub + ".tfam";
	string outFam = out._fileStub + ".tfam";
	
	try {
		_tfamFile.open(inFam.c_str(), ios::in);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open .tfam file " << inFam << " for input: " << error.code().message() << flush;
		perror(" ");
		exit(1);
	}
	
	try {
		out._tfamFile.open(outFam.c_str(), ios::out | ios::trunc);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open .tfam file " << outFam << " for output: " << error.code().message() << flush;
		perror(" ");
		exit(1);
	}
	
	// Copy the .tfam file
	_famCopy(out._tfamFile);
	// Calculate number of SNPs in the input TPED file
	try {
		_varFile.open(_fileName.c_str(), ios::in);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open BED file " << _fileName << ": " << error.code().message() << flush;
		perror(" ");
		exit(1);
	}
	
	uint64_t N = _numLines(); // Number of SNPs
	
	// Test for potential problems
	if (N < n) {
		cerr << "ERROR: requested a sample of " << n << " SNPs that is greater than the number of SNPs (" << N << ") in the input file." << endl;
		return;
	} else if (N == n) {
		cerr << "WARNING: sample size (" << n << ") the same as the number of SNPs (" << N << ") in the input file. Simply copying the files." << endl;
		char *buf;
		const size_t bufSize = BUF_SIZE;
		try {
			buf = new char[bufSize];
		} catch (bad_alloc&) {
			cerr << "ERROR: failed to allocate buffer" << endl;
			exit(4);
		}
		
		// Copy .tped
		try {
			_varFile.open(_fileName.c_str(), ios::in);
			
		} catch (system_error &error) {
			cerr << "ERROR: cannot open .tped file " << _fileName << " for input: " << error.code().message() << flush;
			perror(" ");
			exit(1);
		}
		try {
			out._varFile.open(out._fileName.c_str(), ios::out);
		} catch (system_error &error) {
			cerr << "ERROR: cannot open .tped file " << out._fileName << " for output: " << error.code().message() << flush;
			perror(" ");
			exit(1);
		}
		while (_varFile) {
			_varFile.read(buf, bufSize);
			out._varFile.write(buf, _varFile.gcount());
		}
		_varFile.close();
		out._varFile.close();
		delete [] buf;
		return;
	}
	
	// Passed all the tests, proceed to sampling
	string tpedLine;
	
	try {
		_varFile.open(_fileName.c_str(), ios::in);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open .tped file " << _fileName << " for input: " << error.code().message() << flush;
		perror(" ");
		exit(1);
	}
	try {
		out._varFile.open(out._fileName.c_str(), ios::out| ios::trunc);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open .tped file " << out._fileName << " for output: " << error.code().message() << flush;
		perror(" ");
		exit(1);
	}
	
	uint64_t nloc = n; // local copy of n
	uint64_t S;
	RanDraw *snpSamp;  // so that I can catch RanDraw constructor exceptions
	try {
		snpSamp = new RanDraw();
	} catch (string error) {
		cerr << "ERROR: " << error << endl;
		exit(5);
	}
	while (nloc) {
		try {
			S = snpSamp->vitter(nloc, N); // sample the number of SNPs to skip
		} catch (string error) {
			cerr << "ERROR: " << error << endl;
			exit(5);
		}
		N -= S + 1;
		while (S) { // skipping S lines in the tped file
			_varFile.ignore(numeric_limits<streamsize>::max(), '\n');
			S--;
		}
		getline(_varFile, tpedLine);
		out._varFile << tpedLine << endl;
		nloc--;
		
	}
	_varFile.close();
	out._varFile.close();
	delete snpSamp;

}

void TpedFileO::open(){
	string tfamName = _fileStub + ".tfam";
	
	try {
		_varFile.open(_fileName.c_str(), ios::out);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open TPED file " << _fileName << " for output: " << error.code().message() << flush;
		perror(" ");
		exit(1);
		
	}
	
	try {
		_tfamFile.open(tfamName.c_str(), ios::out);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open .tfam file " << tfamName << " for output: " << error.code().message() << flush;
		perror(" ");
		exit(1);
		
	}
	
}

// VcfFile methods
void VcfFile::close(){
    if (_varFile.is_open()) {
        _varFile.close();
    }
}

uint64_t VcfFileI::_numLines(){
	uint64_t N = 0;
	
	/*
	 * I am using the line-end counting method. It is > 2-fold faster than reading lines with getline().
	 */
	if (_varFile.is_open()) {
        _varFile.close();
    }
    try {
        _varFile.open(_fileName.c_str(), ios::in);
    } catch (system_error &error) {
        cerr << "ERROR: cannot open .vcf file " << _fileName << " for input: " << error.code().message() << flush;
        perror(" ");
        exit(1);
        
    }
	string trash;
	// skip header lines. Will also read the one after the last header line.
	while (1) {
		getline(_varFile, trash);
		if (trash[0] != '#') {
			break;
		}
	}
	
	const size_t bufSize = BUF_SIZE;
	char *buf;
	try {
		buf = new char[bufSize];
	} catch (bad_alloc& error) {
		cerr << "ERROR: failed to allocate buffer in VcfFileI::_numLines(): " << error.what() << endl;
		exit(4);
	}
	
	while (_varFile) {
		_varFile.read(buf, bufSize);
		for (size_t i = 0; i < _varFile.gcount(); i++) {
			if (buf[i] == '\n') {
				N++;
			}
		}
	}
	_varFile.close();
	delete [] buf;
	
	return N + 1; // adding 1 for the line immediately after the header
}

void VcfFileI::open(){
    if (_varFile.is_open()) {
        _varFile.close();
    }
    try {
		_varFile.open(_fileName.c_str(), ios::in);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open VCF file " << _fileName << " for input: " << error.code().message() << flush;
		perror(" ");
		exit(1);
		
	}

}

void VcfFileI::sample(VcfFileO &out, const uint64_t &n){
	if (n == 0) {
		cerr << "WARNING: zero SNPs requested. Nothing to be done." << endl;
		return;
	}
	// start by figuring out the number of SNPs in the file
	if (out._varFile.is_open()) {
		out.close();
	}
	
	
	uint64_t N = _numLines(); // Number of SNPs
	try {
		_varFile.open(_fileName.c_str(), ios::in);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open .vcf file " << _fileName << " for input: " << error.code().message() << flush;
		perror(" ");
		exit(1);
	}
	
	try {
		out._varFile.open(out._fileName.c_str(), ios::out | ios::trunc);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open .vcf file " << out._fileName << " for output: " << error.code().message() << flush;
		perror(" ");
		exit(1);
	}

	// Skip and copy the header
	string headLine;
	while (1) {
		getline(_varFile, headLine);
		if (headLine[0] != '#') {
			fpos_t cur = _varFile.tellg();
			if (cur < headLine.size()) {
				cerr << "ERROR: current position (" << cur << ") smaller than header line length (" << headLine.size() << ")" << endl;
				exit(6);
			}
			_varFile.seekg(cur - headLine.size() - 1); // go back to the start of the line, since this is the first non-header line
			break;
		}
		out._varFile << headLine << endl;
	}
	
	// Test for potential problems
	if (N < n) {
		cerr << "ERROR: requested a sample of " << n << " SNPs that is greater than the number of SNPs (" << N << ") in the input file." << endl;
		return;
	} else if (N == n) {
		cerr << "WARNING: sample size (" << n << ") the same as the number of SNPs (" << N << ") in the input file. Simply copying the files." << endl;
		char *buf;
		const size_t bufSize = BUF_SIZE;
		try {
			buf = new char[bufSize];
		} catch (bad_alloc&) {
			cerr << "ERROR: failed to allocate buffer" << endl;
			exit(4);
		}
		
		// Copy .vcf
		while (_varFile) {
			_varFile.read(buf, bufSize);
			out._varFile.write(buf, _varFile.gcount());
		}
		_varFile.close();
		out._varFile.close();
		delete [] buf;
		return;
	}
	
	// Passed all the tests, proceed to sampling
	string curLine;
	
	uint64_t nloc = n; // local copy of n
	uint64_t S;
	RanDraw *snpSamp;  // so that I can catch RanDraw constructor exceptions
	try {
		snpSamp = new RanDraw();
	} catch (string error) {
		cerr << "ERROR: " << error << endl;
		exit(5);
	}
	while (nloc) {
		try {
			S = snpSamp->vitter(nloc, N); // sample the number of SNPs to skip
		} catch (string error) {
			cerr << "ERROR: " << error << endl;
			exit(5);
		}
		N -= S + 1;
		while (S) { // skipping S lines in the vcf file
			_varFile.ignore(numeric_limits<streamsize>::max(), '\n');
 			S--;
		}
		getline(_varFile, curLine);
		out._varFile << curLine << endl;
		nloc--;
		
	}
	_varFile.close();
	out._varFile.close();
	delete snpSamp;

}

void VcfFileO::open(){
    if (_varFile.is_open()) {
        _varFile.close();
    }
    try {
		_varFile.open(_fileName.c_str(), ios::out);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open VCF file " << _fileName << " for output: " << error.code().message() << flush;
		perror(" ");
		exit(1);
		
	}

}

// HmpFile methods

void HmpFile::close(){
	if (_varFile.is_open()) {
		_varFile.close();
	}
}

HmpFileI::HmpFileI(const string &fileName) : HmpFile(fileName) {
    if (_varFile.is_open()) {
        _varFile.close();
    }
    try {
        _varFile.open(_fileName.c_str(), ios::in);
    } catch (system_error &error) {
        cerr << "ERROR: cannot open HMP file " << _fileName << " for input: " << error.code().message() << flush;
        perror(" ");
        exit(1);
        
    }
    char buf[3];
    _varFile.read(buf, 3);
    if ( (buf[0] == 'r') && (buf[1] == 's') && (buf[2] == '#') ) {
        _head = true;
    }
    _varFile.close();
    
}

void HmpFileI::open(){
    if (_varFile.is_open()) {
        _varFile.close();
    }
    try {
		_varFile.open(_fileName.c_str(), ios::in);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open HMP file " << _fileName << " for input: " << error.code().message() << flush;
		perror(" ");
		exit(1);
	}

}

uint64_t HmpFileI::_numLines() {
    uint64_t N = 0;
    
    /*
     * I am using the line-end counting method. It is > 2-fold faster than reading lines with getline().
     */
    if (_varFile.is_open()) {
        _varFile.close();
    }
    try {
        _varFile.open(_fileName.c_str(), ios::in);
    } catch (system_error &error) {
        cerr << "ERROR: cannot open HMP file " << _fileName << " for input: " << error.code().message() << flush;
        perror(" ");
        exit(1);
        
    }
    // skip header line if present.
    string trash;
    if (_head) {
        getline(_varFile, trash);
    }
    
    const size_t bufSize = BUF_SIZE;
    char *buf;
    try {
        buf = new char[bufSize];
    } catch (bad_alloc& error) {
        cerr << "ERROR: failed to allocate buffer in HmpFileI::_numLines(): " << error.what() << endl;
        exit(4);
    }
    
    while (_varFile) {
        _varFile.read(buf, bufSize);
        for (size_t i = 0; i < _varFile.gcount(); i++) {
            if (buf[i] == '\n') {
                N++;
            }
        }
    }
    _varFile.close();
    delete [] buf;
    
    return N;
}
void HmpFileI::sample(HmpFileO &out, const uint64_t &n){
	if (n == 0) {
		cerr << "WARNING: zero SNPs requested. Nothing to be done." << endl;
		return;
	}
	// start by figuring out the number of SNPs in the file
	if (out._varFile.is_open()) {
		out.close();
	}
	
	uint64_t N = _numLines(); // Number of SNPs
	
	try {
		_varFile.open(_fileName.c_str(), ios::in);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open HMP file " << _fileName << " for input: " << error.code().message() << flush;
		perror(" ");
		exit(1);
	}
	
	try {
		out._varFile.open(out._fileName.c_str(), ios::out | ios::trunc);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open HMP file " << out._fileName << " for output: " << error.code().message() << flush;
		perror(" ");
		exit(1);
	}
	
	// Test for potential problems
	if (N < n) {
		cerr << "ERROR: requested a sample of " << n << " SNPs that is greater than the number of SNPs (" << N << ") in the input file." << endl;
		return;
	} else if (N == n) {
		cerr << "WARNING: sample size (" << n << ") the same as the number of SNPs (" << N << ") in the input file. Simply copying the files." << endl;
		char *buf;
		const size_t bufSize = BUF_SIZE;
		try {
			buf = new char[bufSize];
		} catch (bad_alloc&) {
			cerr << "ERROR: failed to allocate buffer" << endl;
			exit(4);
		}
		
		// Copy .hmp
		while (_varFile) {
			_varFile.read(buf, bufSize);
			out._varFile.write(buf, _varFile.gcount());
		}
		_varFile.close();
		out._varFile.close();
		delete [] buf;
		return;
	}
	
	// Passed all the tests, proceed to sampling
	string curLine;
	char tag[3];
	_varFile.read(tag, 3);
	_varFile.seekg(0);
	if ( (tag[0] == 'r') && (tag[1] == 's') && (tag[2] == '#') ) {
		string header;
		getline(_varFile, header);
		out._varFile << header << endl;
	}
	
	uint64_t nloc = n; // local copy of n
	uint64_t S;
	RanDraw *snpSamp;  // so that I can catch RanDraw constructor exceptions
	try {
		snpSamp = new RanDraw();
	} catch (string error) {
		cerr << "ERROR: " << error << endl;
		exit(5);
	}
	while (nloc) {
		try {
			S = snpSamp->vitter(nloc, N); // sample the number of SNPs to skip
		} catch (string error) {
			cerr << "ERROR: " << error << endl;
			exit(5);
		}
		N -= S + 1;
		while (S) { // skipping S lines in the HMP file
			_varFile.ignore(numeric_limits<streamsize>::max(), '\n');
			S--;
		}
		getline(_varFile, curLine);
		out._varFile << curLine << endl;
		nloc--;
		
	}
	_varFile.close();
	out._varFile.close();
	delete snpSamp;

}

void HmpFileO::open(){
    if (_varFile.is_open()) {
        _varFile.close();
    }
    try {
		_varFile.open(_fileName.c_str(), ios::out);
	} catch (system_error &error) {
		cerr << "ERROR: cannot open HMP file " << _fileName << " for output: " << error.code().message() << flush;
		perror(" ");
		exit(1);
		
	}
}






