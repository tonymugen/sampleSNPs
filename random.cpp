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

/// Random number generation
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2017 Anthony J. Greenberg
 * \version 0.1
 *
 * Class implementation for facilities that generate random draws from various distributions.
 *
 */

#include <string>
#include <cstring>
#include <cstdint>
#include <cmath>

#include "random.hpp"

using std::string;

using namespace sampFiles;

// GenerateHR methods
volatile uint64_t GenerateHR::ranInt(){
	uint64_t rInt = 0;
	
	unsigned char ok  = '\0';
	string error;
	
	unsigned short keepTrying = 10; // number of tries before giving up
	while (keepTrying) {
		keepTrying--;
		
		asm volatile ("rdrand %0; setc %1"
					  : "=r" (rInt), "=qm" (ok)
					  :
					  : "cc"
					  );
		if (ok) {
			break;
		}
	}
	
	if (!ok) {
		error = "RDRAND_failed";
		throw error;
	}
	
	return rInt;
}

// GenerateMT static members
const unsigned short GenerateMT::_n = 312;
const unsigned short GenerateMT::_m = 156;
const uint64_t GenerateMT::_lm      = static_cast<uint64_t>(0xFFFFFFFF80000000);
const uint64_t GenerateMT::_um      = static_cast<uint64_t>(0x7FFFFFFF);
const uint64_t GenerateMT::_b       = static_cast<uint64_t>(0x71D67FFFEDA60000);
const uint64_t GenerateMT::_c       = static_cast<uint64_t>(0xFFF7EEE000000000);
const uint64_t GenerateMT::_d       = static_cast<uint64_t>(0x5555555555555555);
const unsigned int GenerateMT::_l   = 43;
const unsigned int GenerateMT::_s   = 17;
const unsigned int GenerateMT::_t   = 37;
const unsigned int GenerateMT::_u   = 29;
const uint64_t GenerateMT::_alt[2]  = {static_cast<uint64_t>(0), static_cast<uint64_t>(0xB5026F5AA96619E9)};

//GenerateMT methods
GenerateMT::GenerateMT(){
	// start by using RDTSC to set the seed
	unsigned int lo = 0;
	unsigned int hi = 0;
	uint64_t seed   = 0;
	asm volatile ("rdtsc;"
				   : "=a" (lo), "=d" (hi)
	);
	seed = (static_cast<uint64_t>(hi)<<32)|lo;
	uint64_t f = 6364136223846793005ULL;
	_mt[0]     = seed;
	_mti       = 1;
	
	for (; _mti < _n; _mti++) {
		_mt[_mti] = (f * (_mt[_mti - 1]^(_mt[_mti - 1] >> 62)) + _mti);
	}
	
}

volatile uint64_t GenerateMT::ranInt(){
	
	if (_mti == _n) { // do we need to re-run the twister?
		size_t i = 0;
		for (; i < _n - _m; i++) { // first _m words
			_x     = (_mt[i]&_um)|(_mt[i + 1]&_lm);
			_mt[i] = _mt[i + _m]^(_x>>1)^_alt[static_cast<size_t>(_x&1ULL)];
		}
		for (; i < _n - 1; i++) { // rest, except for the last element
			_x     = (_mt[i]&_um)|(_mt[i+1]&_lm);
			_mt[i] = _mt[i + (_m - _n)] ^ (_x>>1)^_alt[static_cast<size_t>(_x&1ULL)];
		}
		// now set the last element
		_x          = (_mt[_n - 1]&_um)|(_mt[0]&_lm);
		_mt[_n - 1] = _mt[_m - 1]^(_x>>1)^_alt[static_cast<size_t>(_x&1ULL)];
		
		_mti = 0;

	}
	// extract pseudo-random number
	_x = _mt[_mti++];
	
	_x ^= ((_x>>_u)&_d);
	_x ^= ((_x<<_s)&_b);
	_x ^= ((_x<<_t)&_c);
	_x ^= (_x>>_l);
	
	return _x;
}

RanDraw::RanDraw(){
	string error;
	
	unsigned int eax;
	unsigned int ebx;
	unsigned int ecx;
	unsigned int edx;
	
	unsigned int leaf    = 0;
	unsigned int subleaf = 0;
	
	// first look at vendor info
	asm volatile ("cpuid;"
				  : "=a" (eax), "=b" (ebx), "=c" (ecx), "=d" (edx)
				  : "a" (leaf), "c" (subleaf)
				  );
	
	if ( !(
		(memcmp(reinterpret_cast<char*>(&ebx), "Genu", 4) && memcmp(reinterpret_cast<char*>(&edx), "ineI", 4) && memcmp(reinterpret_cast<char*>(&ecx), "ntel", 4)) ||
		(memcmp(reinterpret_cast<char*>(&ebx), "Auth", 4) && memcmp(reinterpret_cast<char*>(&edx), "enti", 4) && memcmp(reinterpret_cast<char*>(&ecx), "cAMD", 4))
		)) {
		error = "CPU_unsupported";
		throw error;
	}
	
	
	// now test for RDRAND. bit 30 of ECX when CPUINFO invoked with EAX set to 01H
	leaf = 0x1;
	asm volatile ("cpuid;"
				  : "=a" (eax), "=b" (ebx), "=c" (ecx), "=d" (edx)
				  : "a" (leaf), "c" (subleaf)
				  );
	
	if ( (ecx & 0x40000000) == 0x40000000) { // bit mask 0x40000000 for the 30th bit
		_rng = new GenerateHR(); // have hardware random numbers
	} else {
		_rng = new GenerateMT(); // use MT
	}

}

volatile double RanDraw::runifnz(){
	double rnz = 0.0;
	do {
		rnz = this->runif();
	} while (rnz == 0.0); // simply reject 0.0
	
	return rnz;
}

volatile uint64_t RanDraw::vitterA(const double &n, const double &N){
	// The notation follows Vitter's (1987) as closely as possible
	// Note that my runif() is on [0,1] (Vitter assumes (0,1)), so I have to sometimes adjust accordingly
	uint64_t S  = 0;
	double top  = N - n;
	double quot = top/N;
	double v;
	
	// some trivial conditions first
	if ( (n == 0) || (n > N) ) {
		return S;
	} else if (n == 1) {
		do {
			S = static_cast<uint64_t>(floor(N * this->runif()));
			
		} while (S == static_cast<uint64_t>(N)); // s == N would be a problem because it would overshoot the end of the array/file (in base-0 space); effectively sampling [0,1) uniform to prevent this
		return S;
	}
	// [0,1) uniform
	do {
		v = this->runif();
	} while (v == 1.0); // i.e. only repeat if v == 1.0
	
	double Nloc = N;
	while (quot > v) {
		S++;
		top--;
		Nloc--;
		quot = quot*top/Nloc;
	}
	
	return S;
}

volatile uint64_t RanDraw::vitter(const double &n, const double &N){
	// The notation follows Vitter's (1987) as closely as possible
	// Note that my runif() is on [0,1] (Vitter assumes (0,1)), so I have to sometimes adjust accordingly
	uint64_t S = 0;
	double alphaInv = 13.0;
	if (n >= N/alphaInv) { // if the threshold is not satisfied, use Vitter's A algorithm
		return this->vitterA(n, N);
	} else if (n == 1){ // trivial case
		do {
			S = static_cast<uint64_t>(floor(N * this->runif()));
			
		} while (S == static_cast<uint64_t>(N)); // s == N would be a problem because it would overshoot the end of the array/file (in base-0 space); effectively sampling [0,1) uniform to prevent this
		return S;

	}
	
	// if we pass all thresholds, we use Vitter's rejection scheme
	double nInv     = 1.0/n;
	double nMin1inv = 1.0/(n - 1.0);
	double qu1db    = 1.0 + N - n;
	uint64_t qu1    = static_cast<uint64_t>(qu1db);
	double Vprime;
	double X;
	double Sdb;
	double U;
	double y1;
	
	// outer loop in Vitter (1987) A2
	while (1) {
		
		unsigned int d2tst = 0;
		do {  // step D2; generate U and X
			Vprime = pow(this->runif(), nInv);
			X      = N * (1.0 - Vprime);
			S      = floor(X);
			d2tst++;
		} while (S >= qu1);
		
		U      = this->runif();
		Sdb    = static_cast<double>(S);
		y1     = pow(U*N/qu1db, nMin1inv);
		Vprime = y1 * (1.0 - X/N) * (qu1db/(qu1db - Sdb));
		if (Vprime < 1.0) { // Step D3: accept test 2.8 (Vitter 1987)
			break;
		}
		// moving on to Step D4
		double y2  = 1.0;
		double top = N - 1.0;
		double bottom;
		double limit;
		if (Sdb < n - 1.0) {
			bottom = N - n;
			limit  = N - Sdb;
		} else {
			bottom = N - Sdb - 1.0;
			limit  = qu1db;
		}
		
		// calculate f(|_X_|)
		for (double t = N - 1.0; t >= limit; t--) {
			y2 = y2 * top/bottom;
			top--;
			bottom--;
		}
		if (N/(N - X) >= y1 * pow(y2, nMin1inv)) { // Accept D4 condition
			break;
		}
		// reject everything, go back to the start
	}
	return S;
}



