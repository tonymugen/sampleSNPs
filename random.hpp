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
 * Class definition and interface documentation for facilities that generate random draws from various distributions.
 *
 */

#ifndef random_hpp
#define random_hpp

#include <cstdint>

namespace sampFiles {
	class RanDraw;
	class Generate;
	class GenerateHR;
	class GenerateMT;
	
	/** \brief Abstract base random number class
	 *
	 * Provides the interface for random or pseudorandom (depending on derived class) generation. For internal use by the `RanDraw` interface class.
	 */
	class Generate {
	protected:
		/** \brief Protected default constructor */
		Generate(){};
		/** \brief Protected copy constructor
		 *
		 * \param[in] old object to copy
		 */
		Generate(const Generate &old){};
		/** \brief Protected move constructor
		 *
		 * \param[in] old object to move
		 */
		Generate(Generate &&old){};
		/** \brief Protected copy assignment operator
		 *
		 * \param[in] old object to copy
		 */
		Generate & operator= (const Generate &old) = default;
		/** \brief Protected move assignment
		 *
		 * \param[in] old object to move
		 */
		Generate & operator= (Generate &&old) = default;

	public:
		/** \brief Protected destructor */
		virtual ~Generate(){};
		/** \brief Generate a (pseudo-)random 64-bit unsigned integer
		 *
		 * \return random or pseudo-random 64-bit unsigned integer
		 */
		virtual volatile uint64_t ranInt() = 0;
	};
	
	/** \brief Hardware random number generating class
	 *
	 * Generates random deviates from a number of distributions, using hardware random numbers (_RDRAND_ processor instruction). Health of the RDRAND generator is tested every time a new number is required. Throws a `string` object "RDRAND_failed" if the test fails.
	 * The implementation of random 64-bit integer generation follows [Intel's suggestions](https://software.intel.com/en-us/articles/intel-digital-random-number-generator-drng-software-implementation-guide ).
	 */
	class GenerateHR : public Generate {
	protected:
		// no protected members
	public:
		/** \brief Default constructor */
		GenerateHR(){};
		/** \brief Destructor */
		~GenerateHR(){};
		/** \brief Copy constructor
		 *
		 * \param[in] old object to copy
		 */
		GenerateHR(const GenerateHR &old){};
		/** \brief Move constructor
		 *
		 * \param[in] old object to move
		 */
		GenerateHR(GenerateHR &&old){};
		/** \brief Copy assignment operator
		 *
		 * \param[in] old object to copy
		 */
		GenerateHR & operator= (const GenerateHR &old) = default;
		/** \brief Move assignment
		 *
		 * \param[in] old object to move
		 */
		GenerateHR & operator= (GenerateHR &&old) = default;
		
		/** \brief Generate a random 64-bit unsigned integer
		 *
		 * Monitors the health of the CPU random number generator and throws a `string` object "RDRAND_failed" if a failure is detected after ten tries.
		 *
		 * \return digital random 64-bit unsigned integer
		 */
		volatile uint64_t ranInt();
	};
	
	/** \brief Pseudo-random number generator
	 *
	 * An implementaiton of the 64-bit MT19937 ("Mersenne Twister")  \cite matsumoto98a pseudo-random number generator (PRNG). The constructor automatically seeds the PRNG. The implementation was guided by the reference code [posted by the authors](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt64.html ).
	 */
	class GenerateMT : public Generate {
	protected:
		/** \brief Degree of recurrence */
		static const unsigned short _n;
		/** \brief Middle word */
		static const unsigned short _m;
		/** \brief Most significant 33 bits */
		static const uint64_t _um;
		/** \brief Least significant 31 bits */
		static const uint64_t _lm;
		/** \brief Tempering bitmask */
		static const uint64_t _b;
		/** \brief Tempering bitmask */
		static const uint64_t _c;
		/** \brief Tempering bitmask */
		static const uint64_t _d;
		/** \brief Tempering shift */
		static const unsigned int _l;
		/** \brief Tempering shift */
		static const unsigned int _s;
		/** \brief Tempering shift */
		static const unsigned int _t;
		/** \brief Tempering shift */
		static const unsigned int _u;
		/** \brief Array of alternative values for the twist */
		static const uint64_t _alt[2];
		/** \brief Generator state array */
		uint64_t _mt[312];
		/** \brief State of the array index */
		size_t _mti;
		/** \brief Current state */
		uint64_t _x;
		
	public:
		/** \brief Default constructor 
		 *
		 * Seeds the PRNG with a call to the _RDTSC_ instruction.
		 */
		GenerateMT();
		/** \brief Protected destructor */
		~GenerateMT(){};
		/** \brief Copy constructor
		 *
		 * \param[in] old object to copy
		 */
		GenerateMT(const GenerateMT &old) = default;
		/** \brief Move constructor
		 *
		 * \param[in] old object to move
		 */
		GenerateMT(GenerateMT &&old) = default;
		/** \brief Copy assignment operator
		 *
		 * \param[in] old object to copy
		 */
		GenerateMT & operator= (const GenerateMT &old) = default;
		/** \brief Move assignment
		 *
		 * \param[in] old object to move
		 */
		GenerateMT & operator= (GenerateMT &&old) = default;
		
		/** \brief Generate a pseudo-random 64-bit unsigned integer
		 *
		 *
		 * \return pseudo-random 64-bit unsigned integer
		 */
		volatile uint64_t ranInt();
	};

	/** \brief Random number generating class
	 *
	 * Generates (pseudo-)random deviates from a number of distributions. If hardware random numbers are supported, uses them. Otherwise, falls back to 64-bit MT19937 ("Mersenne Twister").
	 *
	 */
	class RanDraw {
	private:
		Generate *_rng;
		
	public:
		/** \brief Default constructor
		 *
		 * Checks if the processor provides hardware random number support. Seeds the Mersenne Twister if not.
		 * Throws "CPU_unsupported" string object if the CPU is not AMD or Intel.
		 */
		RanDraw();
		
		/** \brief Destructor */
		~RanDraw(){ delete _rng; };
		/** \brief Copy constructor
		 *
		 * \param[in] old pbject to be copied
		 */
		RanDraw(const RanDraw &old) = default;
		/** \brief Move constructor
		 *
		 * \param[in] old pbject to be moved
		 */
		RanDraw(RanDraw &&old) = default;
		
		/** \brief Copy assignment
		 *
		 * \param[in] old pbject to be copied
		 */
		RanDraw & operator= (const RanDraw &old) = default;
		/** \brief Move assignment
		 *
		 * \param[in] old pbject to be moved
		 */
		RanDraw & operator= (RanDraw &&old) = default;
		
		/** \brief Generate random integer
		 *
		 * \return An unsigned random 64-bit integer
		 */
		volatile uint64_t ranInt(){ return _rng->ranInt(); };
		
		/** \brief Generate a uniform deviate
		 *
		 * \return A double-precision value from the \f$ U[0,1]\f$ distribution
		 */
		volatile double runif() {return 5.42101086242752217E-20*static_cast<double>(_rng->ranInt()); };
		/** \brief Generate a non-zero uniform deviate
		 *
		 * \return A double-precision value from the \f$ U(0,1]\f$ distribution
		 */
		volatile double runifnz();
		/** \brief Sample from Vitter's distribution, method A
		 *
		 * Given the number of remaining records in a file \f$N\f$ and the number of records \f$n\f$ remaining to be selected, sample the number of records to skip over. This function implements Vitter's \cite vitter84a \cite vitter87a method A. It is useful for online one-pass sampling of records from a file.
		 * While the inputs are integer, we pass them in as _double_ because that is more efficient for calculations.
		 *
		 * \param[in] n number of records remaining to be picked
		 * \param[in] N number of remaining records in the file
		 *
		 * \return the number of records to skip
		 */
		volatile uint64_t vitterA(const double &n, const double &N);
		/** \brief Sample from Vitter's distribution, method D
		 *
		 * Given the number of remaining records in a file \f$N\f$ and the number of records \f$n\f$ remaining to be selected, sample the number of records to skip over. This function implements Vitter's \cite vitter84a \cite vitter87a method D. It is useful for online one-pass sampling of records from a file.
		 * While the inputs are integer, we pass them in as _double_ because that is more efficient for calculations.
		 *
		 * \param[in] n number of records remaining to be picked
		 * \param[in] N number of remaining records in the file
		 *
		 * \return the number of records to skip
		 */
		volatile uint64_t vitter(const double &n, const double &N);
	};

}

#endif /* random_hpp */
