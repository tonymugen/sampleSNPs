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
 * Definitions and interface documentation for the class that relates individual lines to populations they belong to.
 *
 */

#ifndef populations_hpp
#define populations_hpp

#include <vector>
#include <string>

using std::vector;
using std::string;

namespace sampFiles {
	/** \brief Population index
	 *
	 * For each population, contains indexes of the lines that belong to it.
	 */
	class PopIndex {
		
	private:
		/** \brief Vector of index vectors
		 *
		 * The outside vector is the same length as the number of populations. Each inside vector has the line indexes.
		 */
		vector< vector<size_t> > _popInd;
		/** \brief Sample size
		 *
		 * Total number of lines represented.
		 */
		size_t _N;
		
	public:
		/** \brief Default constructor */
		PopIndex() {};
		/** \brief Array constructor
		 *
		 * The input array has an element for each line, and the value of that element is the population ID in the form of an _int_ that is base-1 (i.e., if line N is in the first population, then `arr[N] == 1`).
		 *
		 * \param[in] arr array of population IDs
		 * \param[in] N array length
		 */
		PopIndex(const int *arr, const size_t &N);
		/** \brief File read constructor
		 *
		 * The input file has an entry for each line (separated by white space), and the value of that entry is the population ID in the form of an _int_ that is base-1 (i.e., if line N is in the first population, then `arr[N] == 1`).
		 *
		 * \param[in] inFileName input file name
		 */
		PopIndex(const string &inFileName);
		
		~PopIndex(){};
		
		/** \brief Vector subscript operator
		 *
		 * Returns the index of population _i_.
		 *
		 * \param[in] i population index
		 * \return index of line IDs
		 */
		vector<size_t> & operator[] (const size_t &i) { return _popInd[i]; };
		/** \brief `const` vector subscript operator
		 *
		 * Returns the index of population _i_.
		 *
		 * \param[in] i population index
		 * \return index of line IDs
		 */
		const vector<size_t> & operator[] (const size_t &i) const { return _popInd[i]; };
		
		/** \brief Population size
		 *
		 * \param[in] i population index
		 * \return size of the _i_th population
		 */
		size_t popSize(const size_t &i) {return _popInd[i].size(); };
		/** \brief `const` population size
		 *
		 * \param[in] i population index
		 * \return size of the _i_th population
		 */
		size_t popSize(const size_t &i) const {return _popInd[i].size(); };
		
		/** \brief Total sample size
		 *
		 * \return total sample size
		 */
		size_t size() {return _N; };
		/** \brief `const` total sample size
		 *
		 * \return total sample size
		 */
		size_t size() const {return _N; };
		
		/** \brief Number of populations
		 *
		 * \return number of populations
		 */
		size_t popNumber() {return _popInd.size(); };
		/** \brief `const` number of populations
		 *
		 * \return number of populations
		 */
		size_t popNumber() const {return _popInd.size(); };
		
	};

}

#endif /* populations_hpp */
