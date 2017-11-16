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
 * Definitions and interface documentation for classes that read and write various genetic variant file formats.
 *
 * Currently supported formats:
 *	- _plink_ BED
 *	- _plink_ TPED
 *	- VCF
 *	- Hapmap (.hmp.txt)
 *
 *
 */

#ifndef varfiles_hpp
#define varfiles_hpp

#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <limits>

#include "populations.hpp"

using std::fstream;
using std::string;
using std::vector;
using std::unordered_map;
using std::numeric_limits;

namespace sampFiles {
	// Forward declarations
	class VarFile;
    class GbinFile;
    class GbinFileI;
    class GbinFileO;
	class BedFile;
	class BedFileI;
	class BedFileO;
    class GtxtFile;
    class GtxtFileI;
    class GtxtFileO;
	class TpedFile;
	class TpedFileI;
	class TpedFileO;
	class VcfFile;
	class VcfFileI;
	class VcfFileO;
	class HmpFile;
	class HmpFileI;
	class HmpFileO;

	/** \brief Buffer size
	 *
	 * Size of the buffer for reading files text files. I use it in functions that count the number of lines, for example. The buffer size (10M) is optimized for a MacBook Pro with an SSD. Other systems may perform better with a different value (e.g., if you have a spinning drive and more RAM you may want to experiemtn with increasing it).
	 *
	 */
	static const size_t BUF_SIZE = 10485760;
	/// Machine \f$\epsilon\f$
	const double EPS = numeric_limits<double>::epsilon();
	/// pi
	const double PI = 3.14159265358979323846264338328;
	
	/** \brief Base variant file class
	 *
	 * Abstract base class for all the input/output formats. Cannot be initialized directly.
	 */
	class VarFile {
	protected:
		/// Variant file stream
		fstream _varFile;
		
		/// Default constructor (protected)
		VarFile() {_varFile.exceptions(fstream::badbit); };
		
	public:
		/// Copy constructor
		VarFile(const VarFile &in) = default;
		/// Copy assignment
		VarFile &operator=(const VarFile &in) = default;
		/// Move constructor
		VarFile(VarFile &&in) = default;
		/// Move assignment
		VarFile &operator=(VarFile &&in) = default;
		/// Destructor
		~VarFile(){if (_varFile.is_open()) _varFile.close(); };
		
		/// \brief Open stream
		virtual void open() = 0;
		/// Close stream
		virtual void close() = 0;
	};
	
    /** \brief Generic binary file base class
     *
     * Sets up streams for binary files. No support for header lines.
     *
     */
    class GbinFile : public VarFile {
    protected:
        /// File name
        string _fileName;
        /// Number of elements in a row
        size_t _nCols;
        /// Size of each element in bytes
        size_t _elemSize;
        
    public:
        /// Default constructor
        GbinFile() : VarFile(), _nCols{0}, _elemSize{sizeof(char)} {};
        /** \brief Constructor with file name
         *
         * Throws ``Number of elements not divisible by the number of columns'' if the total number of elements is not divisible by the number of elements in a column.
         *
         * \param[in] fileName file name
         * \param[in] nCols number of columns, or elements in a row
         * \param[in] elemSize size of each element in bytes
         *
         */
        GbinFile(const string &fileName, const size_t &nCols, const size_t &elemSize) : VarFile(), _fileName{fileName}, _nCols{nCols}, _elemSize{elemSize} {};
        
        /// Copy constructor
        GbinFile(const GbinFile &in) = default;
        /// Copy assignment
        GbinFile &operator=(const GbinFile &in) = default;
        /// Move constructor
        GbinFile(GbinFile &&in) = default;
        /// Move assignment
        GbinFile &operator=(GbinFile &&in) = default;
        /// Destructor
        ~GbinFile(){};
        
        /// \brief Open stream (does nothing)
        virtual void open() {};
        /// Close stream
        virtual void close();
        
    };
    
    /** \brief Binary file input class
     *
     * Reads binary files.
     *
     */
    class GbinFileI : GbinFile {
    protected:
        /** \brief Get number of rows in the binary file
         *
         * Requires knowledge of the number of elements in a row and their size in bytes. Throws a _string_ object ``Number of elements not divisible by row size'' if the total number of elements in the file is not divisible by the product of the number of columns and element size.
         *
         * \return number of rows
         */
        virtual uint64_t _numLines();

    public:
        /// Default constructor
        GbinFileI() : GbinFile() {};
        /** \brief File name constructor
         *
         * \param[in] fileName file name including extension
         * \param[in] nCols number of columns, or elements in a row
         * \param[in] elemSize size of each element in bytes
         *
         */
        GbinFileI(const string &fileName, const size_t &nCols, const size_t &elemSize) : GbinFile(fileName, nCols, elemSize) {};
        /// Copy constructor
        GbinFileI(const GbinFileI &in) = default;
        /// Copy assignment
        GbinFileI &operator=(const GbinFileI &in) = default;
        /// Move constructor
        GbinFileI(GbinFileI &&in) = default;
        /// Move assignment
        GbinFileI &operator=(GbinFileI &&in) = default;
        /// Destructor
        ~GbinFileI(){};
        
        /// \brief Open stream to read
        void open();
        
        /** \brief Sample rows and save to a binary file
         *
         * Sample \f$n\f$ lines without replacement from the file represented by the current object and save to the `out` object. Uses Vitter's \cite vitter87a method. Number of samples has to be smaller that the number of rows in the file.
         *
         * \param[in] out output object
         * \param[in] n number of rows to sample
         *
         */
        void sample(GbinFileO &out, const uint64_t &n);
        /// Number of rows in the object
        uint64_t nlines() { return _numLines(); };
        
    };
    
    /** \brief Generic binary file output class
     *
     * Writes binary files.
     *
     */
    class GbinFileO : public GbinFile {
        friend class GbinFileI;
    protected:
        
    public:
        /// Default constructor
        GbinFileO() : GbinFile() {};
        /** \brief File name constructor
         *
         * \param[in] fileName file name including extension
         * \param[in] nCols number of columns, or elements in a row
         * \param[in] elemSize size of each element in bytes
         *
         */
        GbinFileO(const string &fileName, const size_t &nCols, const size_t &elemSize) : GbinFile(fileName, nCols, elemSize) {};
        /// Copy constructor
        GbinFileO(const GbinFileO &in) = default;
        /// Copy assignment
        GbinFileO &operator=(const GbinFileO &in) = default;
        /// Move constructor
        GbinFileO(GbinFileO &&in) = default;
        /// Move assignment
        GbinFileO &operator=(GbinFileO &&in) = default;
        /// Destructor
        ~GbinFileO(){};
        
        /// \brief Open stream to write
        void open();
        
    };
    /** \brief BED file base class
     *
     * Sets up streams for the auxiliary files.
     */
    class BedFile : public GbinFile {
    protected:
        
        /// Corresponding .fam file stream
        fstream _famFile;
        /// Corresponding .bim file stream
        fstream _bimFile;
        /// File name stub (minus the extension)
        string _fileStub;
        /** \brief Genotype bit masks
         *
         * Used to isolate each of the four genotypes (moving from the last bit pair) from the .bed two-bit genotype coding.
         *
         */
        static const vector<char> _masks;
        /** \brief Genotype bit tests
         *
         * Used to test each of the four genotypes (moving from the last bit pair) from the .bed two-bit genotype coding for the possible states.
         *
         */
        static const unordered_map<char, string> _tests;

        
    public:
        /// Default constructor
        BedFile();
        /** \brief File name constructor
         *
         * \param[in] stubName file name minus the extension
         */
        BedFile(const string &stubName);
        /// Copy constructor
        BedFile(const BedFile &in) = default;
        /// Copy assignment
        BedFile &operator=(const BedFile &in) = default;
        /// Move constructor
        BedFile(BedFile &&in) = default;
        /// Move assignment
        BedFile &operator=(BedFile &&in) = default;
        /// Destructor
        ~BedFile();
        
        /// Open stream (does nothing)
        virtual void open() {};
        /// Close stream
        void close();
        
    };

	/** \brief BED file input class
	 *
	 * Reads BED files and the auxiliary files that come with them (.fam and .bim) as necessary. Only the SNP-major version is supported.
	 *
	 */
	class BedFileI : public BedFile {
	protected:
		
        /** \brief Get number of lines in the `_bimFile`
         *
         * Assumes Unix-like line endings. The result is equal to the number of SNPs.
         *
         * \return number of lines in `_bimFile`
         */
        uint64_t _numLines();
        /** \brief Get number of lines in the `_famFile`
         *
         * Assumes Unix-like line endings. The result is equal to the number of individuals.
         *
         * \return number of lines in `_famFile`
         */
        uint64_t _famLines();
        /** \brief Copy the .fam file and count number of lines
         *
         * Assumes Unix-like line endings. The result is equal to the number of individuals. The current object's .fam file is copied to the provided file stream, which should be open for raading. If not, the function throws a _string_ object ``Output .fam filestream not open''.
         *
         * \param[in] fam .fam file stream
         *
         * \return number of lines in `_famFile`
         */
        uint64_t _famLines(fstream &fam);

		/** \brief Between-SNP linkage disequilibrium (LD)
		 *
		 * Calculates two LD statistics (\f$ r^2 \f$ and \f$ D' \f$)  between two SNPs from a BED file. Missing values are ignored. If there are fewer than three haplotypes with data present at both loci, the return values are -9. This value is also returned if one of the loci is monomorphic after taking out missing data at the other SNP.
		 * Minor (not necessarily derived) allele counts are also reported to enable downstream filtering. Note that the populations are assumed diploid and the counts are of haploid chromosomes (i.e. one homozygote yields count of 2).
		 *
		 * \param[in] snp1 first SNP
		 * \param[in] snp2 second SNP
		 * \param[in] N length of the genotype vector in bytes (four genotypes per byte)
		 * \param[in] pad number of bit pairs of padding in the last byte
		 * \param[out] rSq the \f$ r^2 \f$ estimate
		 * \param[out] Dprime the \f$ D' \f$ estimate
		 * \param[out] dcnt1 minor allele count at locus 1
		 * \param[out] dcnt2 minor allele count at locus 2
		 *
		 */
		void _ld(const char *snp1, const char *snp2, const size_t &N, const unsigned short &pad, double &rSq, double &Dprime, double &dcnt1, double &dcnt2);
		/** \brief Between-SNP LD within populations
		 *
		 * Calculates two LD statistics (\f$ r^2 \f$ and \f$ D' \f$)  between two SNPs from a BED file. Missing values are ignored. If there are fewer than three haplotypes with data present at both loci, the return values are -9. This value is also returned if one of the loci is monomorphic after taking out missing data at the other SNP.
		 * Minor (not necessarily derived) allele counts are also reported to enable downstream filtering. Note that the populations are assumed diploid and the counts are of haploid chromosomes (i.e. one homozygote yields count of 2). The values are calculted within each population as indicated by the `PopIndex` object. The results are returned in the supplied vectors, which are assumed to be of correct size. Since this is an internal function unexposed to the user, this is not chaecked to save on compuation steps.
		 * Care must be taken that the `char` arrays passed to the function have lengths compatible with the number of individuals indexed by `PopIndex`. This is not checked.
		 *
		 * \param[in] snp1 first SNP
		 * \param[in] snp2 second SNP
		 * \param[in] popID population index
		 * \param[out] rSq vector of \f$ r^2 \f$ estimates
		 * \param[out] Dprime vector of \f$ D' \f$ estimates
		 * \param[out] dcnt1 vector of minor allele counts at locus 1
		 * \param[out] dcnt2 vector of minor allele counts at locus 2
		 *
		 */
		void _ld(const char *snp1, const char *snp2, const PopIndex &popID, vector<double> &rSq, vector<double> &Dprime, vector<double> &dcnt1, vector<double> &dcnt2);
	public:
		/// Default constructor
		BedFileI() : BedFile() {};
		/** \brief File name constructor
		 *
		 * \param[in] stubName file name minus the extension
		 */
		BedFileI(const string &stubName) : BedFile(stubName) {};
		/// Copy constructor
		BedFileI(const BedFileI &in) = default;
		/// Copy assignment
		BedFileI &operator=(const BedFileI &in) = default;
		/// Move constructor
		BedFileI(BedFileI &&in) = default;
		/// Move assignment
		BedFileI &operator=(BedFileI &&in) = default;
		/// Destructor
		~BedFileI(){};
		
		/// \brief Open stream to read
		void open();
		
		/** \brief Sample SNPs and save to BED file
		 *
		 * Sample \f$n\f$ SNPs without replacement from the file represented by the current object and save to the `out` object. Uses Vitter's \cite vitter87a method. Number of samples has to be smaller that the number of SNPs in the file.
		 *
		 * \param[in] out output object
		 * \param[in] n number of SNPs to sample
		 *
		 */
		void sample(BedFileO &out, const uint64_t &n);
		
		/** \brief Linkage disequilibrium among sampled sites
		 *
		 * Samples sequential pairs of SNPs and calculates two LD measures (\f$ r^2 \f$  and \f$ D' \f$). Saves to a file with the same name as the one preceding the .bed etc extensions, but adds _LD.tsv at the end. Each line is tab-delimited with the chromosome number (from the .bim file), between-SNP distance, non-reference allele count for each SNP, \f$ r^2 \f$, and \f$ D' \f$.
		 * Missing data are ignored (only pairwise-complete observations are included). If one of the SNPs is monomorphic or if the total number of pairwise present genotypes is fewer than three (exclusive), the LD measures are returned as -9 to indicate missing values.
		 *
		 * \param[in] n number of SNP pairs to sample
		 */
		void sampleLD(const uint64_t &n);
		/** \brief LD among sampled sites within populations
		 *
		 * Samples sequential pairs of SNPs and calculates two LD measures (\f$ r^2 \f$  and \f$ D' \f$) within populations indicated by `PopIndex`. Saves to a file with the same name as the one preceding the .bed etc extensions, but adds _LD.tsv at the end. Each line is tab-delimited with the chromosome number (from the .bim file), between-SNP distance, non-reference allele count for each SNP, \f$ r^2 \f$, and \f$ D' \f$.
		 * Missing data are ignored (only pairwise-complete observations are included). If one of the SNPs is monomorphic or if the total number of pairwise present genotypes is fewer than three (exclusive), the LD measures are returned as -9 to indicate missing values.
		 *
		 * \param[in] popID population index
		 * \param[in] n number of SNP pairs to sample
		 */
		void sampleLD(const PopIndex &popID, const uint64_t &n);
		/// Number of SNPs in the object
		uint64_t nsnp() { return _numLines(); };
		/// Number of individuals in the object
		uint64_t nindiv() { return _famLines(); };
		
	};
	
	/** \brief BED file output class
	 *
	 * Writes to BED files and the auxiliary files that come with them (.fam and .bim) as necessary. Data are written in the SNP-major format.
	 *
	 */
	class BedFileO : public BedFile {
		friend class BedFileI;
	protected:
		
	public:
		/// Default constructor
		BedFileO() : BedFile() {};
		/** \brief File name constructor
		 *
		 * \param[in] stubName file name minus the extension
		 */
		BedFileO(const string &stubName) : BedFile(stubName) {};
		/// Copy constructor
		BedFileO(const BedFileO &in) = default;
		/// Copy assignment
		BedFileO &operator=(const BedFileO &in) = default;
		/// Move constructor
		BedFileO(BedFileO &&in) = default;
		/// Move assignment
		BedFileO &operator=(BedFileO &&in) = default;
		/// Destructor
		~BedFileO(){};
		
		/// \brief Open stream to write
		void open();
		
	};
    /** \brief Generic text file base class
     *
     * Sets up streams for text files. If specified, the first line can be considered a header and treated separately.
     *
     */
    class GtxtFile : public VarFile {
    protected:
        /// File name
        string _fileName;
        /// Is there a header?
        bool _head;
        
    public:
        /// Default constructor
        GtxtFile() : VarFile(), _head{false} {};
        /** \brief Constructor with file name
         *
         * \param[in] fileName file name
         *
         */
        GtxtFile(const string &fileName) : VarFile(), _fileName{fileName}, _head{false} {};
        /** \brief Constructor with file name and header indicator
         *
         * \param[in] fileName file name
         * \param[in] head header presence
         *
         */
        GtxtFile(const string &fileName, const bool &head) : VarFile(), _fileName{fileName}, _head{head} {};
        
        /// Copy constructor
        GtxtFile(const GtxtFile &in) = default;
        /// Copy assignment
        GtxtFile &operator=(const GtxtFile &in) = default;
        /// Move constructor
        GtxtFile(GtxtFile &&in) = default;
        /// Move assignment
        GtxtFile &operator=(GtxtFile &&in) = default;
        /// Destructor
        ~GtxtFile(){};
        
        /// \brief Open stream (does nothing)
        virtual void open() {};
        /// Close stream
        virtual void close();
        
    };
    
    /** \brief Text file input class
     *
     * Reads text files, skipping or copying the header as necessary.
     *
     */
    class GtxtFileI : GtxtFile {
    protected:
        /** \brief Get number of rows in the text file.
         *
         * Assumes Unix-like line endings. Header, if present, is not counted. Is overriden in some, but not all, derived classes.
         *
         * \return number of rows
         */
        virtual uint64_t _numLines();

    public:
        /// Default constructor
        GtxtFileI() : GtxtFile() {};
        /** \brief File name constructor with header specification
         *
         * \param[in] fileName file name including extension
         */
        GtxtFileI(const string &fileName) : GtxtFile(fileName) {};
        /** \brief File name constructor with header specification
         *
         * \param[in] fileName file name including extension
         * \param[in] head header presence
         */
        GtxtFileI(const string &fileName, const bool &head) : GtxtFile(fileName, head) {};
        /// Copy constructor
        GtxtFileI(const GtxtFileI &in) = default;
        /// Copy assignment
        GtxtFileI &operator=(const GtxtFileI &in) = default;
        /// Move constructor
        GtxtFileI(GtxtFileI &&in) = default;
        /// Move assignment
        GtxtFileI &operator=(GtxtFileI &&in) = default;
        /// Destructor
        ~GtxtFileI(){};
        
        /// \brief Open stream to read
        void open();
        
        /** \brief Sample rows and save to a text file
         *
         * Sample \f$n\f$ lines without replacement from the file represented by the current object and save to the `out` object. Uses Vitter's \cite vitter87a method. Number of samples has to be smaller that the number of rows in the file.
         *
         * \param[in] out output object
         * \param[in] n number of rows to sample
         * \param[in] headSkip skip header? Ignored if there is no header
         *
         */
        void sample(GtxtFileO &out, const uint64_t &n, const bool &headSkip);
        /** \brief Sample rows and save export to a vector of strings
         *
         * Sample \f$n\f$ rows without replacement from the file represented by the current object and output a vector of strings. Each field separated by the specified delimiter is stored as an element of the vector. Uses Vitter's \cite vitter87a method. Number of samples has to be smaller that the number of rows in the file. The output vector is erased if it is not empty.
         *
         * \param[in] n number of SNPs to sample
         * \param[in] headSkip skip header? Ignored if there is no header
         * \param[in] delim field delimiter
         * \param[out] out output vector
         *
         */
        void sample(const uint64_t &n, const bool &headSkip, const char &delim, vector<string> &out);
        /// Number of SNPs in the object
        uint64_t nlines() { return _numLines(); };
        
    };
    
    /** \brief Generic text file output class
     *
     * Writes text files.
     *
     */
    class GtxtFileO : public GtxtFile {
        friend class GtxtFileI;
    protected:
        
    public:
        /// Default constructor
        GtxtFileO() : GtxtFile() {};
        /** \brief File name constructor
         *
         * \param[in] fileName file name including the extension
         */
        GtxtFileO(const string &fileName) : GtxtFile(fileName) {};
        /** \brief File name constructor with header specification
         *
         * \param[in] fileName file name including the extension
         * \param[in] head header presence
         */
        GtxtFileO(const string &fileName, const bool &head) : GtxtFile(fileName, head) {};
        /// Copy constructor
        GtxtFileO(const GtxtFileO &in) = default;
        /// Copy assignment
        GtxtFileO &operator=(const GtxtFileO &in) = default;
        /// Move constructor
        GtxtFileO(GtxtFileO &&in) = default;
        /// Move assignment
        GtxtFileO &operator=(GtxtFileO &&in) = default;
        /// Destructor
        ~GtxtFileO(){};
        
        /// \brief Open stream to write
        void open();
        
    };
    

	/** \brief TPED file base class
	 *
	 * Sets up the stream for the corresponding .tfam file.
	 */
	class TpedFile : public GtxtFile {
	protected:
		/// Corresponding .tfam file stream
		fstream _tfamFile;
		/// File name stub (minus the extension)
		string  _fileStub;
		
	public:
		/// Default constructor
		TpedFile() : GtxtFile() {_tfamFile.exceptions(fstream::badbit); };
		/** \brief File name constructor
		 *
		 * \param[in] stubName file name minus the extension
		 */
        TpedFile(const string &stubName) : GtxtFile(stubName + ".tped"), _fileStub{stubName} {_tfamFile.exceptions(fstream::badbit); }; // no headers in .tped
		/// Copy constructor
		TpedFile(const TpedFile &in) = default;
		/// Copy assignment
		TpedFile &operator=(const TpedFile &in) = default;
		/// Move constructor
		TpedFile(TpedFile &&in) = default;
		/// Move assignment
		TpedFile &operator=(TpedFile &&in) = default;
		/// Destructor
		~TpedFile();
		
		/// Open stream (does nothing)
		virtual void open() {};
		/// Close stream
		void close();
		
	};
	
	/** \brief TPED file input class
	 *
	 * Reads TPED files and the corresponding .tfam files as necessary.
	 *
	 */
	class TpedFileI : public TpedFile {
	protected:
        /** \brief Get number of lines in the `_tfamFile`
         *
         * Assumes Unix-like line endings. The result is equal to the number of individuals. The `_tfamFile` should already be open for reading.
         *
         * \return number of lines in `_tfamFile`
         */
        uint64_t _famLines();
        /** \brief Copy the .tfam file and count number of lines
         *
         * Assumes Unix-like line endings. The result is equal to the number of individuals. The current object's .tfam file is copied to the provided file stream, which should already be open for writing. If not, the function throws a _string_ object ``Output .fam filestream not open''.
         *
         * \param[in] fam .tfam file stream
         *
         * \return number of lines in `_tfamFile`
         */
        uint64_t _famLines(fstream &fam);
        /** \brief Copy the .tfam file
         *
         * The current object's .tfam file is copied to the provided file stream, which should already be open for writing. If not, the function throws a _string_ object ``Output .fam filestream not open''.
         *
         * \param[in] fam .tfam file stream
         *
         */
        void _famCopy(fstream &fam);
        /** \brief Get number of rows in the text file.
         *
         * Assumes Unix-like line endings. Header, if present, is not counted. Is overriden in some, but not all, derived classes.
         *
         * \return number of rows
         */
        uint64_t _numLines();

	public:
		/// Default constructor
		TpedFileI() : TpedFile() {};
		/** \brief File name constructor
		 *
		 * \param[in] stubName file name minus the extension
		 */
		TpedFileI(const string &stubName) : TpedFile(stubName) {};
		/// Copy constructor
		TpedFileI(const TpedFileI &in) = default;
		/// Copy assignment
		TpedFileI &operator=(const TpedFileI &in) = default;
		/// Move constructor
		TpedFileI(TpedFileI &&in) = default;
		/// Move assignment
		TpedFileI &operator=(TpedFileI &&in) = default;
		/// Destructor
		~TpedFileI(){};
		
		/// \brief Open stream to read
		void open();
		
		/** \brief Sample SNPs and save to BED file
		 *
		 * Sample \f$n\f$ SNPs without replacement from the file represented by the current object and save to the `out` object. Uses Vitter's \cite vitter87a method. Number of samples has to be smaller that the number of SNPs in the file.
		 *
		 * \param[in] out output object
		 * \param[in] n number of SNPs to sample
		 *
		 */
		void sample(TpedFileO &out, const uint64_t &n);
		/// Number of SNPs in the object
        uint64_t nsnp() { return _numLines(); };
		/// Number of individuals in the object
		uint64_t nindiv() { return _famLines(); };
		
	};
	
	/** \brief TPED file output class
	 *
	 * Writes to TPED files and the corresponding .tfam files as necessary. Data are written in the SNP-major format.
	 *
	 */
	class TpedFileO : TpedFile {
		friend class TpedFileI;
	protected:
		
	public:
		/// Default constructor
		TpedFileO() : TpedFile() {};
		/** \brief File name constructor
		 *
		 * \param[in] stubName file name minus the extension
		 */
		TpedFileO(const string &stubName) : TpedFile(stubName) {};
		/// Copy constructor
		TpedFileO(const TpedFileO &in) = default;
		/// Copy assignment
		TpedFileO &operator=(const TpedFileO &in) = default;
		/// Move constructor
		TpedFileO(TpedFileO &&in) = default;
		/// Move assignment
		TpedFileO &operator=(TpedFileO &&in) = default;
		/// Destructor
		~TpedFileO(){};
		
		/// \brief Open stream to write
		void open();
	};
	
	/** \brief VCF file base class
	 *
	 * Sets up streams for VCF files. Any accompanying .idx files are ignored.
	 *
	 */
	class VcfFile : public GtxtFile {
	protected:
	
    public:
		/// Default constructor
        VcfFile() : GtxtFile() {};
		/** \brief Constructor with file name
		 *
		 * \param[in] fileName file name
		 *
		 */
		VcfFile(const string &fileName) : GtxtFile(fileName) {};
		
		/// Copy constructor
		VcfFile(const VcfFile &in) = default;
		/// Copy assignment
		VcfFile &operator=(const VcfFile &in) = default;
		/// Move constructor
		VcfFile(VcfFile &&in) = default;
		/// Move assignment
		VcfFile &operator=(VcfFile &&in) = default;
		/// Destructor
		~VcfFile(){};
		
		/// \brief Open stream (does nothing)
		void open() {};
		/// Close stream
		void close();
	};
	
	/** \brief VCF file input class
	 *
	 * Reads VCF files, skipping or copying the header as necessary; .idx files are ignored.
	 *
	 */
	class VcfFileI : public VcfFile {
	protected:
        /** \brief Get number of SNPs in the VCF file.
         *
         * Assumes Unix-like line endings. Header is not counted.
         *
         * \return number of SNPs
         */
        uint64_t _numLines();

	public:
		/// Default constructor
		VcfFileI() : VcfFile() {};
		/** \brief File name constructor
		 *
		 * \param[in] fileName file name including extension
		 */
		VcfFileI(const string &fileName) : VcfFile(fileName) {};
		/// Copy constructor
		VcfFileI(const VcfFileI &in) = default;
		/// Copy assignment
		VcfFileI &operator=(const VcfFileI &in) = default;
		/// Move constructor
		VcfFileI(VcfFileI &&in) = default;
		/// Move assignment
		VcfFileI &operator=(VcfFileI &&in) = default;
		/// Destructor
		~VcfFileI(){};
		
		/// \brief Open stream to read
		void open();
		
		/** \brief Sample SNPs and save to VCF file
		 *
		 * Sample \f$n\f$ SNPs without replacement from the file represented by the current object and save to the `out` object. Uses Vitter's \cite vitter87a method. Number of samples has to be smaller that the number of SNPs in the file.
		 *
		 * \param[in] out output object
		 * \param[in] n number of SNPs to sample
		 *
		 */
		void sample(VcfFileO &out, const uint64_t &n);
		/// Number of SNPs in the object
		uint64_t nsnp() { return _numLines(); };
		
	};
	
	/** \brief VCF file output class
	 *
	 * Writes VCF files.
	 *
	 */
	class VcfFileO : public VcfFile {
		friend class VcfFileI;
	protected:
		
	public:
		/// Default constructor
		VcfFileO() : VcfFile() {};
		/** \brief File name constructor
		 *
		 * \param[in] fileName file name including the extension
		 */
		VcfFileO(const string &fileName) : VcfFile(fileName) {};
		/// Copy constructor
		VcfFileO(const VcfFileO &in) = default;
		/// Copy assignment
		VcfFileO &operator=(const VcfFileO &in) = default;
		/// Move constructor
		VcfFileO(VcfFileO &&in) = default;
		/// Move assignment
		VcfFileO &operator=(VcfFileO &&in) = default;
		/// Destructor
		~VcfFileO(){};
		
		/// \brief Open stream to write
		void open();
		
	};
	
	/** \brief Hapmap (HMP) file base class
	 *
	 * Sets up streams for HMP files.
	 *
	 */
	class HmpFile : public GtxtFile {
	protected:
        
	public:
		/// Default constructor
		HmpFile() : GtxtFile() {};
		/** \brief Constructor with file name
		 *
		 * \param[in] fileName file name
		 *
		 */
        HmpFile(const string &fileName) : GtxtFile(fileName) {};
		
		/// Copy constructor
		HmpFile(const HmpFile &in) = default;
		/// Copy assignment
		HmpFile &operator=(const HmpFile &in) = default;
		/// Move constructor
		HmpFile(HmpFile &&in) = default;
		/// Move assignment
		HmpFile &operator=(HmpFile &&in) = default;
		/// Destructor
		~HmpFile(){};
		
		/// \brief Open stream (does nothing)
		virtual void open() {};
		/// Close stream
		virtual void close();
		
	};
	
	/** \brief HMP file input class
	 *
	 * Reads HMP files, skipping or copying the header as necessary.
	 *
	 */
	class HmpFileI : HmpFile {
	protected:
        /** \brief Get number of SNPs in the HMP file.
         *
         * Assumes Unix-like line endings. Header is not counted.
         *
         * \return number of SNPs
         */
        uint64_t _numLines();

	public:
		/// Default constructor
		HmpFileI() : HmpFile() {};
		/** \brief File name constructor
		 *
		 * \param[in] fileName file name including extension
		 */
        HmpFileI(const string &fileName);
		/// Copy constructor
		HmpFileI(const HmpFileI &in) = default;
		/// Copy assignment
		HmpFileI &operator=(const HmpFileI &in) = default;
		/// Move constructor
		HmpFileI(HmpFileI &&in) = default;
		/// Move assignment
		HmpFileI &operator=(HmpFileI &&in) = default;
		/// Destructor
		~HmpFileI(){};
		
		/// \brief Open stream to read
		void open();
		
		/** \brief Sample SNPs and save to HMP file
		 *
		 * Sample \f$n\f$ SNPs without replacement from the file represented by the current object and save to the `out` object. Uses Vitter's \cite vitter87a method. Number of samples has to be smaller that the number of SNPs in the file.
		 *
		 * \param[in] out output object
		 * \param[in] n number of SNPs to sample
		 *
		 */
		void sample(HmpFileO &out, const uint64_t &n);
		/// Number of SNPs in the object
		uint64_t nsnp() { return _numLines(); };
		
	};
	
	/** \brief HMP file output class
	 *
	 * Writes HMP files.
	 *
	 */
	class HmpFileO : public HmpFile {
		friend class HmpFileI;
	protected:
		
	public:
		/// Default constructor
		HmpFileO() : HmpFile() {};
		/** \brief File name constructor
		 *
		 * \param[in] fileName file name including the extension
		 */
		HmpFileO(const string &fileName) : HmpFile(fileName) {};
		/// Copy constructor
		HmpFileO(const HmpFileO &in) = default;
		/// Copy assignment
		HmpFileO &operator=(const HmpFileO &in) = default;
		/// Move constructor
		HmpFileO(HmpFileO &&in) = default;
		/// Move assignment
		HmpFileO &operator=(HmpFileO &&in) = default;
		/// Destructor
		~HmpFileO(){};
		
		/// \brief Open stream to write
		void open();
		
	};
    

}

#endif /* varfiles_hpp */




