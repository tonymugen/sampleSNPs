
# Overview

This distribution consists of three parts. One, a program (_sampleSNPs_) that creates and saves ordered random samples of SNPs from a variety of formats. Two, a program (_sampleLD_) that calculates linkage disequilibrium (LD) among randomly chosen (without replacement) pairs of SNPs from a [binary variant format file](http://zzz.bwh.harvard.edu/plink/data.shtml#bed). Three, a library (_libsampFiles.a_) that allows users to build similar applications taking advantage of the fast sampling algorithms used in the two programs mentioned above.

# Requirements

The software is build for Unix-like systems. Neither compilation nor running was checked under Windows and there are reasons to believe it will not compile on that OS. There are no dependencies other than a compiler that understands the C++11 standard. Random number generation employs the `RDRAND` CPU instruction if supported by the processor, otherwise an implementation of the 64-bit Mersenne Twister \cite matsumoto98a is substituted. Intel Ivy Bridge or later support `RDRAND`. With AMD it is not completely clear. Opteron definitely does not support it. Zen architectures (Ryzen) claim to support it, but I did not have access to one so I cannot personally vouch for it. The RNG choice is made automatically at run time.

# Installation

The simplest way to install everything is to run

	make all
	sudo make install

in the directory with the source code. This will install the executables (_sampleSNPs_ and _sampleLD_) in _/usr/local/bin/_ and the library in the appropriate folders in _/usr/local/_. The included Makefile can be modified to change where things go. The headers _varfiles.hpp_, _populations.hpp_, and _random.hpp_ have to be included in your code as necessary.

# Testing

The _tests/_ directory contains example .bed, .tped, .vcf, and .hmp.txt files to try running the programs on. To keep sizes manageable for distribution, the .bed file has 50,000 loci, while the text files have only 5,000. Make sure your samples do not exceed these values. Uncompress the directory and run, for example,

	./sampleSNPs -i tests/sample -t BED -s 5000

This should sample 5,000 SNPs from the included *sample_ALL.bed* and *sample_ALL.bim* files and save the results into files with the *sample_ALL_s5000* prefix.

Note that _sampleLD_ supports only the .bed format.

Running _sampleSNPs_ and _sampleLD_ without flags will cause these programs to print flag descriptions and exit.

# Timing

Expanding the _timingTrials.tar.gz_ archive will generate a directory with separate software, depending only on random.cpp and .hpp, that performs analyses of execution time using Vitter's Method D and Method S. The README.md file included there explains how to compile and run these analyses.

# Citing this work

The paper describing this work is available, and can be referenced, from [bioRxiv](https://www.biorxiv.org/content/early/2017/11/17/220871) and [arXiv](https://arxiv.org/abs/1711.06325).

