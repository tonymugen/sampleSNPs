###############################################################
#
# Makefile for sampling tools
#
###############################################################
CXX = g++
INSTALLDIR = /usr/local
# sample SNPs tool
SMPOUT = sampleSNPs
# sample LD tool
LDOUT = sampleLD
# library
LIBOUT = libsampFiles.a

CXXFLAGS = -O3 -march=native -std=c++11
SMPOBJ = varfiles.o random.o
LDOBJ = varfiles.o random.o populations.o

all : $(SMPOUT) $(LDOUT) $(LIBOUT)
.PHONY : all

install : $(LIBOUT) $(LDOUT) $(SMPOUT)
	-cp -v $(LIBOUT) $(INSTALLDIR)/lib
	-cp -v *.hpp $(INSTALLDIR)/include
	-cp -v $(LDOUT) $(INSTALLDIR)/bin
	-cp -v $(SMPOUT) $(INSTALLDIR)/bin
.PHONY : install

$(SMPOUT) : sampleSNPs.cpp $(SMPOBJ)
	$(CXX) sampleSNPs.cpp $(SMPOBJ) -o $(SMPOUT) $(CXXFLAGS)

$(LDOUT) : sampleLD.cpp $(LDOBJ)
	$(CXX) sampleLD.cpp $(LDOBJ) -o $(LDOUT) $(CXXFLAGS)

$(LIBOUT) : $(LDOBJ)
	ar crv $(LIBOUT) $(LDOBJ)

random.o : random.cpp random.hpp
	$(CXX) -c random.cpp $(CXXFLAGS)

populations.o : populations.cpp populations.hpp
	$(CXX) -c populations.cpp $(CXXFLAGS)

varfiles.o : varfiles.hpp varfiles.cpp populations.o random.o
	$(CXX) -c varfiles.cpp $(CXXFLAGS)

.PHONY : clean
clean :
	-rm *.o $(SMPOUT) $(LDOUT)

