ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

CXX           = g++
CXXFLAGS      = -g -Wall -ansi -fPIC -I$(ROOTSYS)/include
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared

NGLIBB         = $(ROOTGLIBS)
NGLIBB        += -lMinuit -lTreePlayer -lMLP
#/home/andrea/PhDOxford/OwnRootClasses/CutsArray_cpp.so /home/andrea/PhDOxford/OwnRootClasses/OptCut2_cpp.so

GLIBB          = $(filter-out -lNew, $(NGLIBB))

CXXFLAGS      += $(ROOTCFLAGS)
LIBS           = $(ROOTLIBS)
LIBS          += $(ROOTSYS)/lib/*.sl -lXpm -lX11 -lm -ldld
.SUFFIXES: .cc,.C

Exec_tag:  MakeTrends.C
# -----------------------------------------------------------------------------
	$(CXX) $(CXXFLAGS) -c $<
	$(LD) $(LDFLAGS) -o MakeTrends MakeTrends.o $(GLIBB)

# ================================================================================
clean:
	rm -f *.o MakeTrends
# -----------------------------------------------------------------------------