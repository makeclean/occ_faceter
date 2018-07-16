#LIBS = -L/Users/davisa/opt/moab/lib -lMOAB -L/usr/local/Cellar/cgal/4.11_1/lib -lCGAL_Core -lCGAL -L/usr/local/Cellar/gmp/6.1.2_1/lib/ -lgmp -L/usr/local/Cellar/boost/1.66.0/lib/ -lboost_system #-lboost_thread-mt -lboost_exception-mt -lboost_exception
#INC = -I/Users/davisa/opt/moab/include -I/usr/local/Cellar/cgal/4.11_1/include

OCC_LIBS = -L/usr/lib -lTKSTEP -lTKXSBase -lTKernel -lTKMath -lTKBRep -lTKMesh -lTKTopAlgo -lTKBRep -lTKSTEPBase
OCC_INC = -I/usr/include/oce  

MOAB_LIBS += -L/home/adavis/dev/moab/lib -lMOAB
MOAB_INC += -I/home/adavis/dev/moab/include

CGAL_LIBS = -L/home/adavis/opt/cgal/lib/ -lCGAL -lboost_system -lgmp
CGAL_INCLUDES = -I/home/adavis/opt/cgal/include

INC = $(OCC_INC) $(MOAB_INC)

CXX = g++
CXXFLAGS = -g -std=c++11 $(INC)

.PHONY: default all clean

default:  dagmc_faceter dagmc_slicer dagmc_merge
all: default

HEADERS = MBTool.hpp

mbtool.o: MBTool.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c MBTool.cpp -o mbtool.o

dagmc_faceter.o: dagmc_faceter.cc mbtool.o $(HEADERS)
	$(CXX) $(CXXFLAGS) -c dagmc_faceter.cc -o dagmc_faceter.o

dagmc_faceter: mbtool.o dagmc_faceter.o
	$(CXX) $(CXXFLAGS) dagmc_faceter.o mbtool.o $(OCC_LIBS) $(MOAB_LIBS) -o dagmc_faceter

dagmc_topology.o: dagmc_topology.cc
	$(CXX) $(CXXFLAGS) -c dagmc_topology.cc -o dagmc_topology.o

dagmc_merge: dagmc_topology.o
	$(CXX) $(CXXFLAGS) dagmc_topology.o dagmc_merge.cc $(MOAB_LIBS) -o dagmc_merge

dagmc_slicer.o: dagmc_slicer.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) -c dagmc_slicer.cc -o dagmc_slicer.o

dagmc_slicer: dagmc_slicer.o $(HEADERS)
	$(CXX) $(CXXFLAGS) dagmc_slicer.o $(CGAL_LIBS) $(MOAB_LIBS) -o dagmc_slicer

clean:
	-rm -f *.o
	-rm -f dagmc_faceter dagmc_slicer dagmc_merge
