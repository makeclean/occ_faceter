#LIBS = -L/Users/davisa/opt/moab/lib -lMOAB -L/usr/local/Cellar/cgal/4.11_1/lib -lCGAL_Core -lCGAL -L/usr/local/Cellar/gmp/6.1.2_1/lib/ -lgmp -L/usr/local/Cellar/boost/1.66.0/lib/ -lboost_system #-lboost_thread-mt -lboost_exception-mt -lboost_exception
#INC = -I/Users/davisa/opt/moab/include -I/usr/local/Cellar/cgal/4.11_1/include

LIBS = -L/usr/lib -lTKSTEP -lTKXSBase -lTKernel -lTKMath -lTKBRep -lTKMesh -lTKTopAlgo -lTKBRep
INC = -I/usr/include/oce  

LIBS += -L/home/adavis/dev/moab/lib -lMOAB
INC += -I/home/adavis/dev/moab/include

CXX = g++
CXXFLAGS = -g -std=c++11 $(INC)

.PHONY: default all clean

default:  program
all: default

HEADERS = MBTool.hpp

mbtool.o: MBTool.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c MBTool.cpp -o mbtool.o

program.o: main.cc mbtool.o $(HEADERS)
	$(CXX) $(CXXFLAGS) -c main.cc -o main.o

program: mbtool.o main.o
	$(CXX) $(CXXFLAGS) main.o mbtool.o $(LIBS) -o main

clean:
	-rm -f main.o mbtool.o
	-rm -f main
