#these are the ones that work for your home pc - they wont work for the uni pc
#Inc = /usr/lib/lapackpp/include/lapackpp/
#Lib = /usr/lib/lapackpp/lib/

Inc = /home/k1333466/Work/lib/lapackpp/include/lapackpp/
Lib = /home/k1333466/Work/lib/lapackpp/lib/
CXX = g++
CXXFLAGS = -llapackpp -Wall -g

all: matrix

matrix: complex_matrices.o complex_matrices.h
	$(CXX) complex_matrices.o -I$(Inc) -L$(Lib) $(CXXFLAGS) -o matrix

complex_matrices.o: complex_matrices.cpp complex_matrices.h
	$(CXX) -I$(Inc) -L$(Lib) $(CXXFLAGS) -c complex_matrices.cpp
