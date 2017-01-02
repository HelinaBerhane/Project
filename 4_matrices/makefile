Inc = /usr/lib/lapackpp/include/lapackpp/
Lib = /usr/lib/lapackpp/lib/
CXX = g++
CXXFLAGS = -llapackpp -Wall -g

all: matrix
	
matrix: complex_matrices.o
	$(CXX) complex_matrices.o -I$(Inc) -L$(Lib) $(CXXFLAGS) -o matrix
	
complex_matrices.o: complex_matrices.cpp
	$(CXX) -I$(Inc) -L$(Lib) $(CXXFLAGS) -c complex_matrices.cpp