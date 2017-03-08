#these are the ones that work for your home pc - they wont work for the uni pc
#Inc = /usr/lib/lapackpp/include/lapackpp/
#Lib = /usr/lib/lapackpp/lib/

Inc = /home/k1333466/Work/lib/lapackpp/include/lapackpp/
Lib = /home/k1333466/Work/lib/lapackpp/lib/
CXX = g++
CXXFLAGS = -llapackpp -Wall -g

all: matrix

matrix: double.o double.h pull
	$(CXX) double.o -I$(Inc) -L$(Lib) $(CXXFLAGS) -o double

double.o: double.cpp double.h pull
	$(CXX) -I$(Inc) -L$(Lib) $(CXXFLAGS) -c double.cpp

pull:
	git pull origin master
