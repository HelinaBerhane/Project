#these are the ones that work for your home pc - they wont work for the uni pc
#Inc = /usr/lib/lapackpp/include/lapackpp/
#Lib = /usr/lib/lapackpp/lib/

Inc = /home/k1333466/Work/lib/lapackpp/include/lapackpp/
Lib = /home/k1333466/Work/lib/lapackpp/lib/
CXX = g++
# CXXFLAGS = -llapackpp -Wall -g
CXXFLAGS = -llapackpp -O1

all: matrix

matrix: general.o general.h pull
	$(CXX) general.o -I$(Inc) -L$(Lib) $(CXXFLAGS) -o matrix

general.o: general.cpp general.h pull
	$(CXX) -I$(Inc) -L$(Lib) $(CXXFLAGS) -c general.cpp

pull:
	git pull origin master
