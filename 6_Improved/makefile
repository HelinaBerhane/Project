#these are the ones that work for your home pc - they wont work for the uni pc
#Inc = /usr/lib/lapackpp/include/lapackpp/
#Lib = /usr/lib/lapackpp/lib/

Inc = /home/k1333466/Work/lib/lapackpp/include/lapackpp/
Lib = /home/k1333466/Work/lib/lapackpp/lib/
CXX = g++
# CXXFLAGS = -llapackpp -Wall -g
CXXFLAGS = -llapackpp -O1

all: improv

improv: improv.o improv.h pull
	$(CXX) improv.o -I$(Inc) -L$(Lib) $(CXXFLAGS) -o improv

improv.o: improv.cpp improv.h pull
	$(CXX) -I$(Inc) -L$(Lib) $(CXXFLAGS) -c improv.cpp

pull:
	git pull origin master
