//#include <stdlib.h>
//#include <stdio.h>
#include <iostream> //cout
#include <random>
//#include <cblas.h>

using namespace std;

/*
struct _dcomplex { double re, im; };
typedef struct _dcomplex dcomplex;
*/

float rnd(float x){
    mt19937_64 y;
    uniform_real_distribution<double> doubleDist(0, 1);
    return x * doubleDist(y);
}

int main(){    
    // Source: http://stackoverflow.com/questions/23324480/matrix-operations-in-c-using-blas-lapack-or-some-other-alternative
    cout << rnd(100) << endl;
}

