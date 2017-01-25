#include <cstdlib>  //rand, srand
#include <stdlib.h> //abs
#include <iostream> //cout
#include <ctime>    //time(NULL)
#include "complex_matrices.h"
#include <blas3pp.h>

using namespace std;

/* Printing */
void print_lattice(int array[10][10]){
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            cout.width(3);
            cout << array[i][j];
        }
        cout << endl;
    }
}

/* lattice generation */
int generate_spins(){ /* Do this today! */
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dist(0, 1);
    return (dist(gen) % 2)*2 - 1;
}
void generate_lattice(int array[10][10]){ /* Do this today! */
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            array[i][j] = generate_spins();
            cout << array[i][j];
        }
        cout << endl;
    }
}

void generate_H(const int matrix_size){ /* Do this today! */
    /* initialise everything */
    COMPLEX elements[matrix_volume];
    int n;
    /* generate the matrix */
    for(int i = 0; i < matrix_size; i++){
        for (int j = 0; j < matrix_size; j++) {
            n = (matrix_size * i) + j;
            cout.width(3);
            cout << abs(i-j);
            if(abs(i-j) = 1 || abs(i-j) = matrix_size - 1){
                elements[n].r = -1;
            }else{
                comp[n].r = 0;
            }
            comp[i+j].i = 0;
        }
        cout << endl;
    }
    cout << endl;
    LaGenMatComplex matrix = LaGenMatComplex(elements, 3, 3, false );
    /* print result */
    print_matrix(matrix, "initial matrix");
}

void test_lattice_generation(){
    int lattice [10][10];
    generate_lattice(lattice);
    print_lattice(lattice);
}

int main(){
    /* test lattice generation */
    test_lattice_generation();

    /* test lattice generation */
    generate_H(5);
}

/*
fix the randomisation using http://stackoverflow.com/questions/19665818/generate-random-numbers-using-c11-random-library
*/
