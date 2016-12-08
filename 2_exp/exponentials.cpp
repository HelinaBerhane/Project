/* Exponentials */
// Started at 21:49

#include <iostream>     // cout
#include <stdlib.h>     // srand, rand
#include <time.h>       // time
using namespace std;

void seed_random(){
    srand(time(NULL));
}

int ran(int x){
    return rand() % x;
}

void generate_lattice(int lattice[4][4]){
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){ 
            lattice[i][j] = ran(100);
        }
    }
}

void print_lattice(int lattice[4][4]){
//    cout.length(5);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){ 
            cout << lattice [i][j] << "  ";
        }
        cout << endl;
    }
}

int main(){
    seed_random();
    int lattice [4][4];
    generate_lattice(lattice);
    print_lattice(lattice);
    // this is the matrix to exponentiate
    // to do this, we need to:
       // find the eigenvalues and eigenvectors
       // construct matrices with:
          // the eigen values along the diagonal
          // the corresponding eigenvectors in each column
             // U^\dagger is the matrix with the eigenvectors along the columns
             // U is the conjugate matrix
       // in a for loop, calculate powers of this matrix
          // in this, you only have to calculate powers of the exponential matrix
          
    // then what do you have to do?

}
