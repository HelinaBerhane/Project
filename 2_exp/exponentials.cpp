/* Exponentials */
// Started at 1:45

#include <iostream>     // cout
#include <stdlib.h>     // srand, rand
#include <time.h>       // time
//#include <gmd.h>        // LaGenMatDouble
//#include <blas3pp.h>
#include <cmath>        // temporary for powers
using namespace std;

void seed_random(){
    srand(time(NULL));
}

int ran(int x){
    return rand() % x;
}

int fact(int x){
    int y = x;
    for(int i = x-1; i > 0; i--){
        cout << y << endl;
        y = y * i;
    }
    cout << y << endl;
}

int matrix_exponential(int matrix){
    float partial_result = 0;
    float result = 0;
    // over each iteration, calculate \sum_i{matrix^n/n!}        
    for(float i = 0; i > 10; i++){
        partial_result = matrix;
        result = result + partial_result;
        cout << partial_result << endl;
        //multiply the resulting 
         
        // To Do:
           // find out how long each operation takes
           // or find some way to decide how accurate you want to be
    }
    cout << result << endl;
}

int main(){
    // matrix size
    int n = 5;
    int matrix = 5;
    matrix_exponential(matrix);
    
    // generate the elements
    /*
    seed_random();
    double* elements = new double[n*n];
    for(int i = 0; i = n*n; i++){
        elements[i] = ran(10);
    }
    
    // add the elements to the matrix
    LaGenMatDouble matrix_nn = LaGenMatDouble(elements, n, n, row_ordering = true);
    delete [] elements
    
    // test: print the lattice
    cout << "Matrix:" << endl << matrix_nn << endl;
    
    // diagonalise the matrix
    LaGenMatDouble diag = LaGenMatDouble( n, n );
    diag = spins.diag();
    
    // test: print the diagonalised lattice
    cout << "Diagonal matrix from spins" << endl << diag << endl;
    */
    
    // calculate the exponential

    
//     this is the matrix to exponentiate
//     to do this, we need to:
//        find the eigenvalues and eigenvectors
//        construct matrices with:
//           the eigen values along the diagonal
//           the corresponding eigenvectors in each column
//              U^\dagger is the matrix with the eigenvectors along the columns
//              U is the conjugate matrix
//        in a for loop, calculate powers of this matrix
//           in this, you only have to calculate powers of the exponential matrix
//          
//     then what do you have to do?

}


