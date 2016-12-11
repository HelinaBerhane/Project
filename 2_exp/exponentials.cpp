/* Exponentials */
// Started at 1:45

#include <iostream>     // cout
#include <stdlib.h>     // srand, rand
#include <time.h>       // time
//#include <gmd.h>        // LaGenMatDouble
//#include <blas3pp.h>
using namespace std;

void random_seed(){
    srand(time(NULL));
}

int ran(int x){
    return rand() % x;
}

int fact(int x){
	if(x<=1){
        return 1;
	}else{
        return x*fact(x-1);
	}
}

long double powers(long long int x, long long int y){
    if(y<=1){
        return x;
    }else{
        y--;
        return x * powers(x, y);
    }
}

long double expo(long long int matrix, long long int iterations, long long int n){
    if(iterations == 1){
        iterations--;
        return powers(matrix, n) / fact(n);
    }else{
        cout << powers(matrix, n) << ", " << fact(n) << endl;
        return powers(matrix, n) / fact(n) + expo(matrix, iterations-1, n+1);
    }
    /*
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
    
    int y = 3;
    int powers_test = powers();
    return powers_test;
    */
}

int main(){
    // matrix size
    int iterations = 15;
    int matrix = 3;
    int n = 1;
    cout << expo(matrix, iterations, n) << endl;
//    matrix_exponential(matrix);
    
    // generate the elements
    /*
    random_seed();
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


