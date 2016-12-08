/* Exponentials */
// Started at 1:45

#include <iostream>     // cout
#include <stdlib.h>     // srand, rand
#include <time.h>       // time
#include <gmd.h>        // LaGenMatDouble
#include <blas3pp.h>
using namespace std;

void seed_random(){
    srand(time(NULL));
}

int ran(int x){
    return rand() % x;
}

void generate_lattice(LaGenMatDouble*){
//    for(int i = 0; i < 4; i++){
//        for(int j = 0; j < 4; j++){ 
//            lattice[i][j] = ran(100);
//        }
//    }
    
//    double* elements = new double[n*n]; 
//    LaGenMatDouble spins = LaGenMatDouble( elements, n, n, row_ordering = true);
//    delete [] elements; // we created a new array on the heap, so lets make sure we delete it
}

void print_lattice(int lattice[4][4]){
//    cout.length(5);
//    for(int i = 0; i < 4; i++){
//        for(int j = 0; j < 4; j++){ 
//            cout << lattice [i][j] << "  ";
//        }
//        cout << endl;
//    }
    
//    cout << "Spins matrix:" << endl << spins << endl;
}

void diagonalise_lattice([4][4]){
//    LaGenMatDouble diag = LaGenMatDouble( n, n );
//    diag = spins.diag();
//    cout << "Diagonal matrix from spins" << endl << diag << endl;
}

void multiply_lattice([4][4]){

int main(){
    seed_random();
    
    LaGenMatDouble::lattice(4, 4);
//    int lattice [4][4];
//    generate_lattice(lattice);
//    print_lattice(lattice);
    
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

/*
Lapac notes

#include <gmd.h> to use LaGenMatDouble Class

LaGenMatDouble LaGenMatDouble::diag ( ) const
    Returns a newly allocated column vector of dimension Nx1 that contains the diagonal of the given matrix. (New in lapackpp-2.4.5)

file dsyev.f  dsyev.f plus dependencies
prec double
for  Computes all eigenvalues, and optionally, eigenvectors of a real
,    symmetric matrix.
gams d4a1


*/

double* elements = new double[n*n]; 

LaGenMatDouble spins = LaGenMatDouble( elements, n, n, row_ordering = true);

delete [] elements; // we created a new array on the heap, so lets make sure we delete it

cout << "Spins matrix:" << endl << spins << endl;

LaGenMatDouble diag = LaGenMatDouble( n, n );

diag = spins.diag();

cout << "Diagonal matrix from spins" << endl << diag << endl;

LaGenMatDouble result = LaGenMatDouble( n, n );

// here we will show what power of 2 would look like

Blas_Mat_Mat_Mult ( const &diag, const &diag, &result, alpha = 1.0, beta = 0.0 );
// this will output, effectively, diag^2 in to the result object, to continue doing the power we might do:

Blas_Mat_Mat_Mult ( const &result, const &diag, &result, alpha = 1.0, beta = 0.0 );
// note, this might error because we are passing &result as a variable and a const, if this is the case the best way around would be to create a *new* result variable, and just iterate etc etc'

// i'm not sure if these things will error, you may have to initialize the diag and result (spins should be aite though)

