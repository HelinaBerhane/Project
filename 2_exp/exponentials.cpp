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

//void generate_lattice(LaGenMatDouble*){
//    for(int i = 0; i < 4; i++){
//        for(int j = 0; j < 4; j++){ 
//            lattice[i][j] = ran(100);
//        }
//    }
    
//    double* elements = new double[n*n]; 
//    LaGenMatDouble spins = LaGenMatDouble( elements, n, n, row_ordering = true);
//    delete [] elements; // we created a new array on the heap, so lets make sure we delete it
//}

//void print_lattice(int lattice[4][4]){
////    cout.length(5);
//    for(int i = 0; i < 4; i++){
//        for(int j = 0; j < 4; j++){ 
//            cout << lattice [i][j] << "  ";
//        }
//        cout << endl;
//    }
    
////    cout << "Spins matrix:" << endl << spins << endl;
//}

//void diagonalise_lattice([4][4]){
////    LaGenMatDouble diag = LaGenMatDouble( n, n );
////    diag = spins.diag();
////    cout << "Diagonal matrix from spins" << endl << diag << endl;
//}

//void multiply_lattice([4][4]){

int main(){
    
    // generate the lattice
    int n = 4;
    seed_random();
    double* elements = new double[n*n]; 
    for(int i = 0; i = n; i++){
        cout.width(4);
        elements[i] = ran(10);
        cout << elements[i]
    }
    cout << endl;
    
//    LaGenMatDouble::lattice(4, 4);
        
//    print the lattice

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


