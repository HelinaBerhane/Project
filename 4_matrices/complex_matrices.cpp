#define LA_COMPLEX_SUPPORT

#include <iostream> //cout
#include <cstdlib>	//rand, sran
#include <gmc.h> 	//LaGenMatComplex
#include <string>
#include <laslv.h>  //LUFactorizeIP, LaLUInverseIP, etc.
using namespace std;

/* Testing */
void print_matrix(LaGenMatComplex* matrix){					                    //working
	cout << *matrix << endl;
}
void print_matrix(LaGenMatComplex* matrix, const string name){					//working
	cout << name << ":" << endl << *matrix << endl;
}

/* Matrix generation */
void generate_elements(COMPLEX* elements, int m, int x){	                    //working
	srand(time(NULL));				//seed
    for(int i = 0; i < m*m; i++){
        elements[i].r = rand() % x;	//1 to x
        elements[i].i = rand() % x;
	}
}

/* Matrix manipulation */
void matrix_eigenvstuff (){
    //LaEigSolve: http://lapackpp.sourceforge.net/html/laslv_8h.html#086357d17e9cdcaec69ab7db76998769
//    LaEigSolve();
}                                                   //empty
void matrix_powers(){
	// Blas_Mat_Mat_Mult: http://lapackpp.sourceforge.net/html/blas3pp_8h.html#3b0c7a2c6c951b0d5d9d443981293b0c
}                                                        //empty
void matrix_inverse(LaGenMatComplex* matrix, int len){
    // LaLUInverseIP: http://lapackpp.sourceforge.net/html/laslv_8h.html#a042c82c5b818f54e7f000d068f14189
    LaVectorLongInt PIV = LaVectorLongInt(len);
    LUFactorizeIP(*matrix, PIV);
    LaLUInverseIP(*matrix, PIV);
}                       //untested

/* Main Program */
int main(){
	
    /* initialise stuff */
    
	/* generate the matrix */           //working
	int m = 2, x = 10;                                                     //m = matrix size, x = max
	COMPLEX* elements = new COMPLEX[m*m];                                   //initialise elements
	generate_elements(elements, m, x);                                      //generate random elements
	LaGenMatComplex* initialMatrix = new LaGenMatComplex( elements, m, m, false ); //copy the elements to a matrix
    
    /* save and print the matrix */     //working
    print_matrix(initialMatrix, "initial matrix");   
    
	/* calculate the eigenvalue matrix */
//    LaVectorComplex* eigenvalues ;
//    LaGenMatDouble* eigenvectors;
    
//    matrix_eignevalues(matrix);
    
    /* save and print the eigenvalue matrix */
//    LaGenMatComplex* eigenvalueMatrix;
//    eigenvalueMatrix = matrix->copy();
//    print_matrix(eigenvalueMatrix, "eigenvalue matrix");   
//    matrix = initialMatrix->copy();                      //reset matrix
    
	/* calculate the eigenvector matrix */
    
    /* save and print the eigenvector matrix */
//    LaGenMatComplex* eigenvectorMatrix;
//    eigenvectorMatrix = matrix->copy();
//    print_matrix(eigenvectorMatrix, "eigenvector matrix"); 
//    matrix = initialMatrix->copy();                      //reset matrix
    
	/* calculate the transpose eigenvector matrix */
    
    /* save and print the transpose eigenvector matrix */
//    LaGenMatComplex* tEigenvectorMatrix;
//    tEigenvectorMatrix = matrix->copy();
//    print_matrix(tEigenvectorMatrix, "transpose eigenvector matrix");
//    matrix = initialMatrix->copy();                      //reset matrix
    
	/* calculate the power matrix */
    
    /* save and print the power matrix */
    
	/* calculate the exponential matrix */
    
    /* save and print the exponential matrix */
    
    /* invert the matrix */
    
    /* save and print the inverse matrix */
    LaGenMatComplex inverseMatrix;
    inverseMatrix = initialMatrix->copy();
    matrix_inverse(&inverseMatrix, m);
    print_matrix(&inverseMatrix, "inverse matrix");          
    
	/* delete things on the heap: */
    delete [] elements;
//    delete eigenvalues;
//    delete eigenvectors;
	delete initialMatrix;
	delete initialMatrix;
}
