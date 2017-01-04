#define LA_COMPLEX_SUPPORT

#include <iostream> //cout
#include <cstdlib>	//rand, sran
#include <gmc.h> 	//LaGenMatComplex
#include <string>
#include <laslv.h>  //LUFactorizeIP, LaLUInverseIP, etc.
#include <blas3pp.h>
using namespace std;

/* Testing */
void print_matrix(const LaGenMatComplex& matrix){
	cout << matrix << endl;
}                                      //working
void print_matrix(const LaGenMatComplex& matrix, const string name){
	cout << name << ":" << endl << matrix << endl;
}			        //working
void print_array(const COMPLEX array[], int len){
    for(int i = 0; i < len; i++){
        cout << array[i] << endl;
    }
}                                      //working
void print_array(const COMPLEX array[], int len, const string name){
	cout << name << ":" << endl;
    for(int i = 0; i < len; i++){
        cout << array[i] << endl;
    }
}                   //working

/* Matrix generation */
void generate_elements(COMPLEX elements[], const int len, const int x){
	srand(time(NULL));				//seed
    for(int i = 0; i < len*len; i++){
        elements[i].r = rand() % x;	//1 to x
        elements[i].i = rand() % x;
	}
}                //working

/* Matrix conversion */
void vec_to_array(const LaVectorComplex& vector, const int len, COMPLEX* array){
    for(int i = 0; i < len; i++){
        array[i] = vector(i);
    }
}       //working
void array_to_diag(COMPLEX* array, const int len, LaGenMatComplex& diag){
    diag = 0;
    for(int i = 0; i < len; i++){
        diag(i, i) = array[i];        
    }
}                //checked
void vec_to_diag(const LaVectorComplex& vector, const int len, LaGenMatComplex& diag){
    COMPLEX array[len];
    vec_to_array(vector, len, array);
    array_to_diag(array, len, diag);
}   //checked

/* Scalar manipulation */
void complex_power(const COMPLEX& number, const int power, COMPLEX& answer){
     
}           //empty

/* Matrix manipulation */
void matrix_eigenvstuff(const LaGenMatComplex& matrix, LaVectorComplex& eigenvalues, LaGenMatComplex& eigenvectors){
    //LaEigSolve: http://lapackpp.sourceforge.net/html/laslv_8h.html#086357d17e9cdcaec69ab7db76998769
    LaEigSolve(matrix, eigenvalues, eigenvectors);
} //empty
void matrix_powers(){
	// Blas_Mat_Mat_Mult: http://lapackpp.sourceforge.net/html/blas3pp_8h.html#3b0c7a2c6c951b0d5d9d443981293b0c
}                                                                  //empty
void matrix_inverse(LaGenMatComplex& matrix, int len){
    // LaLUInverseIP: http://lapackpp.sourceforge.net/html/laslv_8h.html#a042c82c5b818f54e7f000d068f14189
    LaVectorLongInt PIV = LaVectorLongInt(len);
    LUFactorizeIP(matrix, PIV);
    LaLUInverseIP(matrix, PIV);
}                                 //working

/* Main Program */
int main(){
    //m = matrix size, x = max                                                            //working
	int m = 3, x = 10;                                                              
    
	/* generate the matrix */                                                             //working
    COMPLEX elements[m*m];
    generate_elements(elements, m, x);
	LaGenMatComplex initialMatrix = LaGenMatComplex(elements, m, m, false );
    print_matrix(initialMatrix, "initial matrix");
    
    /* initialise eigenstuff */                                                           //working
    LaVectorComplex eigenvalueVec = LaVectorComplex(m);
    LaGenMatComplex eigenvalues = LaGenMatComplex::zeros(m, m);
    LaGenMatComplex eigenvectors = LaGenMatComplex::zeros(m, m);
    
	/* calculate the eigenvalues and eigenvectors */                                      //working
    matrix_eigenvstuff(initialMatrix, eigenvalueVec, eigenvectors);
    
    /* check eigenvalues and eigenvectors*/                                               //working
    print_matrix(eigenvalueVec, "eigenvalue vector");
    vec_to_diag(eigenvalueVec, m, eigenvalues);
    print_matrix(eigenvalues, "eigenvalue matrix");
    print_matrix(eigenvectors, "eigenvector matrix");
    
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
    inverseMatrix = initialMatrix.copy();
    matrix_inverse(inverseMatrix, m);
    print_matrix(inverseMatrix, "inverse matrix");          
    
//    delete eigenvalues;
//    delete eigenvectors;
}
