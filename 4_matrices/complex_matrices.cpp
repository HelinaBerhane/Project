#define LA_COMPLEX_SUPPORT

#include <iostream> //cout
#include <cstdlib>	//rand, sran
#include <gmc.h> 	//LaGenMatComplex
#include <string>
#include <laslv.h>  //LUFactorizeIP, LaLUInverseIP, etc.
#include <blas3pp.h>
using namespace std;

/* Total [13/20] */

/* Testing [6/6] */
void print_scalar(const COMPLEX scalar){
    cout << scalar << endl;
}                                                 //working
void print_scalar(const COMPLEX scalar, const string name){
    cout << name << ":" << scalar << endl;
}                              //working
void print_array(const COMPLEX array[], int len){
    for(int i = 0; i < len; i++){
        cout << array[i] << endl;
    }
}                                        //working
void print_array(const COMPLEX array[], int len, const string name){
	cout << name << ":" << endl;
    for(int i = 0; i < len; i++){
        cout << array[i] << endl;
    }
}                     //working
void print_matrix(const LaGenMatComplex& matrix){
	cout << matrix << endl;
}                                        //working
void print_matrix(const LaGenMatComplex& matrix, const string name){
	cout << name << ":" << endl << matrix << endl;
}			          //working

/* Number generation [2/2] */
void generate_scalar(COMPLEX& scalar, const int x){
    scalar.r = rand() % x;	//1 to x
    scalar.i = rand() % x;
}                                      //working
void generate_array(COMPLEX array[], const int len, const int x){
	srand(time(NULL));				//seed
    for(int i = 0; i < len; i++){
        array[i].r = rand() % x;	//1 to x
        array[i].i = rand() % x;
	}
}                        //working

/* Matrix conversion [3/3] */
void vec_to_array(const LaVectorComplex& vector, const int len, COMPLEX array[ ]){
    for(int i = 0; i < len; i++){
        array[i] = vector(i);
    }
}       //working
void array_to_diag(COMPLEX array[], const int len, LaGenMatComplex& diag){
    diag = 0;
    for(int i = 0; i < len; i++){
        diag(i, i) = array[i];        
    }
}               //working
void vec_to_diag(const LaVectorComplex& vector, const int len, LaGenMatComplex& diag){
    COMPLEX array[len];
    vec_to_array(vector, len, array);
    array_to_diag(array, len, diag);
}   //working

/* Scalar manipulation [4/6] */
int factorial(int x){
	if(x <= 1){
        return 1;
	}else{
        return x * factorial(x - 1);
	}
}                                                                    //working
void complex_addition(const COMPLEX& numberA, const COMPLEX& numberB , COMPLEX& result){
    result.r = numberA.r + numberB.r;
    result.i = numberA.i + numberB.i;
} //working
void complex_powers(const COMPLEX& number, const int power, COMPLEX& result){
    la::complex<double> laNumber = la::complex<double>(number);
    la::complex<double> init = la::complex<double>(number);
    for(int i = 1; i < power; i++){
        laNumber *= init;
    }
    result = laNumber.toCOMPLEX();
}            //working
void complex_multiplication(){
    
}                                                           //empty
void complex_division(const COMPLEX& number, const int real, COMPLEX& result){
    result.r = number.r / real;
    result.i = number.i / real;
}           //working
void complex_exponential(const COMPLEX& number, const int iter, COMPLEX& result){
    /*
    for(int i = 0; i < iter; i++){
        powerStep
    }
    */
}        //empty

/* array manipulation [0/1] */
void array_power(COMPLEX array[], const int len, const int power){
    /*
    for(int i = 0; i < len; i++){
        array[i] = complex_power(array[i], power, result);
    }
    */
}                       //empty

/* Matrix manipulation [2/4] */
void matrix_eigenvstuff(const LaGenMatComplex& matrix, LaVectorComplex& eigenvalues, LaGenMatComplex& eigenvectors){
    //LaEigSolve: http://lapackpp.sourceforge.net/html/laslv_8h.html#086357d17e9cdcaec69ab7db76998769
    LaEigSolve(matrix, eigenvalues, eigenvectors);
} //working
void diagonal_power(){
    
}                                                                 //empty
void matrix_inverse(LaGenMatComplex& matrix, int len){
    // LaLUInverseIP: http://lapackpp.sourceforge.net/html/laslv_8h.html#a042c82c5b818f54e7f000d068f14189
    LaVectorLongInt PIV = LaVectorLongInt(len);
    LUFactorizeIP(matrix, PIV);
    LaLUInverseIP(matrix, PIV);
}                                 //working
void matrix_exponential(const LaGenMatComplex& eigenvectors, const LaGenMatComplex& eigenvalues){
    //LaGenMatComplex step = LaGenMatComplex
} //empty

/* Main Program */
int main(){
    //m = matrix size, x = max
	int m = 3, x = 10;                                                              
    
	/* generate the matrix */
    COMPLEX elements[m*m];
    generate_array(elements, m*m, x);
	LaGenMatComplex initialMatrix = LaGenMatComplex(elements, m, m, false );
    print_matrix(initialMatrix, "initial matrix");
    
    /* initialise eigenstuff */
    LaVectorComplex eigenvalueVec = LaVectorComplex(m);
    LaGenMatComplex eigenvalues = LaGenMatComplex::zeros(m, m);
    LaGenMatComplex eigenvectors = LaGenMatComplex::zeros(m, m);
    
	/* calculate the eigenvalues and eigenvectors */
    matrix_eigenvstuff(initialMatrix, eigenvalueVec, eigenvectors);
    
    /* check eigenvalues and eigenvectors*/
    print_matrix(eigenvalueVec, "eigenvalue vector");
    vec_to_diag(eigenvalueVec, m, eigenvalues);
    print_matrix(eigenvalues, "eigenvalue matrix");
    print_matrix(eigenvectors, "eigenvector matrix");
    
    /* test scalar manipulation */
    COMPLEX testA;                             //initialisation
    COMPLEX testB;
    COMPLEX testAdd;
    COMPLEX testPow;
    COMPLEX testDiv;
    for(int i = 1; i < 5; i++){
        cout << "factorial(" << i << "): " << factorial(i) << endl;
    }            //factorials
    cout << endl;
    generate_scalar(testA, x);                //generation
    cout << "testA" << testA << endl;
    generate_scalar(testB, x);
    cout << "testB" << testB << endl;
    cout << endl;
    complex_addition(testA, testB, testAdd);  //addition
    cout << "testAdd" << testAdd << endl;
    cout << endl;
    complex_powers(testA, 3, testPow);        //powers
    cout << "testA cubed" << testPow << endl;
    cout << endl;
    for(int i = 1; i < 4; i++){
        complex_division(testA, factorial(i), testDiv);
        cout << "testA = " << testA << " | ";
        cout << "factorial(" << i << ") = " << factorial(i) << " | ";
        cout << "testA/factorial(" << i << "): " << testDiv << endl;
    }            //division
    cout << endl;
    
    //cout << "division" << testA/4 << endl;
    
    
    /* invert the matrix */
    LaGenMatComplex inverseMatrix;
    inverseMatrix = initialMatrix.copy();
    matrix_inverse(inverseMatrix, m);
    print_matrix(inverseMatrix, "inverse matrix");
}
