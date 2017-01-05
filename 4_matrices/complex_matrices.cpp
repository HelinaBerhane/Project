#define LA_COMPLEX_SUPPORT

#include <iostream> //cout
#include <cstdlib>	//rand, sran
#include <gmc.h> 	//LaGenMatComplex
#include <string>
#include <laslv.h>  //LUFactorizeIP, LaLUInverseIP, etc.
#include <blas3pp.h>
using namespace std;

/* function primitives */
void matrix_eigenvstuff(const LaGenMatComplex& matrix, LaVectorComplex& eigenvalues, LaGenMatComplex& eigenvectors);
void matrix_inverse(LaGenMatComplex& matrix, int len);
void matrix_exponential(const LaGenMatComplex& eigenvectors, const LaGenMatComplex& eigenvalues);
void vec_to_diag(const LaVectorComplex& vector, const int len, LaGenMatComplex& diag);
void print_scalar(const COMPLEX scalar);
void print_scalar(const COMPLEX scalar, const string name);
void print_array(const COMPLEX array[], int len);
void print_array(const COMPLEX array[], int len, const string name);
void print_matrix(const LaGenMatComplex& matrix);
void print_matrix(const LaGenMatComplex& matrix, const string name);
void test_eigenvalues(const LaGenMatComplex& initialMatrix, const int size);
void test_scalar_manipulation(const int max_rand);
void test_inverse(const LaGenMatComplex& initialMatrix, const int size);
void generate_scalar(COMPLEX& A, const int x);
void generate_scalar(int A, const int x);
void generate_array(COMPLEX array[], const int len, const int x);
void vec_to_array(const LaVectorComplex& vector, const int len, COMPLEX array[ ]);
void array_to_diag(COMPLEX array[], const int len, LaGenMatComplex& diag);
void vec_to_diag(const LaVectorComplex& vector, const int len, LaGenMatComplex& diag);
int factorial(int x);
void scalar_addition(const COMPLEX& A, const COMPLEX& B , COMPLEX& result);
void scalar_multiplication(const COMPLEX& A, const int B, COMPLEX& result);
void scalar_multiplication(const COMPLEX& A, const COMPLEX& B, COMPLEX& result);
void scalar_division(const COMPLEX& A, const int B, COMPLEX& result);
void scalar_division(const COMPLEX& A, const COMPLEX& B, COMPLEX& result);
void scalar_powers(const COMPLEX& number, const int power, COMPLEX& result);
void scalar_exponential(const COMPLEX& number, const int iter, COMPLEX& result);
void array_power(COMPLEX array[], const int len, const int power);
void diagonal_matrix_powers();
void matrix_eigenvstuff(const LaGenMatComplex& matrix, LaVectorComplex& eigenvalues, LaGenMatComplex& eigenvectors);
void matrix_inverse(LaGenMatComplex& matrix, int len);
void matrix_exponential(const LaGenMatComplex& eigenvectors, const LaGenMatComplex& eigenvalues);

/* Total [13/20] */

/* Printing [6/9] */
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
void generate_scalar(COMPLEX& A, const int x){
    A.r = rand() % x;	//1 to x
    A.i = rand() % x;
}                                           //working
void generate_scalar(int A, const int x){
    A = rand() % x;	//1 to x
}                                                //working
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

/* Scalar manipulation [5/8] */
int factorial(int x){
	if(x <= 1){
        return 1;
	}else{
        return x * factorial(x - 1);
	}
}                                                                    //working
void scalar_addition(const COMPLEX& A, const COMPLEX& B , COMPLEX& result){
    result.r = A.r + B.r;
    result.i = A.i + B.i;
}              //working
void scalar_multiplication(const COMPLEX& A, const int B, COMPLEX& result){
    result.r = A.r * B;
    result.i = A.i * B;
}              //to test
void scalar_multiplication(const COMPLEX& A, const COMPLEX& B, COMPLEX& result){
    la::complex<double> laA = la::complex<double>(A); //convert to la::complex<double>
    la::complex<double> laB = la::complex<double>(B);
    la::complex<double> laResult = la::complex<double>(result);
    laResult = laA * laB;
    result = laResult.toCOMPLEX();
}         //to test
void scalar_division(const COMPLEX& A, const int B, COMPLEX& result){
    result.r = A.r / B;
    result.i = A.i / B;
}                    //to test
void scalar_division(const COMPLEX& A, const COMPLEX& B, COMPLEX& result){
    la::complex<double> laA = la::complex<double>(A); //convert to la::complex<double>
    la::complex<double> laB = la::complex<double>(B);
    la::complex<double> laResult = la::complex<double>(result);
    laResult = laA / laB;
    result = laResult.toCOMPLEX();
}               //to test
void scalar_powers(const COMPLEX& number, const int power, COMPLEX& result){
    la::complex<double> laResult = la::complex<double>(number);
    la::complex<double> laNumber = la::complex<double>(number);
    for(int i = 1; i < power; i++){
        laResult *= laNumber;
    }
    result = laResult.toCOMPLEX();
}             //to test
void scalar_exponential(const COMPLEX& number, const int iter, COMPLEX& result){
    /*
    for(int i = 0; i < iter; i++){
        powerStep
    }
    */
}         //empty

/* array manipulation [0/1] */
void array_powers(COMPLEX array[], const int len, const int power){
    /*
    for(int i = 0; i < len; i++){
        array[i] = complex_power(array[i], power, result);
    }
    */
}                       //empty

/* Matrix manipulation [2/4] */
void diagonal_matrix_powers(){
    
}                                                           //empty
void matrix_eigenvstuff(const LaGenMatComplex& matrix, LaVectorComplex& eigenvalues, LaGenMatComplex& eigenvectors){
    //LaEigSolve: http://lapackpp.sourceforge.net/html/laslv_8h.html#086357d17e9cdcaec69ab7db76998769
    LaEigSolve(matrix, eigenvalues, eigenvectors);
} //working
void matrix_inverse(LaGenMatComplex& matrix, int len){
    // LaLUInverseIP: http://lapackpp.sourceforge.net/html/laslv_8h.html#a042c82c5b818f54e7f000d068f14189
    LaVectorLongInt PIV = LaVectorLongInt(len);
    LUFactorizeIP(matrix, PIV);
    LaLUInverseIP(matrix, PIV);
}                                   //working
void matrix_exponential(const LaGenMatComplex& eigenvectors, const LaGenMatComplex& eigenvalues){
    //LaGenMatComplex step = LaGenMatComplex
} //empty

/* Testing [3/3] */
void test_eigenvalues(const LaGenMatComplex& initialMatrix, const int size){
    LaVectorComplex eigenvalueVec = LaVectorComplex(size);             //initialise eigenstuff
    LaGenMatComplex eigenvalues = LaGenMatComplex::zeros(size, size);
    LaGenMatComplex eigenvectors = LaGenMatComplex::zeros(size, size);
    matrix_eigenvstuff(initialMatrix, eigenvalueVec, eigenvectors); //calculate eigenstuff
    print_matrix(eigenvalueVec, "eigenvalue vector");               //print eigenstuff
    vec_to_diag(eigenvalueVec, size, eigenvalues);
    print_matrix(eigenvalues, "eigenvalue matrix");
    print_matrix(eigenvectors, "eigenvector matrix");
}             //to test
void test_scalar_manipulation(const int max_rand){
    COMPLEX compA;                                                              //initialisation
    generate_scalar(compA, max_rand);
    COMPLEX compB;
    generate_scalar(compB, max_rand);
    int realB;
    generate_scalar(realB, max_rand);
    COMPLEX result;
    
    for(int i = 1; i < 5; i++){                                                 //factorials
        cout << "factorial(" << i << "): " << factorial(i) << endl;
    }
    
    scalar_addition(compA, compB, result);                                      //addition/subtraction
    cout << "scalar addition: " << result << endl << endl;
    
    scalar_multiplication(compA, realB, result);                                //multiplication
    cout << "scalar multiplication by scalar: " << result << endl << endl;
    scalar_multiplication(compA, compB, result);
    cout << "scalar multiplication by complex: " << result << endl << endl;
    
    scalar_division(compA, realB, result);                                      //division
    cout << "scalar division by scalar: " << result << endl << endl;
    scalar_division(compA, compB, testResult);
    cout << "scalar division by complex: " << result << endl << endl;
    
    scalar_powers(compA, compB, result);
    cout << "scalar powers: " << result << endl << endl;
}                                       //to test
void test_inverse(const LaGenMatComplex& initialMatrix, const int size){
    LaGenMatComplex inverseMatrix;
    inverseMatrix = initialMatrix.copy();
    matrix_inverse(inverseMatrix, size);
    print_matrix(inverseMatrix, "inverse matrix");
}                                                                     //to test

/* Main Program */
int main(){
    //m = matrix size, x = max
	int matrix_size = 3, max_rand = 10;                                                              
    int matrix_volume = matrix_size * matrix_size
	/* generate the matrix */
    COMPLEX elements[matrix_volume];
    generate_array(elements, matrix_volume, x);
	LaGenMatComplex initialMatrix = LaGenMatComplex(elements, m, m, false );
    print_matrix(initialMatrix, "initial matrix");
    
    /* test eigenvalues */
    test_eigenvalues(initialMatrix, matrix_size);
    
    /* test scalar manipulation */
    test_scalar_manipulation(x);

    /* test inversion */
    test_inverse(initialMatrix, matrix_size);
}
