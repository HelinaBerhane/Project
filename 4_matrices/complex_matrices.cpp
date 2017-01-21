#include <iostream> //cout
#include <cstdlib>	//rand, sran
#include <string>
#include "complex_matrices.h"
#include <gmc.h> 	//LaGenMatComplex
#include <laslv.h>  //LUFactorizeIP, LaLUInverseIP, etc.
#include <blas3pp.h>
//#include <random>   //random_device, mt19937

using namespace std;

/* Total [13/20] */

int ran(int a, int b){
//    random_device rd;
//    mt19937 gen(rd());
//    uniform_int_distribution<> dist(a, b);
//    return dist(gen);
    return rand() % b;
}

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
    A.r = ran(1, x);	//1 to x
    A.i = ran(1, x);
}               //working
void generate_scalar(int number, const int x){
    number = ran(1, x);	//1 to x
}                                                //working
void generate_array(COMPLEX array[], const int len, const int x){
    for(int i = 0; i < len; i++){
        array[i].r = ran(1, x);	//1 to x
        array[i].i = ran(1, x);
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
void scalar_addition(COMPLEX& result, const COMPLEX addition){
    result.r += addition.r;
    result.i += addition.i;
}
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
void scalar_product(COMPLEX& total, const COMPLEX& number){
    COMPLEX part;
    part.r = (total.r * number.r) - (total.i * number.i);
    part.i = (total.r * number.i) + (total.i * number.r);
    total = part;
}
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
void scalar_exponential_main(const COMPLEX& number, const int iterations, COMPLEX& result){
    COMPLEX division, total_division;
    result.r = 1;
    result.i = 0 ;
    for(int step = 1; step <= iterations; step++){   //sum (from 1 to n)
        total_division.r = 1;
        total_division.i = 0;
        for(int i = 1; i <= step; i++){        //    ( num^n / n!)
            scalar_division(number, i, division);
            scalar_product(total_division, division);
        }
        scalar_addition(result, total_division);
    }
}
void scalar_exponential(const COMPLEX& number, const int iter, COMPLEX& result){
    //COMPLEX power;
    //COMPLEX division;
    //result.r = 0;
    //result.i = 0;
    //for(int i = 0; i < iter; i++){
    //    scalar_powers(number, i, power);
    //    scalar_division(power, factorial(i), division);
    //    scalar_addition(result, division, result);
    //}
}         //empty
//COMPLEX rec_scalar_exp_step(const COMPLEX& number, const int step){
//    COMPLEX result, division, multiplication;
//	if(step <= 1){
//        result.r = 1;
//        return result;
//	}else{
//        scalar_division(number,step,division);
//        scalar_multiplication(division, rec_scalar_exp_step(step-1),  multiplication);
//        return multiplication;
//	}
//}
void recursive_scalar_exponential(const COMPLEX& number, const int iter, COMPLEX& result){
    //COMPLEX power;
    //COMPLEX division;
    //result.r = 0;
    //result.i = 0;
    //for(int i = 0; i < iter; i++){
    //    scalar_powers(number, i, power);
    //    scalar_division(power, factorial(i), division);
    //    scalar_addition(result, division, result);
    //}
} //empty

/* array manipulation [0/1] */
void array_powers(COMPLEX array[], const int len, const int power){/**/
    /*
    for(int i = 0; i < len; i++){
        array[i] = complex_power(array[i], power, result);
    }
    */
}                       //empty

/* Matrix manipulation [2/4] */
//void diagonal_matrix_powers(){      //empty
  //...
//}
void matrix_eigenvstuff(const LaGenMatComplex& matrix, LaVectorComplex& eigenvalues, LaGenMatComplex& eigenvectors){ //working
    //LaEigSolve: http://lapackpp.sourceforge.net/html/laslv_8h.html#086357d17e9cdcaec69ab7db76998769
    LaEigSolve(matrix, eigenvalues, eigenvectors);
}
void matrix_inverse(LaGenMatComplex& matrix, int len){                                                               //working
    // LaLUInverseIP: http://lapackpp.sourceforge.net/html/laslv_8h.html#a042c82c5b818f54e7f000d068f14189
    LaVectorLongInt PIV = LaVectorLongInt(len);
    LUFactorizeIP(matrix, PIV);
    LaLUInverseIP(matrix, PIV);
}
void matrix_exponential(const LaGenMatComplex& matrix, const int size, LaGenMatComplex& result){
    //initialise eigenstuff
    LaVectorComplex eigenvalueVec = LaVectorComplex(size);
    LaGenMatComplex diagonaleigenexp = LaGenMatComplex::zeros(size, size);
    LaGenMatComplex eigenvectors = LaGenMatComplex::zeros(size, size);
    LaGenMatComplex eigenvectortrans = LaGenMatComplex::zeros(size, size);
    LaGenMatComplex UTD = LaGenMatComplex::zeros(size, size);
    //calculate eigenstuff
    matrix_eigenvstuff(matrix, eigenvalueVec, eigenvectors);
    print_matrix(eigenvalueVec, "eigenvalue vector");
    //calculate exponentials
    COMPLEX eigenexp[size];
    for(int i = 0; i < size; i++){
        test_scalar_exponential(eigenvalueVec(i), 5000, eigenexp[i]);
    }
    cout << endl;
    // save to diagonal
    array_to_diag(eigenexp, size, diagonaleigenexp);
    print_matrix(diagonaleigenexp, "exponential matrix");
    print_matrix(eigenvectors, "eigenvector matrix");
    /* calculate U^T */
    matrix_transpose(eigenvectors, size, eigenvectortrans);
    print_matrix(eigenvectortrans, "transpose eigenvector matrix");
    /* multiply results */
    Blas_Mat_Mat_Mult(eigenvectortrans, diagonaleigenexp, UTD);
    print_matrix(UTD, "UTD");
    Blas_Mat_Mat_Mult(UTD, eigenvectors, result);
    print_matrix(result, "result");
} //empty
void matrix_transpose(const LaGenMatComplex& matrix, const int matrix_size, LaGenMatComplex& result){
    result = LaGenMatComplex::zeros(matrix_size, matrix_size);
    for(int i = 0; i < matrix_size; i++){
        for(int j = 0; j < matrix_size; j++){
            result(i, j) = matrix(j, i);
        }
    }
}

/* Testing [3/3] */
void test_scalar_manipulation(const int max_rand){
    COMPLEX compA;
    generate_scalar(compA, max_rand);
    COMPLEX compB;
    generate_scalar(compB, max_rand);
    int realB;
    generate_scalar(realB, max_rand);
    COMPLEX result;

    cout << "compA = " << compA << endl;
    cout << "compB = " << compB << endl;
    cout << "realB = " << realB << endl;
    cout << endl;

    for(int i = 1; i < 5; i++){                                                 //factorials
        cout << "factorial(" << i << "): " << factorial(i) << endl;
    }

    scalar_addition(compA, compB, result);                                      //addition/subtraction
    cout << "scalar addition: " << result << endl << endl;

    scalar_multiplication(compA, realB, result);                                //multiplication
    cout << "scalar multiplication by scalar: " << result << endl << endl;

    scalar_multiplication(compA, compB, result);
    cout << "scalar multiplication by complex: " << result << endl << endl;

    scalar_division(compA, realB, result);         //division
    cout << "scalar division by scalar: " << result << endl << endl;

    scalar_division(compA, compB, result);
    cout << "scalar division by complex: " << result << endl << endl;

    for(int i = 1; i < 5; i++){
        scalar_powers(compA, i, result);
        cout << "scalar powers - A^" << i << " = " << result << endl;
    }
    cout << endl;

    COMPLEX sum;
    sum.r = 0;
    sum.i = 0;
    for(int i = 0; i < 5; i++){
        cout << "sum(" << i << ") = " << sum << endl;
        scalar_addition(sum, compA);
    }
    cout << "sum = " << sum << endl;

}                                       //to test
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
void test_inverse(const LaGenMatComplex& initialMatrix, const int size){
    LaGenMatComplex inverseMatrix;
    inverseMatrix = initialMatrix.copy();
    matrix_inverse(inverseMatrix, size);
    print_matrix(inverseMatrix, "inverse matrix");
}                 //to test
void test_scalar_sum(const int max_rand, const int iterations){
    COMPLEX number, step;
    generate_scalar(number, max_rand);
    step = number;
    for(int i = 0; i < iterations; i++){
        cout << step << endl;
        scalar_addition(number, step);
    }
    cout << number << endl;
}
void test_scalar_product(const int max_rand, const int iterations){
    //...
}
void test_scalar_exponential(const int iterations, const int max_rand){
    COMPLEX number, result;
    generate_scalar(number, max_rand);
    cout << endl << "scalar exponential test no.: " << number << endl << endl;
    scalar_exponential_main(number, iterations, result);
    cout << "e^" << number << " = " << result << endl;
}
void test_scalar_exponential(COMPLEX& number, const int iterations, COMPLEX& result){
    scalar_exponential_main(number, iterations, result);
    cout << "e^" << number << " = " << result << endl;
}
void test_matrix_multiplication(const int matrix_size, const int max_rand){
    int matrix_volume = matrix_size * matrix_size;
    /* generate matrix A */
    COMPLEX elementsA[matrix_volume];
    generate_array(elementsA, matrix_volume, max_rand);
	LaGenMatComplex matrixA = LaGenMatComplex(elementsA, matrix_size, matrix_size, false );
    print_matrix(matrixA, "Matrix A");
    /* generate matrix B */
    COMPLEX elementsB[matrix_volume];
    generate_array(elementsB, matrix_volume, max_rand);
	LaGenMatComplex matrixB = LaGenMatComplex(elementsB, matrix_size, matrix_size, false );
    print_matrix(matrixB, "Matrix B");
    /* generate matrix A^T */
    LaGenMatComplex transposeA;
    matrix_transpose(matrixA, matrix_size, transposeA);
    print_matrix(transposeA, "transpose A");
    /* generate matrix B^T */
    LaGenMatComplex transposeB;
    matrix_transpose(matrixB, matrix_size, transposeB);
    print_matrix(transposeB, "transpose B");
    /* initialise result */
    LaGenMatComplex result = LaGenMatComplex::zeros(matrix_size, matrix_size);
    /* A * B */
    Blas_Mat_Mat_Mult(matrixA, matrixB, result);
    print_matrix(result, "Matrix A * Matrix B");
    /* A^T * B */
    Blas_Mat_Mat_Mult(transposeA, matrixB, result);
    print_matrix(result, "Matrix A^T * Matrix B");
    /* A * B^T */
    Blas_Mat_Mat_Mult(matrixA, transposeB, result);
    print_matrix(result, "Matrix A * Matrix B^T");
    /* A^T * B^T */
    Blas_Mat_Mat_Mult(transposeA, transposeB, result);
    print_matrix(result, "Matrix A^T * Matrix B^T");
}
void test_matrix_exponential(const int size, const int max_rand){
    int matrix_volume = matrix_size * matrix_size;
    LaGenMatComplex matrix;
    LaGenMatComplex result;
    COMPLEX elements[matrix_volume];
    /* initialise stuff */
    generate_array(elements, matrix_volume, max_rand);
	matrix = LaGenMatComplex(elements, matrix_size, matrix_size, false );
    print_matrix(matrix, "initial matrix");
    result = LaGenMatComplex::zeros(size, size);
    /* calculate exponential */
    matrix_exponential(matrix, size, result);
    print_matrix(result, "e^(matrix)");
}
void test_idenpotent_exponential(){
    // Generate the matrix
    int elements [] = {2, -2, -4, -1, 3, 4, 1, -2, -3};
    COMPLEX comp[9];
    for(int i = 0; i < 9; i++){
        comp[i].r = elements[i];
        comp[i].i = 0;
    }
    LaGenMatComplex initialMatrix = LaGenMatComplex(comp, 3, 3, false );
    cout << "initialMatrix = " << endl << initialMatrix << endl;
    //calculate the exponential
    test_matrix_exponential(initialMatrix, 3);
}

/* Main Program */
int main(){
//	int matrix_size = 2, max_rand = 9;
//    int matrix_volume = matrix_size * matrix_size;

	/* generate the matrix */
//    COMPLEX elements[matrix_volume];
//    generate_array(elements, matrix_volume, max_rand);
//	LaGenMatComplex initialMatrix = LaGenMatComplex(elements, matrix_size, matrix_size, false );
//    print_matrix(initialMatrix, "initial matrix");

    /* test eigenvalues */
//    test_eigenvalues(initialMatrix, matrix_size);

    /* test scalar manipulation */
//    test_scalar_manipulation(4);

    /* test matrix multiplication */
//    test_matrix_multiplication(matrix_size, max_rand);

    /* test scalar exponentials */
//    test_scalar_exponential(5000,40);
//    test_idenpotent_exponential();

    /* test matrix exponentials */
    test_matrix_exponential(initialMatrix, matrix_size);

    /* test inversion */
//    test_inverse(initialMatrix, matrix_size);
}
