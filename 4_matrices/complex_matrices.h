#ifndef COM_MAT_H
#define COM_MAT_H
#define LA_COMPLEX_SUPPORT

#include <gmc.h> 	//LaGenMatComplex
#include <lavc.h> //LaVectorComplex
#include <string>

void print_scalar(const COMPLEX scalar);
void print_scalar(const COMPLEX scalar, const std::string name);
void print_array(const COMPLEX array[], int len);
void print_array(const COMPLEX array[], int len, const std::string name);
void print_matrix(const LaGenMatComplex& matrix);
void print_matrix(const LaGenMatComplex& matrix, const std::string name);
void generate_scalar(COMPLEX& A, const int x);
void generate_scalar(int number, const int x);
void generate_array(COMPLEX array[], const int len, const int x);
void vec_to_array(const LaVectorComplex& vector, const int len, COMPLEX array[ ]);
void array_to_diag(COMPLEX array[], const int len, LaGenMatComplex& diag);
void vec_to_diag(const LaVectorComplex& vector, const int len, LaGenMatComplex& diag);
int factorial(int x);
void scalar_addition(const COMPLEX& A, const COMPLEX& B , COMPLEX& result);
void scalar_addition(COMPLEX& result, const COMPLEX addition);
void scalar_multiplication(const COMPLEX& A, const int B, COMPLEX& result);
void scalar_multiplication(const COMPLEX& A, const COMPLEX& B, COMPLEX& result);
void scalar_product(COMPLEX& total, const COMPLEX& number);
void scalar_division(const COMPLEX& A, const int B, COMPLEX& result);
void scalar_division(const COMPLEX& A, const COMPLEX& B, COMPLEX& result);
void scalar_powers(const COMPLEX& number, const int power, COMPLEX& result);
void scalar_exponential_main(const COMPLEX& number, const int iterations, COMPLEX& result);
void scalar_exponential(const COMPLEX& number, const int iter, COMPLEX& result);
void vector_exponential(const LaVectorComplex& vector, const int matrix_size, const int iterations, LaVectorComplex& result);
void recursive_scalar_exponential(const COMPLEX& number, const int iter, COMPLEX& result);
void array_powers(COMPLEX array[], const int len, const int power);
void matrix_eigenvstuff(const LaGenMatComplex& matrix, LaVectorComplex& eigenvalues, LaGenMatComplex& eigenvectors);
void recombine_diagonalised_matrices(const int matrix_size, const LaGenMatComplex& eigenvectors, const LaVectorComplex& eigenvalues, LaGenMatComplex& result);
void matrix_inverse(LaGenMatComplex& matrix, int matrix_size);
void matrix_exponential(const LaGenMatComplex& matrix, const int matrix_size, const int iterations, LaGenMatComplex& result);
void matrix_transpose(const LaGenMatComplex& matrix, const int matrix_size, LaGenMatComplex& result);
void matrix_product(LaGenMatComplex& product, const LaGenMatComplex& matrix);
void five_matrix_multiplication(const LaGenMatComplex& matrixA, const LaGenMatComplex& matrixB, const LaGenMatComplex& matrixC, const LaGenMatComplex& matrixD, const LaGenMatComplex& matrixE, LaGenMatComplex& result);
void test_scalar_manipulation(const int max_rand);
void test_eigenvalues(const int matrix_size, const int max_rand);
void test_inverse(const LaGenMatComplex& initialMatrix, const int matrix_size);
void test_scalar_sum(const int max_rand, const int iterations);
void test_scalar_product(const int max_rand, const int iterations);
void test_scalar_exponential(const int iterations, const int max_rand);
void test_scalar_exponential(COMPLEX& number, const int iterations, COMPLEX& result);
void test_matrix_multiplication(const int matrix_size, const int max_rand);
void test_matrix_product(const int matrix_size, const int max_rand);
void test_matrix_exponential(const int size, const int max_rand);
void test_five_matrix_multiplication(const int matrix_size, const int max_rand);
void test_matrix_exponential(const int matrix_size, const int max_rand, const int iterations);
void test_idenpotent_exponential(const int iterations);


#endif
