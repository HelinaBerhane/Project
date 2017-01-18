#ifndef COM_MAT_H
#define COM_MAT_H
#define LA_COMPLEX_SUPPORT

#include <gmc.h> 	//LaGenMatComplex
#include <lavc.h> //LaVectorComplex

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

#endif
