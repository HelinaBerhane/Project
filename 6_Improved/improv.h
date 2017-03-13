#ifndef COM_MAT_H
#define COM_MAT_H
#define LA_COMPLEX_SUPPORT

#include <gmc.h> 	//LaGenMatComplex
#include <lavc.h>  //LaVectorComplex
#include <string>

// n.b. string -> std::string
void print_scalar(const COMPLEX scalar);
void print_scalar(const COMPLEX scalar, const std::string name);
void print_scalar(const double scalar, const std::string name);
void print_array(const COMPLEX array[], int array_size, const std::string name);
void print_vector(const LaVectorComplex& vector, const std::string name);
void print_matrix(const LaGenMatComplex& matrix);
void print_matrix(const LaGenMatComplex& matrix, const std::string name);
void print_initial_parameters(double U, double beta, double lambda, double delta_tau, int time_size, int lattice_size);
int random_int(const int max_rand);
double random_double();
void array_to_diag(const COMPLEX array[], const int array_size, LaGenMatComplex& diag);
void vec_to_array(const LaVectorComplex& vector, const int array_size, COMPLEX array[]);
void vec_to_diag(const LaVectorComplex& vector, const int array_size, LaGenMatComplex& diag);
void clear_scalar(COMPLEX& scalar);
void clear_storage(COMPLEX storage[], const int storage_size);
void store_matrix(const LaGenMatComplex& matrix, const int matrix_number, const int matrix_size, COMPLEX storage[], const int storage_size);
void isolate_row(const LaGenMatComplex& matrix, const int matrix_width, const int row, COMPLEX array[]);
void flip_scalar(COMPLEX& spin);
void flip_spin(LaGenMatComplex& lattice, const int t, const int l);
void flip_spin_v(LaGenMatComplex& lattice, const int t, const int l);
void generate_scalar(COMPLEX& scalar, const int max_rand);
int generate_spins();
void generate_slice(const int lattice_size, COMPLEX slice[]);
void generate_lattice(const int lattice_size, const int time_size, LaGenMatComplex& lattice);
void scalar_division(const COMPLEX& A, const int B, COMPLEX& result);
void scalar_division(const COMPLEX& A, const COMPLEX& B, COMPLEX& result);
void scalar_multiplication(const COMPLEX& A, const COMPLEX& B, COMPLEX& result);
COMPLEX scalar_multiple(COMPLEX& A, const COMPLEX& B);
void scalar_product(COMPLEX& product, const double f);
void scalar_product(COMPLEX& product, const COMPLEX& number);
void scalar_sum(COMPLEX& result, const COMPLEX addition);
void scalar_exponential(const COMPLEX& number, COMPLEX& result);
void test_scalar_exponential();
void matrix_sum(const int matrix_size, LaGenMatComplex& sum, const LaGenMatComplex& matrix);
void matrix_multiple(const LaGenMatComplex& matrix, const int matrix_size, const double number, LaGenMatComplex& result);
void matrix_product(LaGenMatComplex& product, const LaGenMatComplex& matrix);
void matrix_inverse(const LaGenMatComplex& matrix, int matrix_size, LaGenMatComplex& result);
void recombine_diagonalised_matrices(const int matrix_size, LaGenMatComplex& eigenvectors, const LaVectorComplex& eigenvalues, LaGenMatComplex& result);
void matrix_exponential(const LaGenMatComplex& matrix, const int matrix_size, LaGenMatComplex& result);
void matrix_exponential_v(const LaGenMatComplex& matrix, const int matrix_size, LaGenMatComplex& result);
void matrix_negative(const int matrix_size, LaGenMatComplex& matrix);
void matrix_negative(const int matrix_size, const LaGenMatComplex& matrix, LaGenMatComplex& result);
void diagonal_matrix_exponential(const LaGenMatComplex& matrix, const int matrix_size, LaGenMatComplex& result);
void matrix_determinant_e(const int matrix_size, const LaGenMatComplex& matrix, COMPLEX& result);
void initial_parameter_calculation(const double U, const double beta, double& lambda, double& delta_tau, int& time_size);
void H_generation(const int lattice_size, LaGenMatComplex& H);
void V_calculation(const COMPLEX slice[], const int lattice_size, const double U, const double lambda, const double sigma, const double delta_tau, LaGenMatComplex& V);
void B_calculation(const COMPLEX slice[], const int lattice_size, const double U, const double lambda, const double sigma, const double delta_tau, LaGenMatComplex& B);
void B_calculation_v(const COMPLEX slice[], const int lattice_size, const double U, const double lambda, const double sigma, const double delta_tau, LaGenMatComplex& B);
void O_calculation(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double sigma, const double delta_tau, LaGenMatComplex& O);
void O_calculation_v(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double sigma, const double delta_tau, LaGenMatComplex& O);
void weight_calculation(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, COMPLEX& weight);
void weight_calculation_v(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, COMPLEX& weight);
void sweep_lattice(LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, const int iterations, double& acceptance, double& rejection);
void sweep_lattice_v(LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, const int iterations, double& acceptance, double& rejection);

#endif
