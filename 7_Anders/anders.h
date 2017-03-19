#ifndef COM_MAT_H
#define COM_MAT_H
#define LA_COMPLEX_SUPPORT

#include <gmc.h> 	//LaGenMatComplex
#include <lavc.h>  //LaVectorComplex
#include <string>

// n.b. std::string -> std::string
void print_space(const std::string file);
void print_text(const std::string text, const std::string file);
void print_scalar(const COMPLEX scalar);
void print_scalar(const COMPLEX scalar, const std::string name);
void print_scalar(const double scalar, const std::string name);
void print_scalar(const COMPLEX scalar, const std::string name, const std::string file);
void print_scalar_f(const double scalar, const std::string file);
void print_scalar(const double scalar, const std::string name, const std::string file);
void print_array(const COMPLEX array[], int array_size, const std::string name);
void print_array(const COMPLEX array[], int array_size, const std::string name, const std::string file);
void print_vector(const LaVectorComplex& vector, const std::string name);
void print_matrix(const LaGenMatComplex& matrix);
void print_matrix(const LaGenMatComplex& matrix, const std::string name);
void print_matrix(const LaGenMatComplex& matrix, const std::string name, const std::string file);
void print_matrix_f(const LaGenMatComplex& matrix, const std::string file);
void print_initial_parameters(const double U, const double beta, const double lambda, const double delta_tau, const double mu, const int time_size, const int lattice_size, const int iterations);
void print_initial_parameters(const double U, const double beta, const double lambda, double delta_tau, const double mu, const int time_size, const int lattice_size, const int iterations, const std::string file);
int random_int(const int max_rand);
double random_double();
void vec_to_diag(const LaVectorComplex& vector, const int array_size, LaGenMatComplex& diag);
void clear_scalar(COMPLEX& scalar);
void isolate_row(const LaGenMatComplex& matrix, const int matrix_width, const int row, COMPLEX array[]);
void clear_array(COMPLEX array[], const int array_size);
void flip_scalar(COMPLEX& spin);
void flip_spin(LaGenMatComplex& lattice, const int t, const int l);
int generate_spins();
void generate_slice(const int lattice_size, COMPLEX slice[]);
void generate_lattice(const int lattice_size, const int time_size, LaGenMatComplex& lattice);
std::string generate_file_name(const double U, const double beta, const int iterations, const std::string test);
int check_size(const double scalar);
void scalar_division(const COMPLEX& A, const int B, COMPLEX& result);
void scalar_multiplication(const COMPLEX& A, const COMPLEX& B, COMPLEX& result);
COMPLEX scalar_multiple(COMPLEX& A, const COMPLEX& B);
void scalar_product(COMPLEX& product, const double f);
void scalar_product(COMPLEX& product, const COMPLEX& number);
void scalar_sum(COMPLEX& result, const COMPLEX addition);
void scalar_exponential(const COMPLEX& number, COMPLEX& result);
void matrix_sum(const int matrix_size, LaGenMatComplex& sum, const LaGenMatComplex& matrix);
void matrix_multiple(const LaGenMatComplex& matrix, const int matrix_size, const double number, LaGenMatComplex& result);
void matrix_product(LaGenMatComplex& product, const LaGenMatComplex& matrix);
void matrix_negative(const int matrix_size, LaGenMatComplex& matrix);
void matrix_negative(const int matrix_size, const LaGenMatComplex& matrix, LaGenMatComplex& result);
void matrix_inverse(const LaGenMatComplex& matrix, int matrix_size, LaGenMatComplex& result);
void recombine_diagonalised_matrices(const int matrix_size, LaGenMatComplex& eigenvectors, const LaVectorComplex& eigenvalues, LaGenMatComplex& result);
void matrix_exponential(const LaGenMatComplex& matrix, const int matrix_size, LaGenMatComplex& result);
void diagonal_matrix_exponential(const LaGenMatComplex& matrix, const int matrix_size, LaGenMatComplex& result);
void triangle_matrix_v(const LaGenMatComplex& matrix, const int matrix_size, LaGenMatComplex& triangle);
void triangle_matrix(const LaGenMatComplex& matrix, const int matrix_size, LaGenMatComplex& triangle);
COMPLEX matrix_determinant(const LaGenMatComplex& matrix, const int matrix_size);
void initial_parameter_calculation(const double U, const double beta, double& lambda, double& delta_tau, int& time_size);
void H_calculation(const int lattice_size, LaGenMatComplex& H);
void V_calculation(const COMPLEX slice[], const int lattice_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& V);
void V_calculation_v(const COMPLEX slice[], const int lattice_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& V);
void V_calculation_f(const COMPLEX slice[], const int lattice_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& V, const std::string file);
void B_calculation(const COMPLEX slice[], const int lattice_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& B);
void B_calculation_v(const COMPLEX slice[], const int lattice_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& B);
void B_calculation_f(const COMPLEX slice[], const int lattice_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& B, const std::string file);
void O_calculation(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& O);
void O_calculation_v(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& O);
void O_calculation_f(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& O, const std::string file);
void weight_calculation(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, const double mu, COMPLEX& weight);
void weight_calculation_v(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, const double mu, COMPLEX& weight);
void weight_calculation_f(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, const double mu, COMPLEX& weight, const std::string file);
void weight_calculation_O(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, const double mu, LaGenMatComplex& OUP, COMPLEX& weight);
int judge_acceptance(const double probability);
void sweep_lattice(LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double beta, const double lambda, const double delta_tau, const double mu, const int iterations);
void measure_execution_time(const int iterations, const int start_s, const int stop_s, const std::string file);
void measure_result(const int count, const int acceptance, const int rejection, const std::string result, const double probability, const std::string file);
void measure_acceptance(const int acceptance, const int rejection, const int total_count, const std::string file);
void calculate_total_spin_f(const LaGenMatComplex& lattice, const int time_size, const int lattice_size);
void calculate_total_spin(const LaGenMatComplex& lattice, const int time_size, const int lattice_size, const std::string file);
void measure_spin(const LaGenMatComplex& lattice, const int time_size, const int lattice_size, const std::string file);
void measure_weight(const int count, const double probability, const COMPLEX weightBefore, const COMPLEX weightAfter, const std::string file);
void measure_av_weight();
void measure_double_occcupancy(const LaGenMatComplex& O, const int lattice_size, const std::string file);
double double_occupancy_ii(const int i, const LaGenMatComplex& O, const int lattice_size);
void measure_double_occcupancy_ii(const int i, const LaGenMatComplex& O, const int lattice_size, const std::string file);
double n(const LaGenMatComplex& O, const int lattice_size);
void measure_n(const LaGenMatComplex& O, const int lattice_size, const std::string file);

#endif
