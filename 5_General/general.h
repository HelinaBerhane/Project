#ifndef COM_MAT_H
#define COM_MAT_H
#define LA_COMPLEX_SUPPORT

#include <gmc.h> 	//LaGenMatComplex
#include <lavc.h> //LaVectorComplex
#include <string>

// n.b. string -> std::string

/* -- Output -- */
void print_vector(const LaVectorComplex& vector, const std::string name);
void print_matrix(const LaGenMatComplex& matrix);
void print_matrix(const LaGenMatComplex& matrix, const std::string name);
void print_initial_parameters(double U, double beta, double lambda, double delta_tau, int time_size, int lattice_size);
/* -- Processing -- */
int generate_spins();
void generate_lattice(const int lattice_size, const int time_size, LaGenMatComplex& lattice);
void initial_parameter_calculation(const double U, const double beta, double& lambda, double& delta_tau, int& time_size);
/* -- Testing -- */
void test_initial_parameters();
void test_generate_lattice();

#endif
