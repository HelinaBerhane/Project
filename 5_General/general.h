#ifndef COM_MAT_H
#define COM_MAT_H
#define LA_COMPLEX_SUPPORT

#include <gmc.h> 	//LaGenMatComplex
#include <lavc.h> //LaVectorComplex
#include <string>

/* -- Output -- */
void print_initial_parameters(double U, double beta, double lambda, double delta_tau, int time_size, int lattice_size);
/* -- Generation -- */
// Random Numbers
// Structure
/* -- Calculation -- */
// Initial Parameters
void initial_parameter_calculation(const double U, const double beta, double& lambda, double& delta_tau, int& time_size);
// Matrix Operations
// Weights
// Testing
void test_initial_parameters();


#endif
