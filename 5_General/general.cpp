#include <iostream> //cout
#include <string>
#include "general.h"
#include <gmd.h> 	//LaGenMatDouble
#include <laslv.h>  //LUFactorizeIP, LaLUInverseIP, etc.
#include <blas3pp.h>
#include <random>   //random_device, mt19937
#include <cstdlib>	//rand, srand
#include <math.h>

using namespace std;
						/* ------ NOTES ------ */
// be careful with row ordering
	// true = row ordered + DOES NOT link to the array it was made with
	// false = column ordered + LINKS to the array it was made with
// check that everything is reset to 0 or I where needed

						/* ------ WORKING ------ */

/* -- Output -- */
void print_initial_parameters(double U, double beta, double lambda, double delta_tau, int time_size, int lattice_size){
	cout << "no of lattice points = " << lattice_size << endl;
	cout << "no of time slices = " << time_size << endl;
	cout << "U = " << U << endl;
	cout << "beta = " << beta << endl;
	cout << "lambda = " << lambda << endl;
	cout << "delta tau = " << delta_tau << endl;
}

/* -- Generation -- */
// Random Numbers
// Structure

/* -- Calculation -- */
// Initial Parameters
void initial_parameter_calculation(const double U, const double beta, double& lambda, double& delta_tau, int& time_size){
    lambda = acoshf(exp(sqrt(0.125*U)/2));  // by definition
    time_size = ceil(beta / lambda);        // by definition
    delta_tau = beta / time_size;           // = sqrt(0.125 / U) by convension
}
// Matrix Operations
// Weights
// Testing
void test_initial_parameters(){
    double U = 1, beta = 10, lambda, delta_tau;
    int lattice_size = 5, time_size;
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, time_size, lattice_size);
}

						/* ------ TO TEST ------ */
//...

void print_matrix(const LaGenMatComplex& matrix){
	cout << matrix << endl;
}
void print_matrix(const LaGenMatComplex& matrix, const string name){
	cout << name << ":" << endl << matrix << endl;
}
void print_vector(const LaVectorComplex& vector, const string name){
    cout << name << ":" << endl << vector << endl;
}
int generate_spins(){
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dist(0, 1);
    return (dist(gen) % 2)*2 - 1;
}
void generate_lattice(const int lattice_size, const int time_size, LaGenMatComplex& lattice){
    int matrix_volume = lattice_size * time_size;
    COMPLEX elements[matrix_volume];
    for(int row = 0; row < time_size; row++){
        for(int column = 0; column < lattice_size; column++){
            int i = row * lattice_size + time_size;
            elements[i].r = generate_spins();
            elements[i].i = 0;
        }
    }
    lattice = LaGenMatComplex(elements, lattice_size, time_size, false);
}
void test_generate_lattice(){
    int lattice_size = 5, time_size = 17;
    LaGenMatComplex lattice;
    generate_lattice(lattice_size, time_size, lattice);
    print_matrix(lattice, "lattice");
}
void generate_H(const int lattice_size, LaGenMatComplex& H){
    /* initialise everything */
    int matrix_volume = lattice_size * lattice_size, n;
    COMPLEX elements[matrix_volume];

    /* generate the matrix */
    for(int i = 0; i < lattice_size; i++){
        for (int j = 0; j < lattice_size; j++) {
            n = (lattice_size * i) + j;
            //cout.width(3);
            //cout << abs(i-j);
            if(abs(i-j) == 1 || abs(i-j) == lattice_size - 1){
                elements[n].r = -1;
            }else{
                elements[n].r = 0;
            }
            elements[n].i = 0;
        }
    }

    H = LaGenMatComplex(elements, lattice_size, lattice_size, false);
}
void test_H(const int matrix_size){
    /* initialise everything */
    LaGenMatComplex H;
    LaVectorComplex eigenvalues = LaVectorComplex(matrix_size);
    LaGenMatComplex eigenvectors = LaGenMatComplex::zeros(matrix_size, matrix_size);
    /* generate matrix */
    generate_H(matrix_size, H);
    print_matrix(H);
    /* calculate eigenstuff */
    LaEigSolve(H, eigenvalues, eigenvectors);
    print_vector(eigenvalues, "eigenvalues");
    // eigenvalues are 2 cos(n pi / q), where q = the matrix size
}

						/* ------ TO CONVERT ------ */
// matrix_size          -> lattice_size or time_size
// len                  -> array_size
// matrix_eigenvstuff   -> LaEigSolve
// float                -> double

/* ------ Main QMC Program ------ */
int main(){
    test_generate_lattice();
}
