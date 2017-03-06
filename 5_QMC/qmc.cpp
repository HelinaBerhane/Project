#include <iostream> //cout
#include <string>
#include "qmc.h"
#include <gmd.h> 	//LaGenMatDouble
#include <laslv.h>  //LUFactorizeIP, LaLUInverseIP, etc.
#include <blas3pp.h>
#include <random>   //random_device, mt19937
#include <cstdlib>	//rand, srand
#include <math.h>

using namespace std;

/* ------ WORKING ------ */

/* -- Output -- */
void print_matrix(const LaGenMatDouble& matrix){
	cout << matrix << endl;
}
void print_matrix(const LaGenMatDouble& matrix, const string name){
	cout << name << ":" << endl << matrix << endl;
}
void print_sites(const double array[], const int array_size){
    for(int i = 0; i < array_size; i++){
        cout.width(7);
        cout << array[i];
    }
    cout << endl;
}
void print_vector(const LaVectorComplex& vector, const string name){
    cout << name << ":" << endl << vector << endl;
}
/* -- Generation -- */
// Random Numbers
int random_spin(){
    random_device rd;
    mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, 1);
    return dist(gen)*2 - 1;
}
double random_probability(){
    random_device rd;
    mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0, 1);
    return dist(gen);
}
// Structures
void generate_lattice_array(double array[], const int array_size){
    for(int i = 0; i < array_size; i++){
        array[i] = random_spin();
    }
}


/* -- Calculation -- */
// Initial Parameters
void initial_parameter_calculation(const double U, const double beta, double& lambda, double& delta_tau, double& time_slices){
    delta_tau = sqrt(0.125 / U);
    lambda = acosh(exp(sqrt(0.125*U)/2));
    time_slices = beta / lambda;
}
// Weights
void generate_H(const int matrix_size, LaGenMatDouble& H){
    /* initialise everything */
    int matrix_volume = matrix_size * matrix_size, i;
    double elements[matrix_volume];

    /* generate the matrix */
    for(int row = 0; row < matrix_size; row++){
        for (int column = 0; column < matrix_size; column++) {
            i = (matrix_size * row) + column;
            if(abs(row - column) == 1 || abs(row - column) == matrix_size - 1){
                elements[i] = -1;
            }else{
                elements[i] = 0;
            }
        }
    }
    H = LaGenMatDouble(elements, matrix_size, matrix_size, false);

    /* print the matrix */
    print_matrix(H, "H");
}
void test_H(){
    /* initialise everything */
    int matrix_size = 5;
    LaGenMatDouble H;
    LaVectorComplex eigenvalues = LaVectorComplex(matrix_size);
    LaGenMatDouble eigenvectors = LaGenMatDouble::zeros(matrix_size, matrix_size);
    /* generate matrix */
    generate_H(matrix_size, H);
    print_matrix(H);
    /* calculate eigenstuff */
    LaEigSolve(H, eigenvalues, eigenvectors);
    print_vector(eigenvalues, "eigenvalues");
    // eigenvalues are 2 cos(n pi / q), where q = the matrix size
}

/* -- Testing -- */
void test_spin_generation(){
    for(int i = 0; i < 5; i++){
        cout.width(10);
        cout << random_spin();
    }
    cout << endl;
}
void test_probability_generation(){
    for(int i = 0; i < 5; i++){
        cout.width(5);
        cout << random_probability();
    }
    cout << endl;
}
void test_lattice_array_generation(){
    int array_size = 10;
    double array[array_size];
    generate_lattice_array(array, array_size);
    print_sites(array, array_size);
}
void test_parameter_calculation(){
    double beta = 10, U = 1, lambda, delta_tau, time_slices;
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_slices);
    cout << "lambda = " << lambda << endl;
    cout << "delta tau = " << delta_tau << endl;
    cout << "time slices = " << time_slices << endl;
}
void test_matrix_storage(){
    /* initialise everything */
    int array_width = 5, array_size = array_width * array_width;
    double array[array_size];
    /* generate the lattice */
    for(int i = 0; i < array_size; i++){
        array[i] = i;
    }
    LaGenMatDouble matrix(array, array_width, array_width, true);
    /* test the results */
    print_matrix(matrix, "true");
    array[5] = -array[5];
        // row_ordering = false does change the matrix, so they're linked
        // row_ordering = tru does not change the matrix, so they're copied
    print_matrix(matrix, "true (changed)");

}


/* ------ TO TEST ------ */
// COMPLEX -> double
// float -> double
// LaVectorComplex -> LaGenMatDouble
// .r ->
// .i ->

void array_to_diag(const double array[], const int len, LaGenMatDouble& diag){
    diag = 0;
    for(int i = 0; i < len; i++){
        diag(i, i) = array[i];
    }
}
void V_calculation(const double time_slice[], const int time_size, const double U, const double lambda, const double sigma, const double delta_tau, LaGenMatDouble& V){//should be working
    /* initialise everything */
    double V_elements[time_size], mu = 0;
    /* V_ii = lambda sigma s_l / delta_tau + (mu - U/2) */
    for(int i = 0; i < time_size; i++){
        V_elements[i] = lambda * sigma * time_slice[i] / delta_tau + mu - U/2;
    }
    /* save to diagonal matrix */
    array_to_diag(V_elements, time_size, V);
}
void test_V_generation(){//should work

    /* initialise everything */
    int matrix_size = 5, placeholder_time_slices = 5;
    double time_slice[matrix_size], U = 1, beta = 10, lambda, sigma, delta_tau, time_slices;
    LaGenMatDouble V = LaGenMatDouble::zeros(matrix_size, matrix_size);

    /* calculate initial parameters */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_slices);

    /* generate the time_slice */
    generate_lattice_array(matrix_size, placeholder_time_slices);
    V_calculation(time_slice, matrix_size, U, lambda, sigma, delta_tau, V);

    /* print result */
    print_matrix(V, "V");
}


/* ------ Main QMC Program ------ */
int main(){
    test_parameter_calculation();
}
