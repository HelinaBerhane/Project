#include <iostream> //cout
#include <string>
#include "qmc.h"
#include <gmd.h> 	//LaGenMatComplex
#include <laslv.h>  //LUFactorizeIP, LaLUInverseIP, etc.
#include <blas3pp.h>
#include <random>   //random_device, mt19937
#include <cstdlib>	//rand, srand
#include <math.h>

using namespace std;



/* ---- WORKING ---- */
/* -- Output -- */
void print_sites(const int array[], const int array_size){
    for(int i = 0; i < array_size; i++){
        cout.width(7);
        cout << array[i];
    }
    cout << endl;
}


/* -- Generation -- */
/* Random Numbers */
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
/* Structures */
void generate_lattice_array(const int array_size, int array[]){
    for(int i = 0; i < array_size; i++){
        array[i] = random_spin();
    }
}

/* -- Calculation -- */
/* Initial Parameters */
/* Weights */

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
    int array_size = 10, array[array_size];
    generate_lattice_array(array_size, array);
    print_sites(array, array_size);
}

/* ---- TESTING ----*/
/* Output */


/* Randomisation */


/* Generation */


/* -- Calculation -- */
double lambda_calculation(const double U){
    return acoshf(exp(sqrt(0.125*U)/2));
}
double delta_tau_calculation(const double U){
    return sqrt(0.125 / U);
}
void initial_parameter_calculation(const double beta, const double U, double lambda, double delta_tau, int time_slices){
    delta_tau = delta_tau_calculation(U);
    lambda = lambda_calculation(U);
    time_slices = beta / lambda;
}

// void V_calculation(const double lattice[], const int time_size, const double U, const double lambda, const int sigma, const double delta_tau, LaGenMatComplex& V){//should be working
//     /* initialise everything */
//     double elements[time_size];
//     // double mu = 0;
//
//     /* calculate lambda sigma s_l */
//     for(int i = 0; i < time_size; i++){
//         scalar_multiplication(lattice[i], lambda / delta_tau, elements[i]);
//         elements[i].r = elements[i].r + mu - U / 2;
//     }
//     /* given a lattice */
//     array_to_diag(elements, time_size, V);
// }


/* Testing */
void test_parameter_calculation(){
    double beta = 10, U = 1, lambda, delta_tau, time_slices;
    initial_parameter_calculation(beta, U, lambda, delta_tau, time_slices);
    cout << "lambda = " << lambda << endl;
    cout << "delta tau = " << delta_tau << endl;
    cout << "time slices = " << time_slices << endl;
}


void test_sweep(){
    /* Plan */
    /* [ ] Input */
        // [ ] lattice_size - int
        // [ ] time_slices  - int
        // [ ] iterations   - int
        // [ ] lattice      - LaGenMatComplex
        // [ ] U            - double

    /* [ ] Processing */
        // [ ] generate a lattice of spins
        // [ ] sweep the lattice

    /* [ ] Output */
        // [ ] average spins
        // [ ] acceptance probabilities

    // int lattice_size = 5, time_slices = 4;
    // int matrix_size = lattice_size * time_slices;
    // int lattice_array[matrix_size];
    // generate_lattice_array(matrix_size, lattice_array);

}
void test_increasing_U(){
    /* Plan */

    /* [ ] Input */
        // [ ] matrix_size  - int
        // [ ] iterations   - int
        // [ ] lattice      - LaGenMatComplex
        // [ ] U            - double

    /* [ ] Processing */
        // [ ] for n increasing values of U
            // [ ] generate a lattice of spins
            // [ ] sweep the lattice

    /* [ ] Output */
        // [ ] acceptance probabilities
}


/* --- Main QMC Program --- */
int main(){
    test_parameter_calculation();
}
