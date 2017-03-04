#include <iostream> //cout
#include <string>
#include "qmc.h"
#include <gmc.h> 	//LaGenMatComplex
#include <laslv.h>  //LUFactorizeIP, LaLUInverseIP, etc.
#include <blas3pp.h>
#include <random>   //random_device, mt19937
#include <cstdlib>	//rand, srand
#include <math.h>

using namespace std;

/* Randomisation [1/1]*/
int random_spin(){
    random_device rd;
    mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, 1);
    return dist(gen)*2 - 1;
}
float random_probability(){
    random_device rd;
    mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0, 1);
    return dist(gen);
}


void generate_spin_array(){
    // for()
}

void test_spin_generation(){
    for(int i = 0; i < 5; i++){
        cout.width(5);
        cout << random_spin();
    }
}
void test_sweep(){
    /* Plan */
    /* [ ] Input */
        // [ ] lattice_size - int
        // [ ] time_slices  - int
        // [ ] iterations   - int
        // [ ] lattice      - LaGenMatComplex
        // [ ] U            - float

    /* [ ] Processing */
        // [ ] generate a lattice of spins
        // [ ] sweep the lattice

    /* [ ] Output */
        // [ ] average spins
        // [ ] acceptance probabilities

    int lattice_size = 5, time_slices = 4;
    int matrix_size = lattice_size * time_slices;

}
void test_increasing_U(){
    /* Plan */

    /* [ ] Input */
        // [ ] matrix_size  - int
        // [ ] iterations   - int
        // [ ] lattice      - LaGenMatComplex
        // [ ] U            - float

    /* [ ] Processing */
        // [ ] for n increasing values of U
            // [ ] generate a lattice of spins
            // [ ] sweep the lattice

    /* [ ] Output */
        // [ ] acceptance probabilities
}

/* --- Main QMC Program --- */
int main(){
    test_spin_generation();
}
