#ifndef COM_MAT_H
#define COM_MAT_H
#define LA_COMPLEX_SUPPORT

#include <gmc.h> 	//LaGenMatComplex
#include <lavc.h> //LaVectorComplex
#include <string>
#include <iostream> //cout
#include <string>
#include "complex_matrices.h"
#include <gmc.h> 	//LaGenMatComplex
#include <laslv.h>  //LUFactorizeIP, LaLUInverseIP, etc.
#include <blas3pp.h>
#include <random>   //random_device, mt19937
#include <cstdlib>	//rand, srand
#include <math.h>

using namespace std;

/* Total [35/35] - QMC [3/3] */

/* Randomisation [1/1]*/
int basic_random_int(int max_rand){
    return rand() % (max_rand+1);
}//working
float basic_random_float(){//fix later
    return basic_random_int(1000)/1000;
}
float random_float(float min, float max){
    random_device rd;
    mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}
/* QMC */
float random_probability(){
    random_device rd;
    mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    return dis(gen);
}

/* Printing [7/7] */
void print_scalar(const COMPLEX scalar){
    cout << scalar << endl;
}//working
void print_scalar(const COMPLEX scalar, const string name){
    cout << name << ": " << scalar << endl;
}//working
void print_array(const COMPLEX array[], int len){
    for(int i = 0; i < len; i++){
        cout.width(7);
        cout << array[i] << " ";
    }
    cout << endl;
}//working
void print_array(const COMPLEX array[], int len, const string name){
	cout << name << ":" << endl;
    for(int i = 0; i < len; i++){
        cout << array[i] << endl;
    }
    cout << endl;
}//working
void print_vector(const LaVectorComplex& vector, const string name){
    cout << name << ":" << endl << vector << endl;
}//working
void print_matrix(const LaGenMatComplex& matrix){
	cout << matrix << endl;
}//working
void print_matrix(const LaGenMatComplex& matrix, const string name){
	cout << name << ":" << endl << matrix << endl;
}//working

/* Generation [5/5]*/
void generate_scalar(COMPLEX& scalar, const int max_rand){
    scalar.r = basic_random_int(max_rand);	//1 to x
    scalar.i = basic_random_int(max_rand);
}//working
void generate_scalar(int scalar, const int max_rand){
    scalar = basic_random_int(max_rand);	//1 to x
}//working
void generate_array(COMPLEX array[], const int array_length, const int max_rand){
    for(int i = 0; i < array_length; i++){
        array[i].r = basic_random_int(max_rand);	//1 to x
        array[i].i = basic_random_int(max_rand);
	}
}//working
void generate_matrix(const int matrix_size, const int max_rand, LaGenMatComplex& matrix){
    int matrix_volume = matrix_size*matrix_size;
    COMPLEX elements[matrix_volume];
    generate_array(elements, matrix_volume, max_rand);
    matrix = LaGenMatComplex(elements, matrix_size, matrix_size, false);
}//working
void generate_cofactor_matrix(const int matrix_size, const LaGenMatComplex& matrix, const int element, LaGenMatComplex& cofactorMatrix){
    for(int r = 1; r < matrix_size; r++){ // skip first row
        int newC = 0;
        for(int c = 0; c < matrix_size; c++){
            if(c != element){ // slip column
                cofactorMatrix(r - 1, newC).r = matrix(r, c).r;
                cofactorMatrix(r - 1, newC).i = matrix(r, c).i;
                newC++;
            }
        }
    }
}//working
// QMC - [4/4]
int generate_spins(){
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dist(0, 1);
    return (dist(gen) % 2)*2 - 1;
}//working
void generate_lattice(const int matrix_size, LaGenMatComplex& lattice){
    int matrix_volume = matrix_size * matrix_size;
    COMPLEX elements[matrix_volume];
    for(int i = 0; i < matrix_volume; i++){
        elements[i].r = generate_spins();
        elements[i].i = 0;
    }
    lattice = LaGenMatComplex(elements, matrix_size, matrix_size, false);
}//working
void generate_H(const int matrix_size, LaGenMatComplex& H){
    int matrix_volume = matrix_size * matrix_size;
    COMPLEX elements[matrix_volume];
    int n;
    /* generate the matrix */
    for(int i = 0; i < matrix_size; i++){
        for (int j = 0; j < matrix_size; j++) {
            n = (matrix_size * i) + j;
            //cout.width(3);
            //cout << abs(i-j);
            if(abs(i-j) == 1 || abs(i-j) == matrix_size - 1){
                elements[n].r = -1;
            }else{
                elements[n].r = 0;
            }
            elements[n].i = 0;
        }
        //cout << endl;
    }
    //cout << endl;
    H = LaGenMatComplex(elements, matrix_size, matrix_size, false );
    /* print result */
}//working
void generate_lattice_array(const int matrix_size, COMPLEX elements[]){
    for(int i = 0; i < matrix_size; i++){   // for each element,
        elements[i].r = generate_spins();   // generate real random spin
        elements[i].i = 0;
    }
}//working
void generate_lattice_matrix(const int matrix_size, LaGenMatComplex& lattice){
    lattice = LaGenMatComplex::zeros(matrix_size, matrix_size);
    for(int i = 0; i < matrix_size; i++){
        for(int j = 0; j < matrix_size; j++){
            lattice(i, j).r = generate_spins();
        }
    }
}//working

/* Matrix conversion [5/5] */
void vec_to_array(const LaVectorComplex& vector, const int len, COMPLEX array[]){
    for(int i = 0; i < len; i++){
        array[i] = vector(i);
    }
}//working
void array_to_diag(const COMPLEX array[], const int len, LaGenMatComplex& diag){
    diag = 0;
    for(int i = 0; i < len; i++){
        diag(i, i) = array[i];
    }
}//working
void vec_to_diag(const LaVectorComplex& vector, const int len, LaGenMatComplex& diag){
    COMPLEX array[len];
    vec_to_array(vector, len, array);
    array_to_diag(array, len, diag);
}//working
void copy_array(const int len, const COMPLEX array[], COMPLEX copy[]){//in progress
    for(int i = 0; i < len; i++){
        //
    }
}
void isolate_row(const LaGenMatComplex& matrix, const int len, const int row, COMPLEX array[]){
    for(int i = 0; i < len; i++){
        array[i] = matrix(row, i);
    }
}

/* Scalar manipulation [14/14] */
int factorial(int x){
	if(x <= 1){
        return 1;
	}else{
        return x * factorial(x - 1);
	}
}//working
void copy_scalar(const COMPLEX& scalar, COMPLEX& copy){//should work
    copy.r = scalar.r;
    copy.i = scalar.i;
}
void copy_negative_scalar(const COMPLEX& scalar, COMPLEX& copy){//should work
    copy.r = -scalar.r;
    copy.i = -scalar.i;
}
void flip_scalar(COMPLEX& spin){//should work
    spin.r = -spin.r;
    spin.i = -spin.i;
}
void scalar_addition(const COMPLEX& A, const COMPLEX& B, COMPLEX& result){
    result.r = A.r + B.r;
    result.i = A.i + B.i;
}//working
void scalar_sum(COMPLEX& result, const COMPLEX addition){//probably working
    result.r += addition.r;
    result.i += addition.i;
}//working
void scalar_multiplication(const COMPLEX& A, const int B, COMPLEX& result){//to test
    result.r = A.r * B;
    result.i = A.i * B;
}//working
void scalar_multiplication(const COMPLEX& A, const COMPLEX& B, COMPLEX& result){
    la::complex<double> laA = la::complex<double>(A); //convert to la::complex<double>
    la::complex<double> laB = la::complex<double>(B);
    la::complex<double> laResult = la::complex<double>(result);
    laResult = laA * laB;
    result = laResult.toCOMPLEX();
}//working
void scalar_product(COMPLEX& product, const COMPLEX& number){
    COMPLEX part;
    part.r = (product.r * number.r) - (product.i * number.i);
    part.i = (product.r * number.i) + (product.i * number.r);
    product = part;
}//working
COMPLEX scalar_multiple(COMPLEX& A, const COMPLEX& B){
    COMPLEX part;
    part.r = (A.r * B.r) - (A.i * B.i);
    part.i = (A.r * B.i) + (A.i * B.r);
    return part;
}
void scalar_division(const COMPLEX& A, const int B, COMPLEX& result){
    result.r = A.r / B;
    result.i = A.i / B;
}//working
void scalar_division(const COMPLEX& A, const COMPLEX& B, COMPLEX& result){
    la::complex<double> laA = la::complex<double>(A); //convert to la::complex<double>
    la::complex<double> laB = la::complex<double>(B);
    la::complex<double> laResult = la::complex<double>(result);
    laResult = laA / laB;
    result = laResult.toCOMPLEX();
}//working
void scalar_powers(const COMPLEX& number, const int power, COMPLEX& result);
void scalar_exponential_main(const COMPLEX& number, const int iterations, COMPLEX& result);
void vector_exponential(const LaVectorComplex& vector, const int matrix_size, const int iterations, LaVectorComplex& result);
void matrix_negative(const int matrix_size, LaGenMatComplex& matrix);
void matrix_negative(const int matrix_size, const LaGenMatComplex& matrix, LaGenMatComplex& result);
void matrix_sum(const int matrix_size, LaGenMatComplex& sum, const LaGenMatComplex& matrix);
void matrix_eigenvstuff(const LaGenMatComplex& matrix, LaVectorComplex& eigenvalues, LaGenMatComplex& eigenvectors);
void recombine_diagonalised_matrices(const int matrix_size, LaGenMatComplex& eigenvectors, const LaVectorComplex& eigenvalues, LaGenMatComplex& result);
void matrix_inverse(LaGenMatComplex& matrix, int matrix_size);
void matrix_exponential(const LaGenMatComplex& matrix, const int matrix_size, const int iterations, LaGenMatComplex& result);
void diagonal_matrix_exponential(const LaGenMatComplex& matrix, const int matrix_size, const int iterations, LaGenMatComplex& result);
void matrix_transpose(const LaGenMatComplex& matrix, const int matrix_size, LaGenMatComplex& result);
void matrix_product(LaGenMatComplex& product, const LaGenMatComplex& matrix);
void five_matrix_multiplication(const LaGenMatComplex& matrixA, const LaGenMatComplex& matrixB, const LaGenMatComplex& matrixC, const LaGenMatComplex& matrixD, const LaGenMatComplex& matrixE, LaGenMatComplex& result);
COMPLEX simple_matrix_determinant(const LaGenMatComplex& matrix);
COMPLEX determinant_coefficient(const LaGenMatComplex& matrix, const int element);
COMPLEX my_matrix_determinant(const int matrix_size, const LaGenMatComplex& matrix);
void matrix_determinant(const int matrix_size, const LaGenMatComplex& matrix, COMPLEX& result);
float lambda_calculation(const float U);
float delta_tau_calculation(const float U);
void V_calculation(const COMPLEX lattice[], const int time_size, const float U, const float lambda, const float delta_tau, LaGenMatComplex& V);
void B_calculation(LaGenMatComplex& H, LaGenMatComplex& V, LaGenMatComplex& B, const int matrix_size, const int iterations);
void O_calculation(const int matrix_size, const LaGenMatComplex& BA, const LaGenMatComplex& BB, const LaGenMatComplex& BC, const LaGenMatComplex& BD, const LaGenMatComplex&BE, LaGenMatComplex& O);
void detO_calculation(const int matrix_size, const LaGenMatComplex& O, COMPLEX& detO);
void calculate_weight(const int matrix_size, const COMPLEX latticeUP[], const float U, const float lambda, const float delta_tau, COMPLEX& weight);
void sweep_lattice(const int matrix_size, LaGenMatComplex& lattice, const float U, const int iterations);
void test_random_int();
void test_random_float();
void test_negative_scalar();
void test_scalar_manipulation(const int max_rand);
void test_eigenvalues(const int matrix_size, const int max_rand);
void test_inverse(const LaGenMatComplex& initialMatrix, const int matrix_size);
void test_scalar_sum(const int max_rand, const int iterations);
void test_scalar_exponential(const int max_rand, const int iterations);
void test_scalar_exponential(COMPLEX& number, const int iterations, COMPLEX& result);
void test_matrix_subtraction(const int matrix_size, const int max_rand);
void test_matrix_multiplication(const int matrix_size, const int max_rand);
void test_matrix_product(const int matrix_size, const int max_rand);
void test_five_matrix_multiplication(const int matrix_size, const int max_rand);
void test_matrix_exponential(const int matrix_size, const int max_rand, const int iterations);
void test_idenpotent_exponential(const int iterations);
void test_diagonal_exponential(const int iterations);
void test_simple_matrix_determinant(const int max_rand);
void test_determinant_coefficient();
void test_reduced_matrix();
void test_matrix_determinant();
void test_isolate_row();
void test_random_probability();
void test_lattice_generation();
void test_parameter_calculation();
void test_H(const int matrix_size);
void test_V_generation();
void test_B_generation();
void test_O_generation(const int time_size, const int iterations);
void test_detO();
void test_weight();
void test_sweep();
void test_increasing_U();

/* --- Main QMC Program --- */
int main(){

    cout << "---- TESTING INCREASING U ----" << endl;
    test_increasing_U();
    /* notes */

}

#endif
