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
void print_array(const double array[], const int array_size){
    for(int i = 0; i < array_size; i++){
        cout.width(7);
        cout << array[i];
    }
    cout << endl;
}
void print_array(const double array[], int array_size, const string name){
	cout << name << ":" << endl;
    for(int i = 0; i < array_size; i++){
		cout.width(7);
        cout << array[i];
    }
    cout << endl;
}
void print_vector(const LaVectorComplex& vector, const string name){
    cout << name << ":" << endl << vector << endl;
}
void print_matrix(const LaGenMatDouble& matrix){
	cout << matrix << endl;
}
void print_matrix(const LaGenMatDouble& matrix, const string name){
	cout << name << ":" << endl << matrix << endl;
}
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
void array_to_diag(const double array[], const int array_size, LaGenMatDouble& diag){
    for(int i = 0; i < array_size; i++){
        diag(i, i) = array[i];
    }
}

/* -- Calculation -- */

// Initial Parameters
void initial_parameter_calculation(const double U, const double beta, double& lambda, double& delta_tau, int& time_size){
    // delta_tau = sqrt(0.125 / U);            // by convension
    lambda = acosh(exp(sqrt(0.125*U)/2));   // by definition
    time_size = ceil(beta / lambda);      // by definition
    delta_tau = beta / time_size;
}

// Weights
void H_generation(const int lattice_size, LaGenMatDouble& H){
    /* initialise everything */
    int matrix_volume = lattice_size * lattice_size, i;
    double elements[matrix_volume];

    /* generate the matrix */
    for(int row = 0; row < lattice_size; row++){
        for (int column = 0; column < lattice_size; column++) {
            i = (lattice_size * row) + column;
            if(abs(row - column) == 1 || abs(row - column) == lattice_size - 1){
                elements[i] = -1;
            }else{
                elements[i] = 0;
            }
        }
    }
    H = LaGenMatDouble(elements, lattice_size, lattice_size, false);

    /* print the matrix */
    print_matrix(H, "H");
}
void V_calculation(const double time_slice[], const int lattice_size, const double U, const double lambda, const double sigma, const double delta_tau, LaGenMatDouble& V){

    /* initialise everything */
    double V_elements[lattice_size];
	double mu = 0, beta = 0, time_size = 0;

    /* calculate V */
    for(int l = 0; l < lattice_size; l++){
        V_elements[l] = lambda * sigma * time_slice[l] / delta_tau + mu - U/2;

		// /* Testing */
		// print_initial_parameters(U, beta, lambda, delta_tau, time_size, lattice_size);
		//
		// cout << "lattice point = " << time_slice[l] << endl;
		// cout << "V_" << l << l << " = lambda * sigma * lattice point / delta_tau + mu - U/2 = " << endl;
		//
		// cout << "     = " << lambda << " * " << sigma << " * " << time_slice[l];
		// cout << " / " << delta_tau << " + " << mu << " - " << U << " / " << 2;
		// cout << " = " << V_elements[l] << endl << endl;
		//
		// cout << "     = " << lambda * sigma * time_slice[l] / delta_tau;
		// cout << " + " << mu << " - " << U / 2 << " = ";
		// cout << V_elements[l] << endl << endl;
    }

    /* save to diagonal matrix */
    array_to_diag(V_elements, lattice_size, V);
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
    print_array(array, array_size);
}
void test_parameter_calculation(){
    double beta = 10, U = 1, lambda, delta_tau;
    int time_size;
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    cout << "lambda = " << lambda << endl;
    cout << "delta tau = " << delta_tau << endl;
    cout << "no of time slices = " << time_size << endl;
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
void test_H(){
    /* initialise everything */
    int lattice_size = 5;
    LaGenMatDouble H;
    LaVectorComplex eigenvalues = LaVectorComplex(lattice_size);
    LaGenMatDouble eigenvectors = LaGenMatDouble::zeros(lattice_size, lattice_size);
    /* generate matrix */
    H_generation(lattice_size, H);
    print_matrix(H);
    /* calculate eigenstuff */
    LaEigSolve(H, eigenvalues, eigenvectors);
    print_vector(eigenvalues, "eigenvalues");
    // eigenvalues are 2 cos(n pi / q), where q = the matrix size
}
void test_V_generation(){//should work

    /* initialise everything */
    int lattice_size = 5, time_size;
    double time_slice[lattice_size], U = 1, beta = 10, lambda, sigma = 1, delta_tau;
    LaGenMatDouble V = LaGenMatDouble::zeros(lattice_size, lattice_size);

    /* calculate initial parameters */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);

    /* generate the time_slice */
    generate_lattice_array(time_slice, lattice_size);
    V_calculation(time_slice, lattice_size, U, lambda, sigma, delta_tau, V);

    /* print result */
	print_array(time_slice, lattice_size, "slice");
    print_matrix(V, "V");
}


						/* ------ TO TEST ------ */

// LaGenMatComplex 	  -> LaGenMatDouble
// COMPLEX 			  -> double
// float			  -> double
// .r and .i		  ->
// matrix_size 		  -> lattice_size or time_size
// len 				  -> array_size
// matrix_eigenvstuff -> LaEigSolve
// generate_H		  -> H_generation


void matrix_to_array(const LaGenMatDouble& matrix, const int matrix_size, double array[]){
	/* initialise everything */
	int matrix_volume = matrix_size * matrix_size, i = j * matrix_size + k;

	/* convert everything */
	for(int j = 0; j < matrix_size; j++){
		for(int k = 0; k < matrix_size; k++){
			array[e] = matrix(j, k);
		}
	}
}





						/* ------ TO WRITE ------ */
void generate_matrix(const int matrix_size, const int max_rand, LaGenMatDouble& matrix){
    int matrix_volume = matrix_size * matrix_size;
    COMPLEX elements[matrix_volume];
    generate_array(elements, matrix_volume, max_rand);
    matrix = LaGenMatDouble(elements, matrix_size, matrix_size, false);
}
void test_eigenvalues(const int matrix_size, const int max_rand){
    /* initialise everything */
    LaGenMatDouble matrix;
    LaVectorComplex eigenvalues = LaVectorComplex(matrix_size);
    LaGenMatDouble eigenvectors = LaGenMatDouble::zeros(matrix_size, matrix_size);
    LaGenMatDouble result;
    /* generate matrix */
    generate_matrix(matrix_size, max_rand, matrix);
    print_matrix(matrix, "initial matrix");
    /* calculate eigenstuff */
    LaEigSolve(matrix, eigenvalues, eigenvectors);
    /* multiply them back together to get the matrix */
    recombine_diagonalised_matrices(matrix_size, eigenvectors, eigenvalues, result);
    /* wolfram test */
        // 2x2 real:
        // 3x3 complex: {{1+7i, 1+3i, 5+7i},{7i, 6+i, 5+4i},{5+7i, 5+4i, 6}}
}//working

void test_matrix_product(const int matrix_size, const int max_rand){
    /* initialise everything */
    LaGenMatDouble matrixA;
    LaGenMatDouble matrixB;
    /* generate everything */
    generate_matrix(matrix_size, max_rand, matrixA);
    generate_matrix(matrix_size, max_rand, matrixB);
    /* print everything */
    print_matrix(matrixA, "Matrix A");
    print_matrix(matrixB, "Matrix B");
    /* matrix product */
    matrix_product(matrixA, matrixB);
    print_matrix(matrixB, "result");
}//working

void n_matrix_product(const double matrices[], const int matrix_size, const int n, LaGenMatDouble& product){

	/* Plan */
		// I want to multiply an arbitrary number of matrices
		// to do this, i'll store each matrix as an array
		// and convert from an array to matrix as needed
		// ideally, if i could do this with lapackpp
			// but lapack only does this with column ordered matrices
			// i only need to do this with the B matrices
				// B = exp(-H) * exp(-V)
				// H doesn't care if it's row ordered or column ordered
				// V also doesnt care
				// so neither does B?
		// so the steps are

		/* [ ] Input */
			// [ ] no of matrices		- int
			// [ ] point[matrix_size] 	- double

		/* [ ] Processing */
			// [ ] convert the first pointer to a matrix
			// [ ] the product = the first matrix
			// [ ] for each subsequent pointer (matrix)
				// [ ] convert the pointer to a matrix
				// [ ] multiply it with the ongoing product

		/* [ ] Output */
			// [ ] return the product

	/* initialise everything */

	/* multiply the matrices */





    // if(n <= 0){
    //     return;
    // }
    // Blas_Mat_Mat_Mult(product, * matrices[0], product);
    // n_matrix_product(product, matrices + 1, n - 1 );

}
void test_n_matrix_product(){

    /* initialise everything */
    int n = 3, matrix_size = 5, max_rand = 9;
    LaGenMatDouble* matrices[n];
        // this is an array of pointers
    LaGenMatDouble product = LaGenMatDouble::eye(2, 2);

    for(int i = 0; i < n; i++){

        /* generate everything */
        generate_matrix(matrix_size, max_rand, *matrices[n]);

        /* print everything */
        cout << "(" << n << ")" << endl;
        print_matrix(*matrices[n]);
    }

    /* multiply everything */
    n_matrix_product(product, matrices, n);

    /* print everything */
    print_matrix(product, "result");
}

void matrix_negative(const int matrix_size, LaGenMatDouble& matrix){

    LaGenMatDouble result = LaGenMatDouble::zeros(matrix_size, matrix_size);

    for(int i = 0; i < matrix_size; i++){
        for(int j = 0; j < matrix_size; j++){
            result(i, j) -= matrix(i, j);
        }
    }

    matrix = result.copy();
}
void matrix_negative(const int matrix_size, const LaGenMatDouble& matrix, LaGenMatDouble& result){
    result = LaGenMatDouble::zeros(matrix_size, matrix_size);
    for(int i = 0; i < matrix_size; i++){
        for(int j = 0; j < matrix_size; j++){
            result(i, j) -= matrix(i, j);
        }
    }
}

void scalar_exponential_main(const double& number, const int iterations, double& result){

    double division, total_division;
    result.r = 1;
    for(int step = 1; step <= iterations; step++){   //sum (from 1 to n)
        total_division.r = 1;
        total_division.i= 0;
        for(int i = 1; i <= step; i++){        //    ( num^n / n!)
            scalar_division(number, i, division);
            scalar_product(total_division, division);
        }
        scalar_sum(result, total_division);
    }
}//probably working

void vector_exponential(const LaVectorComplex& vector, const int matrix_size, const int iterations, LaVectorComplex& result){
    for(int i = 0; i < matrix_size; i++){
        scalar_exponential_main(vector(i), iterations, result(i));
    }
}

void matrix_exponential(const LaGenMatDouble& matrix, const int matrix_size, const int iterations, LaGenMatDouble& result){

    /* initialise everything */
    LaVectorComplex eigenvalues = LaVectorComplex(matrix_size);
    LaGenMatDouble eigenvectors = LaGenMatDouble::zeros(matrix_size, matrix_size);
    LaGenMatDouble diagonalEigenExp = LaGenMatDouble::zeros(matrix_size, matrix_size);
    LaVectorComplex eigenExponential = LaVectorComplex(matrix_size);

    /* calculate eigenstuff */
    LaEigSolve(matrix, eigenvalues, eigenvectors);

    /* calculate exponentials */
    vector_exponential(eigenvalues, matrix_size, iterations, eigenExponential);

    /* multiply them back together to get the matrix */
    recombine_diagonalised_matrices(matrix_size, eigenvectors, eigenExponential, result);
}
void test_matrix_exponential(const int matrix_size, const int max_rand, const int iterations){
    /* initialise everything */
    LaGenMatDouble matrix;
    LaGenMatDouble result;
    result = LaGenMatDouble::zeros(matrix_size, matrix_size);
    /* generate matrix */
    generate_matrix(matrix_size, max_rand, matrix);
    print_matrix(matrix, "initial matrix");
    /* calculate exponential */
    matrix_exponential(matrix, matrix_size, iterations, result);
    print_matrix(result, "e^(matrix)");
}
void test_idenpotent_exponential(const int iterations){
    /* generate the matrix */
    int numbers [] = {2, -2, -4, -1, 3, 4, 1, -2, -3};
    COMPLEX elements[9];
    for(int i = 0; i < 9; i++){
        elements[i].r = numbers[i];
        elements[i].i = 0;
    }
    LaGenMatDouble matrix = LaGenMatDouble(elements, 3, 3, false );
    LaGenMatDouble result = LaGenMatDouble::zeros(3, 3);
    //print_matrix(matrix, "initial matrix");
    /* calculate the exponential */
    for(int j = 1; j <= 5; j++){
        matrix_exponential(matrix, 3, j, result);
        cout << j << " iterations:" << endl;
        print_matrix(result);
    }
    matrix_exponential(matrix, 3, iterations, result);
    print_matrix(result, "idenpotent exponential");
}


void B_calculation(LaGenMatDouble& H, LaGenMatDouble& V, LaGenMatDouble& B, const int lattice_size, const int iterations){

	/* Plan */
		// Theres a B matrix for each time slice
		//B = exp(-H) * exp(-V)

    /* initialise everything */
    LaGenMatDouble negH;
    LaGenMatDouble negV;
    LaGenMatDouble expH;
    LaGenMatDouble expV;

    /* negate matrices (not in place) */
    matrix_negative(lattice_size, H, negH);
    matrix_negative(lattice_size, V, negV);

    /* calculate exponentials */
    matrix_exponential(negH, lattice_size, iterations, expH);
    matrix_exponential(negV, lattice_size, iterations, expV);

	/* calculate B */
    B = expH.copy();
    matrix_product(B, expV);

}
void test_B_generation(){

    /* initialise everything */
    int time_size = 5, max_rand = 9, iterations = 1000;
    LaGenMatDouble H = LaGenMatDouble::eye(time_size, time_size);
    LaGenMatDouble V = LaGenMatDouble::eye(time_size, time_size);
    LaGenMatDouble B = LaGenMatDouble::zeros(time_size, time_size);

    /* generate matrices */
    for(int i = 0; i < time_size; i++){
        H(i,i) = basic_random_int(max_rand);
        V(i,i) = basic_random_int(max_rand);
    }

    /* print matrices */
    print_matrix(H, "H");
    print_matrix(V, "V");

    /* calculate B */
    B_calculation(H, V, B, time_size, iterations);

    /* print result */
    print_matrix(B,"B = e^-H e^-V");
}

// void O_stuff(){
// 	/* initialise everything */
// 	int storage_size = lattice_size * n, element;
// 	double storage[storage_size];
// 	LaGenMatDouble B;
//
// 	// for each time slice, calculate a B lattice_size
// 	for(int t = 0; t < time_size; t++){
// 		// B_calculation(...);
//
// 		// for each lattice point, append to storage
// 		for(int l = 0; l < lattice_size, l++)
// 		element = t * lattice_size + l;
// 		storage[element] = B(t, l);
// 	}
// }

/* ------ Main QMC Program ------ */
int main(){
    test_V_generation();
}
