#include <iostream> //cout
#include <string>
#include "double.h"
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
int random_int(const int max_rand){
    random_device rd;
    mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, max_rand);
    return dist(gen);
}
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
void generate_array(double array[], const int array_length, const int max_rand){
    for(int i = 0; i < array_length; i++){
        array[i] = random_int(max_rand);
	}
}
void generate_lattice_array(double array[], const int array_size){
    for(int i = 0; i < array_size; i++){
        array[i] = random_spin();
    }
}
void generate_matrix(const int matrix_size, const int max_rand, LaGenMatDouble& matrix){
    int matrix_volume = matrix_size * matrix_size;
    double elements[matrix_volume];
    generate_array(elements, matrix_volume, max_rand);
    matrix = LaGenMatDouble(elements, matrix_size, matrix_size, true);
}
void array_to_diag(const double array[], const int array_size, LaGenMatDouble& diag){
	diag = LaGenMatDouble::eye(array_size, array_size);
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

// Matrix Operations

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
	double mu = 0;

    /* calculate V */
    for(int l = 0; l < lattice_size; l++){
        V_elements[l] = lambda * sigma * time_slice[l] / delta_tau + mu - U/2;

		// /* Testing */
		// print_initial_parameters(U, 10, lambda, delta_tau, 1, lattice_size);
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
void test_matrix_multiplication(){

	/* initialise everything */
	int matrix_size = 2, max_rand = 9;
    int matrix_volume = matrix_size * matrix_size;
	LaGenMatDouble result = LaGenMatDouble::zeros(matrix_size, matrix_size);
	LaGenMatDouble matrixA = LaGenMatDouble::zeros(matrix_size, matrix_size);
	LaGenMatDouble matrixB = LaGenMatDouble::zeros(matrix_size, matrix_size);
	double elements[matrix_volume];

    /* generate the matrices */
	generate_matrix(matrix_size, max_rand, matrixA);
	generate_matrix(matrix_size, max_rand, matrixB);

    /* A * B */
    print_matrix(matrixA, "Matrix A");
    print_matrix(matrixB, "Matrix B");
    Blas_Mat_Mat_Mult(matrixA, matrixB, result);
    print_matrix(result, "result = Matrix A * Matrix B");
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
void matrix_to_array(const LaGenMatDouble& matrix, const int matrix_size, double array[]){
	for(int j = 0; j < matrix_size; j++){
		for(int k = 0; k < matrix_size; k++){
			int i = j * matrix_size + k;
			array[i] = matrix(j, k);
		}
	}
}
void vec_to_array(const LaVectorComplex& vector, const int array_size, double array[]){
    for(int i = 0; i < array_size; i++){
        array[i] = vector(i).r;
		if(vector(i).i > 0){
			cout << "array(" << i << ") = " << vector(i).i << endl;
			cout << "Check this!" << endl;
		}
    }
}
void vec_to_diag(const LaVectorComplex& vector, const int array_size, LaGenMatDouble& diag){
    double array[array_size];
    vec_to_array(vector, array_size, array);
    array_to_diag(array, array_size, diag);
}

void matrix_inverse(LaGenMatDouble& matrix, int matrix_size){
    LaVectorLongInt PIV = LaVectorLongInt(matrix_size);
    LUFactorizeIP(matrix, PIV);
    LaLUInverseIP(matrix, PIV);
}
void matrix_product(LaGenMatDouble& product, const LaGenMatDouble& matrix){
    LaGenMatDouble result = matrix.copy();
    Blas_Mat_Mat_Mult(product, matrix, result);
    product = result.copy();
}

void test_inverse(){
	/* initialise everything */
	int matrix_size = 3, max_rand = 9;
	LaGenMatDouble matrix = LaGenMatDouble::zeros(matrix_size, matrix_size);
	LaGenMatDouble product =  LaGenMatDouble::eye(matrix_size, matrix_size);
	/* generate the matrix */
	generate_matrix(matrix_size, max_rand, matrix);
	/* calculate the inverse + test*/
	print_matrix(matrix, "initial matrix");
	matrix_product(product, matrix);
    matrix_inverse(matrix, matrix_size);
    print_matrix(matrix, "inverse matrix");
	matrix_product(product, matrix);
    print_matrix(product, "I");
}

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
}

void recombine_diagonalised_matrices(const int matrix_size, LaGenMatDouble& eigenvectors, const LaVectorComplex& eigenvalues, LaGenMatDouble& result){
    /* initialise  everything */
    LaGenMatDouble eigenvalueMatrix = LaGenMatDouble::zeros(matrix_size, matrix_size);
    LaGenMatDouble transposeEigenvectors;
    /* process matrices */
    result = eigenvectors.copy();
    vec_to_diag(eigenvalues, matrix_size, eigenvalueMatrix);
    matrix_inverse(eigenvectors, matrix_size);
    /* print matrices */
    //print_matrix(result, "U - eigenvectors (check if column based?)");
    //print_matrix(eigenvalueMatrix, "D - eigenvalues (vector)");
    //print_matrix(eigenvectors, "U^-1 - inverse eigenvectors");
    /* multiply results */
    matrix_product(result, eigenvalueMatrix);
    matrix_product(result, eigenvectors);
    /* print results */
    //print_matrix(result, "U D U^-1");
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

void test_scalar_exponential(double& number, const int iterations, double& result){
    cout << "e^" << number << " = " << exp(number) << endl;
}

void vector_exponential(const LaVectorComplex& vector, const int matrix_size, LaVectorComplex& result){
    for(int i = 0; i < matrix_size; i++){
		result(i).r = exp(vector(i).r);
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
    vector_exponential(eigenvalues, matrix_size, eigenExponential);

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

	/* initialise everything */
    double numbers[] = {2, -2, -4, -1, 3, 4, 1, -2, -3};
    LaGenMatDouble matrix = LaGenMatDouble(numbers, 3, 3, false);		// test this row ordering!!!
    LaGenMatDouble result = LaGenMatDouble::zeros(3, 3);
    print_matrix(matrix, "initial matrix");

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
        H(i,i) = random_int(max_rand);
        V(i,i) = random_int(max_rand);
    }

    /* print matrices */
    print_matrix(H, "H");
    print_matrix(V, "V");

    /* calculate B */
    B_calculation(H, V, B, time_size, iterations);

    /* print result */
    print_matrix(B,"B = e^-H e^-V");
}

							/* ------ TO CONVERT ------ */

							// LaGenMatComplex 	  -> LaGenMatDouble
							// COMPLEX 			  -> double
							// float			  -> double
							// .r and .i		  ->
							// matrix_size 		  -> lattice_size or time_size
							// len 				  -> array_size
							// matrix_eigenvstuff -> LaEigSolve
							// generate_H		  -> H_generation
							// scalar_...		  -> + - / ...
							// basic_random_int   -> random_int
							// generate_scalar	  -> random_int
							// scalar_exponential -> exp()
							// generate_spins	  -> random_spin
//

// void n_matrix_product(const double matrices[], const int matrix_size, const int n, LaGenMatDouble& product){
//
// 	/* Plan */
// 		// I want to multiply an arbitrary number of matrices
// 		// to do this, i'll store each matrix as an array
// 		// and convert from an array to matrix as needed
// 		// ideally, if i could do this with lapackpp
// 			// but lapack only does this with column ordered matrices
// 			// i only need to do this with the B matrices
// 				// B = exp(-H) * exp(-V)
// 				// H doesn't care if it's row ordered or column ordered
// 				// V also doesnt care
// 				// so neither does B?
// 		// so the steps are
//
// 		/* [ ] Input */
// 			// [ ] no of matrices		- int
// 			// [ ] point[matrix_size] 	- double
//
// 		/* [ ] Processing */
// 			// [ ] convert the first pointer to a matrix
// 			// [ ] the product = the first matrix
// 			// [ ] for each subsequent pointer (matrix)
// 				// [ ] convert the pointer to a matrix
// 				// [ ] multiply it with the ongoing product
//
// 		/* [ ] Output */
// 			// [ ] return the product
//
// 	/* initialise everything */
//
// 	/* multiply the matrices */
//
//
//
//
//
//     if(n <= 0){
//         return;
//     }
//     Blas_Mat_Mat_Mult(product, * matrices[0], product);
//     n_matrix_product(product, matrices + 1, n - 1 );
//
// }
// void test_n_matrix_product(){
//
//     /* initialise everything */
//     int n = 3, matrix_size = 5, max_rand = 9;
//     LaGenMatDouble* matrices[n];
//         // this is an array of pointers
//     LaGenMatDouble product = LaGenMatDouble::eye(2, 2);
//
//     for(int i = 0; i < n; i++){
//
//         /* generate everything */
//         generate_matrix(matrix_size, max_rand, *matrices[n]);
//
//         /* print everything */
//         cout << "(" << n << ")" << endl;
//         print_matrix(*matrices[n]);
//     }
//
//     /* multiply everything */
//     n_matrix_product(product, matrices, n);
//
//     /* print everything */
//     print_matrix(product, "result");
// }

//18:49

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
    test_inverse();
}
