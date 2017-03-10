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
void print_array(const COMPLEX array[], int array_size, const string name){
	cout << name << ":" << endl;
    for(int i = 0; i < array_size; i++){
        cout << array[i] << endl;
    }
    cout << endl;
}
void print_vector(const LaVectorComplex& vector, const string name){
    cout << name << ":" << endl << vector << endl;
}
void print_matrix(const LaGenMatComplex& matrix){
	cout << matrix << endl;
}
void print_matrix(const LaGenMatComplex& matrix, const string name){
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

/* -- Processing -- */
// Randomisation
int random_int(const int max_rand){
    random_device rd;
    mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, max_rand);
    return dist(gen);
}
// Manipulation
void array_to_diag(const COMPLEX array[], const int array_size, LaGenMatComplex& diag){
    diag = 0;
    for(int i = 0; i < array_size; i++){
        diag(i, i) = array[i];
    }
}
void vec_to_array(const LaVectorComplex& vector, const int array_size, COMPLEX array[]){
    for(int i = 0; i < array_size; i++){
        array[i] = vector(i);
    }
}
void vec_to_diag(const LaVectorComplex& vector, const int array_size, LaGenMatComplex& diag){
    COMPLEX array[array_size];
    vec_to_array(vector, array_size, array);
    array_to_diag(array, array_size, diag);
}
// Generation
void generate_scalar(COMPLEX& scalar, const int max_rand){
    scalar.r = random_int(max_rand);
    scalar.i = random_int(max_rand);
}
int generate_spins(){
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dist(0, 1);
    return (dist(gen) % 2)*2 - 1;
}
void generate_slice(const int lattice_size, COMPLEX slice[]){
    for(int i = 0; i < lattice_size; i++){
        slice[i].r = generate_spins();
        slice[i].i = 0;
    }
}
void generate_lattice(const int lattice_size, const int time_size, LaGenMatComplex& lattice){
    int matrix_volume = lattice_size * time_size;
    COMPLEX elements[matrix_volume];
    for(int row = 0; row < time_size; row++){
        for(int column = 0; column < lattice_size; column++){
            int i = row * lattice_size + column;
            elements[i].r = generate_spins();
            elements[i].i = 0;
        }
    }
    lattice = LaGenMatComplex(elements, time_size, lattice_size, false);
}
// Calculation
// - generic
void scalar_division(const COMPLEX& A, const int B, COMPLEX& result){
    result.r = A.r / B;
    result.i = A.i / B;
}
void scalar_division(const COMPLEX& A, const COMPLEX& B, COMPLEX& result){
    la::complex<double> laA = la::complex<double>(A);
    la::complex<double> laB = la::complex<double>(B);
    la::complex<double> laResult = la::complex<double>(result);
    laResult = laA / laB;
    result = laResult.toCOMPLEX();
}
void scalar_product(COMPLEX& product, const double f){
    product.r = product.r * f;
    product.i = product.i * f;
}
void scalar_product(COMPLEX& product, const COMPLEX& number){
    COMPLEX part;
    part.r = (product.r * number.r) - (product.i * number.i);
    part.i = (product.r * number.i) + (product.i * number.r);
    product = part;
}
void scalar_sum(COMPLEX& result, const COMPLEX addition){
    result.r += addition.r;
    result.i += addition.i;
}
void scalar_exponential(const COMPLEX& number, COMPLEX& result){
    /* initialise everything */
    int iterations = 100;
    COMPLEX division, total_division;
    /* reset everything */
    result.r = 1;
    result.i = 0;
    /* calculate e^n */
    for(int step = 1; step <= iterations; step++){
        total_division.r = 1;
        total_division.i = 0;
        for(int i = 1; i <= step; i++){
            scalar_division(number, i, division);
            scalar_product(total_division, division);
        }
        scalar_sum(result, total_division);
    }
}
void test_scalar_exponential(){
    int max_rand = 9;
    COMPLEX number, result;
    generate_scalar(number, max_rand);
    cout << endl << "scalar exponential test no.: " << number << endl << endl;
    scalar_exponential(number, result);
    cout << "e^" << number << " = " << result << endl;
}
void matrix_inverse(const LaGenMatComplex& matrix, int matrix_size, LaGenMatComplex& result){
    result = matrix.copy();
    LaVectorLongInt PIV = LaVectorLongInt(matrix_size);
    LUFactorizeIP(result, PIV);
    LaLUInverseIP(result, PIV);
}
void matrix_product(LaGenMatComplex& product, const LaGenMatComplex& matrix){
    LaGenMatComplex result = matrix.copy();
    Blas_Mat_Mat_Mult(product, matrix, result);
    product = result.copy();
}
void recombine_diagonalised_matrices(const int matrix_size, LaGenMatComplex& eigenvectors, const LaVectorComplex& eigenvalues, LaGenMatComplex& result){
    /* initialise  everything */
    LaGenMatComplex eigenvalueMatrix = LaGenMatComplex::zeros(matrix_size, matrix_size);
    LaGenMatComplex inverseEigenvectors = LaGenMatComplex::zeros(matrix_size, matrix_size);
    result = LaGenMatComplex::eye(matrix_size, matrix_size);
    /* process matrices */
    vec_to_diag(eigenvalues, matrix_size, eigenvalueMatrix);
    matrix_inverse(eigenvectors, matrix_size, inverseEigenvectors);
    /* multiply results */
    matrix_product(result, eigenvectors);
    matrix_product(result, eigenvalueMatrix);
    matrix_product(result, inverseEigenvectors);
}
void matrix_exponential(const LaGenMatComplex& matrix, const int matrix_size, LaGenMatComplex& result){
    /* initialise everything */
    LaVectorComplex eigenvalues = LaVectorComplex(matrix_size);
    LaGenMatComplex eigenvectors = LaGenMatComplex::zeros(matrix_size, matrix_size);
    LaGenMatComplex diagonalEigenExp = LaGenMatComplex::zeros(matrix_size, matrix_size);
    LaVectorComplex eigenExponential = LaVectorComplex(matrix_size);
    /* calculate eigenstuff */
    LaEigSolve(matrix, eigenvalues, eigenvectors);
    /* calculate exponentials */
    for(int i = 0; i < matrix_size; i++){
        scalar_exponential(eigenvalues(i), result(i,i));
    }
    /* multiply them back together */
    recombine_diagonalised_matrices(matrix_size, eigenvectors, eigenExponential, result);
}

void matrix_negative(const int matrix_size, LaGenMatComplex& matrix){
    LaGenMatComplex result = LaGenMatComplex::zeros(matrix_size, matrix_size);
    for(int i = 0; i < matrix_size; i++){
        for(int j = 0; j < matrix_size; j++){
            result(i, j).r -= matrix(i, j).r;
            result(i, j).i -= matrix(i, j).i;
        }
    }
    matrix = result.copy();
}
void matrix_negative(const int matrix_size, const LaGenMatComplex& matrix, LaGenMatComplex& result){
    result = LaGenMatComplex::zeros(matrix_size, matrix_size);
    for(int i = 0; i < matrix_size; i++){
        for(int j = 0; j < matrix_size; j++){
            result(i, j).r -= matrix(i, j).r;
            result(i, j).i -= matrix(i, j).i;
        }
    }
}
void diagonal_matrix_exponential(const LaGenMatComplex& matrix, const int matrix_size, LaGenMatComplex& result){
    result = 0;
    for(int i = 0; i < matrix_size; i++){
        scalar_exponential(matrix(i,i), result(i,i));
    }
}
// - qmc
void initial_parameter_calculation(const double U, const double beta, double& lambda, double& delta_tau, int& time_size){
    lambda = acoshf(exp(sqrt(0.125*U)/2));  // by definition
    time_size = ceil(beta / lambda);        // by definition
    delta_tau = beta / time_size;           // = sqrt(0.125 / U) by convension
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
void V_calculation(const COMPLEX slice[], const int lattice_size, const double U, const double lambda, const double sigma, const double delta_tau, LaGenMatComplex& V){
    /* initialise everything */
    V = LaGenMatComplex::zeros(lattice_size, lattice_size);
    COMPLEX V_ii[lattice_size];

    /* V_ii = (lambda sigma s_l / delta_tau) + mu - U / 2 */
    for(int i = 0; i < lattice_size; i++){
        V_ii[i].r = lambda * sigma * slice[i].r / delta_tau - U / 2;
        V_ii[i].i = 0;
    }

    /* plot to diagonal */
    array_to_diag(V_ii, lattice_size, V);
}

/* -- Testing -- */
// - generic
void test_inverse(){
    int matrix_size = 3;
    LaGenMatComplex initialMatrix = LaGenMatComplex::rand(matrix_size, matrix_size, 0, 9);
    LaGenMatComplex inverseMatrix = LaGenMatComplex::zeros(matrix_size, matrix_size);
    matrix_inverse(initialMatrix, matrix_size, inverseMatrix);
    print_matrix(initialMatrix, "inverse matrix");
    print_matrix(inverseMatrix, "inverse matrix");
}
void test_matrix_product(){
    /* initialise everything */
    int matrix_size = 5;
    LaGenMatComplex matrixA = LaGenMatComplex::rand(matrix_size, matrix_size, 0, 9);
    LaGenMatComplex matrixB = LaGenMatComplex::rand(matrix_size, matrix_size, 0, 9);
    /* print everything */
    print_matrix(matrixA, "Matrix A");
    print_matrix(matrixB, "Matrix B");
    /* matrix product */
    matrix_product(matrixA, matrixB);
    print_matrix(matrixB, "result");
}
void test_recombine_diagonalised_matrices(){
    /* initialise everything */
    int matrix_size = 5;
    LaGenMatComplex matrix = LaGenMatComplex::rand(matrix_size, matrix_size, 0, 9);
    LaVectorComplex eigenvalues = LaVectorComplex(matrix_size);
    LaGenMatComplex eigenvectors = LaGenMatComplex::zeros(matrix_size, matrix_size);
    LaGenMatComplex result = LaGenMatComplex::zeros(matrix_size, matrix_size);
    /* generate matrix */
    print_matrix(matrix, "initial matrix");
    /* calculate eigenstuff */
    LaEigSolve(matrix, eigenvalues, eigenvectors);
    /* multiply them back together */
    recombine_diagonalised_matrices(matrix_size, eigenvectors, eigenvalues, result);
    print_matrix(result, "final matrix");
}
void test_matrix_exponential(){
    int matrix_size = 5;
    /* initialise everything */
    LaGenMatComplex matrix = LaGenMatComplex::rand(matrix_size, matrix_size, 0, 9);
    LaGenMatComplex result = LaGenMatComplex::zeros(matrix_size, matrix_size);
    print_matrix(matrix, "initial matrix");
    /* calculate exponential */
    matrix_exponential(matrix, matrix_size, result);
    print_matrix(result, "e^(matrix)");
}
void test_matrix_negative(){
    int matrix_size = 3;
    LaGenMatComplex matrix = LaGenMatComplex::rand(matrix_size, matrix_size, 0, 9);
    LaGenMatComplex result = LaGenMatComplex::zeros(matrix_size, matrix_size);
    print_matrix(matrix, "Matrix");
    matrix_negative(matrix_size, matrix, result);
    print_matrix(result, "- Matrix");
    matrix_negative(matrix_size, matrix);
    print_matrix(matrix, "- Matrix (in place)");
}
void test_matrix_equals_(){
    int matrix_size = 5;
    LaGenMatComplex matrix = LaGenMatComplex::rand(matrix_size, matrix_size, 0, 5);
    print_matrix(matrix, "initial matrix");
    matrix = 0;
    print_matrix(matrix, "matrix = 0");
    matrix = 1;
    print_matrix(matrix, "matrix = 1");
}
void test_diagonal_exponential(){
    /* initialise everything */
    int matrix_size = 3;
    LaGenMatComplex test = LaGenMatComplex::zeros(matrix_size, matrix_size);
    LaGenMatComplex result = LaGenMatComplex::zeros(matrix_size, matrix_size);
    for(int i = 0; i < matrix_size; i++){
        test(i,i).r = i;
    }
    /* calculate exponential */
    diagonal_matrix_exponential(test, matrix_size, result);
    print_matrix(test, "test");
    print_matrix(result);
}
// - qmc
void test_initial_parameters(){
    double U = 1, beta = 10, lambda, delta_tau;
    int lattice_size = 5, time_size;
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, time_size, lattice_size);
}
void test_generate_lattice(){
    int lattice_size = 5, time_size = 17;
    LaGenMatComplex lattice;
    generate_lattice(lattice_size, time_size, lattice);
    print_matrix(lattice, "lattice");
}
void test_H(){
    /* initialise everything */
    int lattice_size = 5;
    LaGenMatComplex H;
    LaVectorComplex eigenvalues = LaVectorComplex(lattice_size);
    LaGenMatComplex eigenvectors = LaGenMatComplex::zeros(lattice_size, lattice_size);
    /* generate matrix */
    generate_H(lattice_size, H);
    print_matrix(H);
    /* calculate eigenstuff */
    LaEigSolve(H, eigenvalues, eigenvectors);
    print_vector(eigenvalues, "eigenvalues");
    // eigenvalues are 2 cos(n pi / q), where q = the matrix size
}
void test_V(){
    /* initialise everything */
    int lattice_size = 5, time_size;
    LaGenMatComplex V = LaGenMatComplex::zeros(lattice_size, lattice_size);
    COMPLEX slice[lattice_size];
    double U = 1, beta = 10, lambda, delta_tau;

    /* calculate initial parameters */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, time_size, lattice_size);
    cout << endl;

    /* generate the lattice */
    generate_slice(lattice_size, slice);
    print_array(slice, lattice_size, "slice");

    /* calculate V */
    V_calculation(slice, lattice_size, U, lambda, 1, delta_tau, V);
    print_matrix(V, "V");
}


						/* ------ TO TEST ------ */
//...
void print_scalar(const COMPLEX scalar){
    cout << scalar << endl;
}
void print_scalar(const COMPLEX scalar, const string name){
    cout << name << ": " << scalar << endl;
}
void matrix_exponential_v(const LaGenMatComplex& matrix, const int matrix_size, LaGenMatComplex& result){
    /* initialise everything */
    LaVectorComplex eigenvalues = LaVectorComplex(matrix_size);
    LaGenMatComplex eigenvectors = LaGenMatComplex::zeros(matrix_size, matrix_size);
    LaGenMatComplex diagonalEigenExp = LaGenMatComplex::zeros(matrix_size, matrix_size);
    LaVectorComplex eigenExponential = LaVectorComplex(matrix_size);
    /* calculate eigenstuff */
    LaEigSolve(matrix, eigenvalues, eigenvectors);
    print_matrix(matrix, "initial matrix");
    print_vector(eigenvalues, "eigenvalues");
    print_matrix(eigenvectors, "eigenvectors");
        // matrix       = Input LaGenMatComplex
        // eigenvalues  = LaVectorComplex
        // eigenvectors = LaGenMatComplex
    /* calculate exponentials */
    for(int i = 0; i < matrix_size; i++){
        scalar_exponential(eigenvalues(i), eigenExponential(i));
        print_scalar(eigenExponential(i), "exp_ii");
    }
    cout << endl;
    print_matrix(eigenExponential, "exponential eigenvalues");
    print_matrix(result, "exponential eigenvalues - r");
    /* multiply them back together */
    recombine_diagonalised_matrices(matrix_size, eigenvectors, eigenExponential, result);
}
void test_negH_exponential(){
    /* initialise everything */
    int lattice_size = 5;
    LaGenMatComplex H = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex negH = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex expH = LaGenMatComplex::zeros(lattice_size, lattice_size);
    /* calculate H */
    generate_H(lattice_size, H);
    print_matrix(H, "H");
    /* calculate -H */
    matrix_negative(lattice_size, H, negH);
    print_matrix(negH, "-H");
    /* calculate exponentials */
    matrix_exponential_v(H, lattice_size, expH);
}
void B_calculation(const COMPLEX slice[], const int lattice_size, const double U, const double lambda, const double sigma, const double delta_tau, LaGenMatComplex& B){
    /* initialise everything */
    LaGenMatComplex H = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex V = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex negH = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex negV = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex expH = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex expV = LaGenMatComplex::zeros(lattice_size, lattice_size);

    /* calculate H and V */
    generate_H(lattice_size, H);
    V_calculation(slice, lattice_size, U, lambda, sigma, delta_tau, V);
    print_matrix(H, "H");
    print_matrix(V, "V");

    /* calculate -H and -V */
    matrix_negative(lattice_size, H, negH);
    matrix_negative(lattice_size, V, negV);
    print_matrix(negH, "-H");
    print_matrix(negV, "-V");

    /* calculate exponentials */
    matrix_exponential_v(negH, lattice_size, expH);
    diagonal_matrix_exponential(negV, lattice_size, expV);
    print_matrix(expH, "e^(-H)");
    print_matrix(expV, "e^(-V)");

    /* multiply exponentials */
    B = expH.copy();
    matrix_product(B, expV);
}
void test_B_generation(){
    /* initialise everything */
    int lattice_size = 5, time_size;
    double U = 1, beta = 10, lambda, sigma = 0, delta_tau;
    LaGenMatComplex B = LaGenMatComplex::zeros(lattice_size, lattice_size);
    COMPLEX slice[lattice_size];
    /* generate initial conditions */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    /* generate time slice */
    generate_slice(lattice_size, slice);
    print_initial_parameters(U, beta, lambda, delta_tau, time_size, lattice_size);
    /* calculate B */
    B_calculation(slice, lattice_size, U, lambda, sigma, delta_tau, B);
    /* print result */
    print_matrix(B,"B = e^-H e^-V");
}

						/* ------ TO CONVERT ------ */
// matrix_size                      -> lattice_size or time_size
// len                              -> array_size
// matrix_eigenvstuff               -> LaEigSolve
// float                            -> double
// scalar_exponential_main(n,i,r)   -> scalar_exponential(n,r)
// matrix_exponential(m,ms,i,r)     -> matrix_exponential(m,ms,r)
// generate_matrix                  -> LaGenMatComplex::rand
// scalar_product_f                 -> scalar_product
// basic_random_int                 -> random_int
// generate_lattice_array           -> generate_slice
// iterations                       -> ...
// void test_...(...)               -> void test_...( );
// five_matrix_multiplication       -> n_matrix_product

void store_matrix(const LaGenMatComplex& matrix, const int matrix_number, const int matrix_size, COMPLEX storage[], const int storage_size){
    /* initialise everything */
    int matrix_volume = matrix_size * matrix_size;
    /* store the matrix in storage */
    for(int r = 0; r < matrix_size; r++){
        for(int c = 0; c < matrix_size; c++){
            int i = (matrix_number * matrix_volume) + (r * matrix_size) + c;
            storage[i].r = matrix(r,c).r;
            storage[i].i = matrix(r,c).i;
        }
    }
}
void test_store_matrix(){
    /* initialise everything */
    int matrix_size = 5, storage_size = matrix_size * matrix_size * 3;
    COMPLEX storage[storage_size];
    LaGenMatComplex A = LaGenMatComplex::rand(matrix_size, matrix_size, 0, 9);
    LaGenMatComplex B = LaGenMatComplex::rand(matrix_size, matrix_size, 0, 9);
    LaGenMatComplex C = LaGenMatComplex::rand(matrix_size, matrix_size, 0, 9);
    /* store the matrices in storage */
    store_matrix(A, 0, matrix_size, storage, storage_size);
    print_array(storage, storage_size, "storage");
    store_matrix(B, 1, matrix_size, storage, storage_size);
    print_array(storage, storage_size, "storage");
    store_matrix(C, 2, matrix_size, storage, storage_size);
    print_array(storage, storage_size, "storage");
}
// void n_matrix_product(const COMPLEX storage[], const int matrix_size, const int n, LaGenMatComplex& result){
//     /* initialise everything */
//     LaGenMatComplex matrix;
//     int matrix_volume = matrix_size * matrix_size;
//     /* reset variables */
//     result = LaGenMatComplex::eye(matrix_size, matrix_size);
//     //for each matrix
//     for(int m = 0; m < n; m++){
//         // reset variables
//         matrix = LaGenMatComplex::eye(matrix_size, matrix_size);
//         // convert the storage to a matrix
//         for(int r = 0; r < matrix_size; r++){
//             for(int c = 0; c < matrix_size; c++){
//                 int e = r * matrix_size + c;
//                 int i = m * matrix_volume + e;
//                 matrix(r, c).r = storage[i].r;
//                 matrix(r, c).i = storage[i].i;
//             }
//         }
//         // multiply with the result
//         matrix_product(result, matrix);
//     }
// }
// void test_n_matrix_product(){
//
//     /* initialise everything */
//     int n = 4, matrix_size = 3, max_rand = 5;
//     int storage_size = matrix_size * matrix_size * n;
//     COMPLEX storage[storage_size];
//     LaGenMatComplex result = LaGenMatComplex::eye(matrix_size, matrix_size);
//
//
//     /* generate matrices (skip to storage) */
//     generate_real_array(storage, storage_size, max_rand);
//     print_array(storage, storage_size, "storage");
//
//     /* multiply everything */
//     n_matrix_product(storage, matrix_size, n, result);
//
//     print_matrix(result, "result");
// }
// void O_calculation(const int matrix_size, const LaGenMatComplex& BA, const LaGenMatComplex& BB, const LaGenMatComplex& BC, const LaGenMatComplex& BD, const LaGenMatComplex&BE, LaGenMatComplex& O){
//     //O = 1 + B(m) B(m-1) B(...) B(1)
//     /* initialise everything */
//     LaGenMatComplex I = LaGenMatComplex::eye(matrix_size, matrix_size);
//     //LaGenMatComplex multiplication;
//     /* multiply exponentials */
//     n_matrix_product(BA, BB, BC, BD, BE, O);
//     /* add I */
//     matrix_sum(matrix_size, O, I);
// }
// void test_O(){
//     /* initialise everything */
//     int time_size = 17;
//     COMPLEX elements[time_size];
//     LaGenMatComplex H;
//     LaGenMatComplex V = LaGenMatComplex::zeros(time_size, time_size);
//     LaGenMatComplex BA = LaGenMatComplex::zeros(time_size, time_size);
//     LaGenMatComplex BB = LaGenMatComplex::zeros(time_size, time_size);
//     LaGenMatComplex BC = LaGenMatComplex::zeros(time_size, time_size);
//     LaGenMatComplex BD = LaGenMatComplex::zeros(time_size, time_size);
//     LaGenMatComplex BE = LaGenMatComplex::zeros(time_size, time_size);
//     LaGenMatComplex O = LaGenMatComplex::zeros(time_size, time_size);
//     float U = 1, lambda = lambda_calculation(U), delta_tau = delta_tau_calculation(U);
//
//     /* generate matrices */
//     generate_H(time_size, H);
//     for(int i = 0; i < time_size; i++){
//         /* generate matrices */
//         generate_slice(time_size, elements);
//         V_calculation(elements, time_size, U, lambda, 1, delta_tau, V);
//         /* calculate B */
//         if(i == 0){
//             B_calculation(slice, lattice_size, U, lambda, 1, delta_tau, BA);
//         }else if(i == 1){
//             B_calculation(slice, lattice_size, U, lambda, 1, delta_tau, BB);
//         }else if(i == 2){
//             B_calculation(slice, lattice_size, U, lambda, 1, delta_tau, BC);
//         }else if(i == 3){
//             B_calculation(slice, lattice_size, U, lambda, 1, delta_tau, BD);
//         }else if(i == 4){
//             B_calculation(slice, lattice_size, U, lambda, 1, delta_tau, BE);
//         }
//     }
//     O_calculation(time_size, BA, BB, BC, BD, BE, O);
//     /* print result */
//     print_matrix(BA, "BA");
//     print_matrix(BB, "BB");
//     print_matrix(BC, "BC");
//     print_matrix(BD, "BD");
//     print_matrix(BE, "BE");
//     print_matrix(O, "O");
// }

/* ------ Main QMC Program ------ */
int main(){
    test_negH_exponential();
}
