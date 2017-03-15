#include <iostream> //cout
#include <fstream>
#include <string>
#include "improv.h"
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
void print_scalar(const COMPLEX scalar){
    cout << scalar << endl;
}
void print_scalar(const COMPLEX scalar, const string name){
    cout << name << ": " << scalar << endl;
}
void print_scalar(const double scalar, const string name){
    cout << name << ": " << scalar << endl;
}
void print_scalar(const COMPLEX scalar, const string name, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* print the scalar */
    myfile << name << ": " << scalar << endl;
    /* close the file */
    myfile.close();
}
void print_scalar(const double scalar, const string name, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* print the scalar */
    myfile << name << ": " << scalar << endl;
    /* close the file */
    myfile.close();
}
void print_array(const COMPLEX array[], int array_size, const string name){
	cout << name << ": ";
    for(int i = 0; i < array_size; i++){
        cout << array[i] << " ";
    }
    cout << endl;
}
void print_array(const COMPLEX array[], int array_size, const string name, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* print stuff */
	myfile << name << ": ";
    for(int i = 0; i < array_size; i++){
        myfile << array[i] << " ";
    }
    myfile << endl;
    /* close the file */
    myfile.close();
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
void print_matrix(const LaGenMatComplex& matrix, const string name, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);

    /* print the matrix */
	myfile << name << ":" << endl << matrix << endl;
    /* close the file */
    myfile.close();
}
void print_initial_parameters(double U, double beta, double lambda, double delta_tau, double mu, int time_size, int lattice_size){
	cout << "no of lattice points = " << lattice_size << endl;
	cout << "no of time slices = " << time_size << endl;
	cout << "U = " << U << endl;
	cout << "beta = " << beta << endl;
	cout << "lambda = " << lambda << endl;
	cout << "delta tau = " << delta_tau << endl;
    cout << "mu = " << mu << endl;
}
void print_initial_parameters(double U, double beta, double lambda, double delta_tau, double mu, int time_size, int lattice_size, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* print the initial parameters */
	myfile << "no of lattice points = " << lattice_size << endl;
	myfile << "no of time slices = " << time_size << endl;
	myfile << "U = " << U << endl;
	myfile << "beta = " << beta << endl;
	myfile << "lambda = " << lambda << endl;
	myfile << "delta tau = " << delta_tau << endl;
    myfile << "mu = " << mu << endl;
    myfile << endl;
    /* close the file */
    myfile.close();
}

/* -- Processing -- */
// Randomisation
int random_int(const int max_rand){
    random_device rd;
    mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, max_rand);
    return dist(gen);
}
double random_double(){
    random_device rd;
    mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    return dis(gen);
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
void clear_scalar(COMPLEX& scalar){
    scalar.r = 0;
    scalar.i = 0;
}
void clear_array(COMPLEX array[], const int array_size){
    for(int i = 0; i < array_size; i++){
        array[i].r = 0;
        array[i].i = 0;
    }
}
void isolate_row(const LaGenMatComplex& matrix, const int matrix_width, const int row, COMPLEX array[]){
    for(int i = 0; i < matrix_width; i++){
        array[i] = matrix(row, i);
    }
}
void flip_scalar(COMPLEX& spin){
    spin.r = -spin.r;
    spin.i = -spin.i;
}
void flip_spin(LaGenMatComplex& lattice, const int t, const int l){
    lattice(t,l).r = -lattice(t,l).r;
    lattice(t,l).i = -lattice(t,l).i;
}
void flip_spin_v(LaGenMatComplex& lattice, const int t, const int l){
    cout << "flipped ("<<t<<", "<<l<<"): " << lattice(t,l);
    lattice(t,l).r = -lattice(t,l).r;
    lattice(t,l).i = -lattice(t,l).i;
    cout << " -> " << lattice(t,l) << endl;
}
void flip_spin(LaGenMatComplex& lattice, const int t, const int l, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* print stuff */
    myfile << "flipped ("<<t<<", "<<l<<"): " << lattice(t,l);
    lattice(t,l).r = -lattice(t,l).r;
    lattice(t,l).i = -lattice(t,l).i;
    myfile << " -> " << lattice(t,l) << endl;
    /* close the file */
    myfile.close();
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
double check_size(const double number){
    return floor(log10(number));
}
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
void scalar_multiplication(const COMPLEX& A, const COMPLEX& B, COMPLEX& result){
    la::complex<double> laA = la::complex<double>(A);
    la::complex<double> laB = la::complex<double>(B);
    la::complex<double> laResult = la::complex<double>(result);
    laResult = laA * laB;
    result = laResult.toCOMPLEX();
}
COMPLEX scalar_multiple(COMPLEX& A, const COMPLEX& B){
    COMPLEX result;
    result.r = (A.r * B.r) - (A.i * B.i);
    result.i = (A.r * B.i) + (A.i * B.r);
    return result;
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
void matrix_sum(const int matrix_size, LaGenMatComplex& sum, const LaGenMatComplex& matrix){
    for(int i = 0; i < matrix_size; i++){
        for(int j = 0; j < matrix_size; j++){
            sum(i, j).r += matrix(i, j).r;
            sum(i, j).i += matrix(i, j).i;
        }
    }
}
void matrix_multiple(const LaGenMatComplex& matrix, const int matrix_size, const double number, LaGenMatComplex& result){
    for(int r = 0; r < matrix_size; r++){
        for(int c = 0; c < matrix_size; c++){
            result(r,c).r = matrix(r,c).r * number;
            result(r,c).i = matrix(r,c).i * number;
        }
    }
}
void matrix_product(LaGenMatComplex& product, const LaGenMatComplex& matrix){
    LaGenMatComplex result = matrix.copy();
    Blas_Mat_Mat_Mult(product, matrix, result);
    product = result.copy();
}
void matrix_inverse(const LaGenMatComplex& matrix, int matrix_size, LaGenMatComplex& result){
    result = matrix.copy();
    LaVectorLongInt PIV = LaVectorLongInt(matrix_size);
    LUFactorizeIP(result, PIV);
    LaLUInverseIP(result, PIV);
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
        scalar_exponential(eigenvalues(i), eigenExponential(i));
    }
    /* multiply them back together */
    recombine_diagonalised_matrices(matrix_size, eigenvectors, eigenExponential, result);
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
    /* calculate exponentials */
    for(int i = 0; i < matrix_size; i++){
        scalar_exponential(eigenvalues(i), eigenExponential(i));
        print_scalar(eigenExponential(i), "exp_ii");
    }
    cout << endl;
    print_matrix(eigenExponential, "exponential eigenvalues");
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
void matrix_determinant_e(const int matrix_size, const LaGenMatComplex& matrix, COMPLEX& result){
    /* initialise everything */
    LaVectorComplex eigenvalues = LaVectorComplex(matrix_size);
    LaGenMatComplex eigenvectors = LaGenMatComplex::zeros(matrix_size, matrix_size);
    result.r = 1;
    result.i = 0;
    /* calculate eigenvectors */
    LaEigSolve(matrix, eigenvalues, eigenvectors);
    /* calculate determinant */
    for(int i = 0; i < matrix_size; i++){
        scalar_product(result, eigenvalues(i));
    }
}
COMPLEX simple_matrix_determinant(const LaGenMatComplex& matrix){
    /* initialise everything */
    COMPLEX AD;
    COMPLEX BC;
    COMPLEX det;
    int max_size = 1000000;
    /* multiply opposite corners */
    scalar_multiplication(matrix(0,0), matrix(1,1), AD);
    scalar_multiplication(matrix(0,1), matrix(1,0), BC);
    /* - B */
    det.r = AD.r - BC.r;
    det.i = AD.i - BC.i;
    return det;
}
COMPLEX determinant_coefficient(const LaGenMatComplex& matrix, const int i){
    COMPLEX coefficient;
    if(i % 2 == 1){
        coefficient.r = - matrix(0, i).r;
        coefficient.i = - matrix(0, i).i;
    }else{
        coefficient.r = matrix(0, i).r;
        coefficient.i = matrix(0, i).i;
    }
    return coefficient;
}
void generate_cofactor_matrix(const int matrix_size, const LaGenMatComplex& matrix, const int i, LaGenMatComplex& cofactorMatrix){
    for(int r = 1; r < matrix_size; r++){
        int newC = 0;
        for(int c = 0; c < matrix_size; c++){
            if(c != i){
                cofactorMatrix(r - 1, newC).r = matrix(r, c).r;
                cofactorMatrix(r - 1, newC).i = matrix(r, c).i;
                newC++;
            }
        }
    }
}
COMPLEX matrix_determinant(const int matrix_size, const LaGenMatComplex& matrix){
    /* initialise everything */
    COMPLEX determinant;
    LaGenMatComplex cofactorMatrix;
    COMPLEX coefficient;
    cofactorMatrix = 0;
    /* calculate determinant */
    if(matrix_size == 2){
        return simple_matrix_determinant(matrix);
    }else{
        clear_scalar(determinant);
        clear_scalar(coefficient);
        //for each i in the first row
        for(int i = 0; i < matrix_size; i++){
            /* initialise everything */
            int cofactor_size = matrix_size - 1;
            cofactorMatrix = LaGenMatComplex::zeros(cofactor_size, cofactor_size);
            /* determine the coefficient */
            coefficient = determinant_coefficient(matrix, i);
            /* calculate the cofactor */
            generate_cofactor_matrix(matrix_size, matrix, i, cofactorMatrix);
            /* finish calculation */
            scalar_sum(determinant, scalar_multiple(coefficient, matrix_determinant(cofactor_size, cofactorMatrix)));
        }
        return determinant;
    }
}
COMPLEX matrix_determinant_v(const int matrix_size, const LaGenMatComplex& matrix){
    /* initialise everything */
    LaGenMatComplex scaled_matrix = LaGenMatComplex::zeros(matrix_size, matrix_size);
    LaGenMatComplex cofactorMatrix;
    COMPLEX determinant;
    COMPLEX coefficient;
    cofactorMatrix = 0;
    /* test size of elements */
    double scale = check_size(matrix(0,0).r);
    /* scale matrix */
    matrix_multiple(matrix, matrix_size, 1 / scale, scaled_matrix);
    print_matrix(scaled_matrix, "scaled matrix");
    /* do stuff */
    if(matrix_size == 2){
        return simple_matrix_determinant(matrix);
    }else{
        clear_scalar(determinant);
        clear_scalar(coefficient);
        //for each i in the first row
        for(int i = 0; i < matrix_size; i++){
            /* initialise everything */
            int cofactor_size = matrix_size - 1;
            cofactorMatrix = LaGenMatComplex::zeros(cofactor_size, cofactor_size);
            /* determine the coefficient */
            coefficient = determinant_coefficient(matrix, i);
            print_scalar(coefficient, "coefficient");
            /* calculate the cofactor */
            generate_cofactor_matrix(matrix_size, matrix, i, cofactorMatrix);
            print_matrix(cofactorMatrix, "cofactorMatrix");
            /* finish calculation */
            scalar_sum(determinant, scalar_multiple(coefficient, matrix_determinant(cofactor_size, cofactorMatrix)));
        }
        cout << pow (10, scale * matrix_size) << endl;
        scalar_product(determinant, pow (10, scale * matrix_size));
        return determinant;
    }
}
// - qmc
void initial_parameter_calculation(const double U, const double beta, double& lambda, double& delta_tau, int& time_size){
    // mu = U / 2;                             // by definition
    lambda = acoshf(exp(sqrt(0.125*U)/2));  // by definition
    time_size = ceil(beta / lambda);        // by definition
    delta_tau = beta / time_size;           // = sqrt(0.125 / U) by convension
}
void H_generation(const int lattice_size, LaGenMatComplex& H){
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
void V_calculation(const COMPLEX slice[], const int lattice_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& V){
    /* initialise everything */
    V = 0;
    COMPLEX V_ii[lattice_size];

    /* V_ii = (lambda sigma s_l / delta_tau) + mu - U / 2 */
    for(int i = 0; i < lattice_size; i++){
        V_ii[i].r = lambda * sigma * slice[i].r / delta_tau + mu - U / 2;
        V_ii[i].i = 0;
    }

    /* plot to diagonal */
    array_to_diag(V_ii, lattice_size, V);
}
void B_calculation(const COMPLEX slice[], const int lattice_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& B){
    /* initialise everything */
    B = 0;
    LaGenMatComplex H = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex V = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex sum = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex product = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex negative = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex exponential = LaGenMatComplex::zeros(lattice_size, lattice_size);

    /* calculate H and V */
    H_generation(lattice_size, H);
    V_calculation(slice, lattice_size, U, lambda, sigma, delta_tau, mu, V);
    sum = H.copy();

    /* calculate H + V */
    matrix_sum(lattice_size, sum, V);

    /* calculate delta_tau * (H + V) */
    matrix_multiple(sum, lattice_size, delta_tau, product);

    /* calculate - delta_tau * (H + V) */
    matrix_negative(lattice_size, product, negative);

    /* calculate exp(- delta_tau * (H + V)) */
    matrix_exponential(negative, lattice_size, B);
}
void B_calculation_v(const COMPLEX slice[], const int lattice_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& B){
    /* initialise everything */
    B = 0;
    LaGenMatComplex H = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex V = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex sum = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex product = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex negative = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex exponential = LaGenMatComplex::zeros(lattice_size, lattice_size);

    /* calculate H and V */
    H_generation(lattice_size, H);
    V_calculation(slice, lattice_size, U, lambda, sigma, delta_tau, mu, V);
    sum = H.copy();
    print_matrix(sum, "H");
    print_matrix(V, "V");

    /* calculate H + V */
    matrix_sum(lattice_size, sum, V);
    print_matrix(sum, "H + V");

    /* calculate delta_tau * (H + V) */
    matrix_multiple(sum, lattice_size, delta_tau, product);
    print_matrix(product, "delta_tau * (H + V)");

    /* calculate - delta_tau * (H + V) */
    matrix_negative(lattice_size, product, negative);
    print_matrix(negative, "- delta_tau * (H + V)");

    /* calculate exp(- delta_tau * (H + V)) */
    matrix_exponential(negative, lattice_size, B);
    print_matrix(B, "B = exp(- delta_tau * (H + V))");
}
void B_calculation(const COMPLEX slice[], const int lattice_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& B, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* initialise everything */
    B = 0;
    LaGenMatComplex H = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex V = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex sum = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex product = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex negative = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex exponential = LaGenMatComplex::zeros(lattice_size, lattice_size);

    /* calculate H and V */
    H_generation(lattice_size, H);
    V_calculation(slice, lattice_size, U, lambda, sigma, delta_tau, mu, V);
    print_matrix(H, "H", file);
    print_matrix(V, "V", file);

    /* calculate H + V */
    sum = H.copy();
    matrix_sum(lattice_size, sum, V);
    print_matrix(sum, "H + V", file);

    /* calculate delta_tau * (H + V) */
    matrix_multiple(sum, lattice_size, delta_tau, product);
    print_matrix(product, "delta_tau * (H + V)", file);

    /* calculate - delta_tau * (H + V) */
    matrix_negative(lattice_size, product, negative);
    print_matrix(negative, "- delta_tau * (H + V)", file);

    /* calculate exp(- delta_tau * (H + V)) */
    matrix_exponential(negative, lattice_size, B);
    print_matrix(B, "B = exp(- delta_tau * (H + V))", file);
    /* close the file */
    myfile.close();
}
void O_calculation(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& O){
    /* initialise everything */
    LaGenMatComplex B = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex I = LaGenMatComplex::eye(lattice_size, lattice_size);
    COMPLEX slice[lattice_size];
    O = LaGenMatComplex::eye(lattice_size, lattice_size);
    /* calculate B matrices */
    for(int x = 0; x < time_size; x++){
        clear_array(slice, lattice_size);
        int t = time_size - x - 1;
        isolate_row(lattice, lattice_size, t, slice);
        B_calculation(slice, lattice_size, U, lambda, sigma, delta_tau, mu, B);
        matrix_product(O, B);
    }
    /* add I */
    matrix_sum(lattice_size, O, I);
}
void O_calculation_v(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& O){
    /* initialise everything */
    LaGenMatComplex B = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex I = LaGenMatComplex::eye(lattice_size, lattice_size);
    COMPLEX slice[lattice_size];
    O = LaGenMatComplex::eye(lattice_size, lattice_size);
    print_matrix(O, "product");
    /* calculate B matrices */
    for(int x = 0; x < time_size; x++){
        clear_array(slice, lattice_size);
        int t = time_size - x - 1;
        cout << "t = " << t << ": ";
        isolate_row(lattice, lattice_size, t, slice);
        // print_array(slice, lattice_size, "slice");
        B_calculation_v(slice, lattice_size, U, lambda, sigma, delta_tau, mu, B);
        // print_matrix(B, "B");
        matrix_product(O, B);
        print_matrix(O, "product");
    }
    /* add I */
    matrix_sum(lattice_size, O, I);
}
void O_calculation(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& O, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* initialise everything */
    LaGenMatComplex B = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex I = LaGenMatComplex::eye(lattice_size, lattice_size);
    COMPLEX slice[lattice_size];
    O = LaGenMatComplex::eye(lattice_size, lattice_size);
    print_matrix(O, "initial product", file);
    /* calculate B matrices */
    for(int t = time_size - 1; t >= 0 ; t--){
        myfile << "t = " << t << ": " << endl;
        /* isolate the time slice */
        clear_array(slice, lattice_size);
        isolate_row(lattice, lattice_size, t, slice);
        print_array(slice, lattice_size, "slice", file);
        /* calculate the B matrix */
        B_calculation(slice, lattice_size, U, lambda, sigma, delta_tau, mu, B);
        print_matrix(B, "B", file);
        /* calculate the O matrix */
        matrix_product(O, B);
        print_matrix(O, "product", file);
    }
    /* add I */
    matrix_sum(lattice_size, O, I);
    /* close the file */
    myfile.close();
}
void weight_calculation(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, const double mu, COMPLEX& weight){
    /* initialise everything */
    LaGenMatComplex OUP = LaGenMatComplex::zeros(lattice_size,lattice_size);
    LaGenMatComplex ODN = LaGenMatComplex::zeros(lattice_size,lattice_size);
    COMPLEX detOUP;
    COMPLEX detODN;
    clear_scalar(weight);
    clear_scalar(detOUP);
    clear_scalar(detODN);
    /* calculate O */
    O_calculation(lattice, lattice_size, time_size, U, lambda, 1, delta_tau, mu, OUP);
    O_calculation(lattice, lattice_size, time_size, U, lambda, -1, delta_tau, mu, ODN);
    /* calculate det(O) */
    matrix_determinant_e(lattice_size, OUP, detOUP);
    matrix_determinant_e(lattice_size, ODN, detODN);
    /* calculate weight */
    weight = scalar_multiple(detOUP, detODN);
}
void weight_calculation_v(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, const double mu, COMPLEX& weight){
    /* initialise everything */
    LaGenMatComplex OUP = LaGenMatComplex::zeros(lattice_size,lattice_size);
    LaGenMatComplex ODN = LaGenMatComplex::zeros(lattice_size,lattice_size);
    COMPLEX detOUP;
    COMPLEX detODN;
    clear_scalar(weight);
    clear_scalar(detOUP);
    clear_scalar(detODN);
    /* calculate O */
    cout << "sigma = 1" << endl;
    O_calculation_v(lattice, lattice_size, time_size, U, lambda, 1, delta_tau, mu, OUP);
    cout << "sigma = -1" << endl;
    O_calculation_v(lattice, lattice_size, time_size, U, lambda, -1, delta_tau, mu, ODN);
    print_matrix(OUP, "O UP");
    print_matrix(ODN, "O DN");
    /* calculate det(O) */
    matrix_determinant_e(lattice_size, OUP, detOUP);
    matrix_determinant_e(lattice_size, ODN, detODN);
    print_scalar(detOUP, "det(O UP)");
    print_scalar(detODN, "det(O DN)");
    /* calculate weight */
    weight = scalar_multiple(detOUP, detODN);
    print_scalar(weight, "weight");
}
void weight_calculation_vv(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, const double mu, COMPLEX& weight, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* initialise everything */
    LaGenMatComplex OUP = LaGenMatComplex::zeros(lattice_size,lattice_size);
    LaGenMatComplex ODN = LaGenMatComplex::zeros(lattice_size,lattice_size);
    COMPLEX detOUP;
    COMPLEX detODN;
    COMPLEX detOUP_e;
    COMPLEX detODN_e;
    clear_scalar(weight);
    clear_scalar(detOUP);
    clear_scalar(detODN);
    clear_scalar(detOUP_e);
    clear_scalar(detODN_e);

    /* calculate O */
    myfile << "sigma = 1" << endl;
    O_calculation(lattice, lattice_size, time_size, U, lambda,  1, delta_tau, mu, OUP, file);
    print_matrix(OUP, "O UP", file);
    myfile << "sigma = -1" << endl;
    O_calculation(lattice, lattice_size, time_size, U, lambda, -1, delta_tau, mu, ODN, file);
    print_matrix(ODN, "O DN", file);

    /* calculate det(O) */
    myfile << "my determinants:" << endl;
    detOUP = matrix_determinant(lattice_size, OUP);
    detODN = matrix_determinant(lattice_size, ODN);
    print_scalar(detOUP, "det(O UP)", file);
    print_scalar(detODN, "det(O DN)", file);
    myfile << "eigen determinants:" << endl;
    matrix_determinant_e(lattice_size, OUP, detOUP_e);
    matrix_determinant_e(lattice_size, ODN, detODN_e);
    print_scalar(detOUP_e, "det(O UP) (e)", file);
    print_scalar(detODN_e, "det(O DN) (e)", file);
    myfile << endl;

    // figure out what to do about O being too big !!!
        // det(O) should be ~ 10^44

    /* calculate weight */
    weight = scalar_multiple(detOUP, detODN);
    print_scalar(weight, "weight", file);
    if(weight.i != 0){
        print_scalar(weight, "weight");
        print_matrix(lattice, "lattice", "imag_weight.txt");
        print_matrix(OUP, "O UP", "imag_weight.txt");
        print_matrix(ODN, "O DN", "imag_weight.txt");
        print_scalar(detOUP, "det(O UP)", "imag_weight.txt");
        print_scalar(detODN, "det(O DN)", "imag_weight.txt");
        print_scalar(weight, "weight", "imag_weight.txt");
        myfile << endl;
    }
    myfile << endl;
    /* close the file */
    myfile.close();
}
void weight_calculation_v(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, const double mu, COMPLEX& weight, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* initialise everything */
    LaGenMatComplex OUP = LaGenMatComplex::zeros(lattice_size,lattice_size);
    LaGenMatComplex ODN = LaGenMatComplex::zeros(lattice_size,lattice_size);
    COMPLEX detOUP;
    COMPLEX detODN;
    clear_scalar(weight);
    clear_scalar(detOUP);
    clear_scalar(detODN);
    /* calculate O */
    // myfile << "sigma = 1" << endl;
    O_calculation(lattice, lattice_size, time_size, U, lambda,  1, delta_tau, mu, OUP);
    print_matrix(OUP, "O UP", file);
    // myfile << "sigma = -1" << endl;
    O_calculation(lattice, lattice_size, time_size, U, lambda, -1, delta_tau, mu, ODN);
    print_matrix(ODN, "O DN", file);

    /* calculate det(O) */
    detOUP = matrix_determinant(lattice_size, OUP);
    detODN = matrix_determinant(lattice_size, ODN);
    print_scalar(detOUP, "det(O UP)", file);
    print_scalar(detODN, "det(O DN)", file);
    myfile << endl;

    /* calculate weight */
    weight = scalar_multiple(detOUP, detODN);
    print_scalar(weight, "weight", file);
    if(weight.i != 0){
        print_scalar(weight, "weight");
        print_matrix(lattice, "lattice", "imag_weight.txt");
        print_matrix(OUP, "O UP", "imag_weight.txt");
        print_matrix(ODN, "O DN", "imag_weight.txt");
        print_scalar(detOUP, "det(O UP)", "imag_weight.txt");
        print_scalar(detODN, "det(O DN)", "imag_weight.txt");
        print_scalar(weight, "weight", "imag_weight.txt");
        myfile << endl;
    }
    myfile << endl;
    /* close the file */
    myfile.close();
}
void weight_calculation(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, const double mu, COMPLEX& weight, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* initialise everything */
    LaGenMatComplex OUP = LaGenMatComplex::zeros(lattice_size,lattice_size);
    LaGenMatComplex ODN = LaGenMatComplex::zeros(lattice_size,lattice_size);
    COMPLEX detOUP;
    COMPLEX detODN;
    clear_scalar(weight);
    clear_scalar(detOUP);
    clear_scalar(detODN);
    /* calculate O */
    O_calculation(lattice, lattice_size, time_size, U, lambda,  1, delta_tau, mu, OUP);
    O_calculation(lattice, lattice_size, time_size, U, lambda, -1, delta_tau, mu, ODN);
    /* calculate det(O) */
    detOUP = matrix_determinant(lattice_size, OUP);
    detODN = matrix_determinant(lattice_size, ODN);
    /* calculate weight */
    weight = scalar_multiple(detOUP, detODN);
    if(weight.i != 0){
        print_scalar(weight, "weight");
        print_matrix(lattice, "lattice", "imag_weight.txt");
        print_matrix(OUP, "O UP", "imag_weight.txt");
        print_matrix(ODN, "O DN", "imag_weight.txt");
        print_scalar(detOUP, "det(O UP)", "imag_weight.txt");
        print_scalar(detODN, "det(O DN)", "imag_weight.txt");
        print_scalar(weight, "weight", "imag_weight.txt");
        myfile << endl;
    }
    /* close the file */
    myfile.close();
}

void sweep_lattice(LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, const double mu, const int iterations, double& acceptance, double& rejection){

    /* initialise everything */
    COMPLEX weightBefore;
    COMPLEX weightAfter;
    clear_scalar(weightBefore);
    clear_scalar(weightAfter);
    double probability = 0;
    string result;
    acceptance = 0;
    rejection = 0;

    /* sweep through the lattice */
    for(int i = 0; i < iterations; i++){
        for(int t = 0; t < time_size; t++){
            for(int l = 0; l < lattice_size; l++){
                /* calculate the weight before the flip */
                weight_calculation(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, weightBefore);

                /* propose the flip */
                flip_spin(lattice, t, l);

                /* calculate the weight after the flip */
                weight_calculation(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, weightAfter);

                /* calculate the ratio of weights */
                probability = weightAfter.r / weightBefore.r;

                /* accept or reject the flip */
                double prob = random_double();
                if(abs(probability) >= 1){
                    acceptance++;
                }else{
                    if(probability > prob){
                        acceptance++;
                    }else{
                        flip_spin(lattice, t, l);
                        rejection++;
                    }
                }
            }
        }
    }
}
void sweep_lattice_v(LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, const double mu, const int iterations, double& acceptance, double& rejection){

    /* initialise everything */
    COMPLEX weightBefore;
    COMPLEX weightAfter;
    clear_scalar(weightBefore);
    clear_scalar(weightAfter);
    double probability = 0;
    string result;
    int count = 0;
    acceptance = 0;
    rejection = 0;

    /* output headings */
    cout.width(11);
    cout << "weight";
    cout << " lattice" << endl;

    /* sweep through the lattice */
    for(int i = 0; i < iterations; i++){
        for(int t = 0; t < time_size; t++){
            for(int l = 0; l < lattice_size; l++){
                /* calculate the weight before the flip */
                weight_calculation(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, weightBefore);

                /* propose the flip */
                flip_spin(lattice, t, l);

                /* calculate the weight after the flip */
                weight_calculation(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, weightAfter);

                /* calculate the ratio of weights */
                probability = weightAfter.r / weightBefore.r;

                /* accept or reject the flip */
                double prob = random_double();
                if(abs(probability) >= 1){
                    result = "accepted";
                    acceptance++;
                }else{
                    if(probability > prob){
                        result = "accepted";
                        acceptance++;
                    }else{
                        flip_spin(lattice, t, l);
                        result = "rejected";
                        rejection++;
                    }
                }
                /* comments */
                    //for negative values, we do some integration
                    //P\to\tilde{P} = |P| and  F\to \tilde
                    //you have to multiply each quan you measure bu the sign
                count++;
                if(count%1000 == 0){
                    cout << " (" << count <<") " << "[" << acceptance << "/" << rejection << "] " << result << " - probability: " << probability;
                    cout.width(15);
                    cout << " - weightBefore: " << weightBefore << ", weightAfter: " << weightAfter << endl;
                }
                // if(result == "accepted"){
                //     print_matrix(lattice);
                // }else{
                //     cout << endl;
                // }
            }
            /* Comments */
                //when you take measurements, there is noise
                //we're doing marcov chain
                //the simplest quan we measure is double occupancy \bra n_up n_down \ket
        }
    }
    //results
        // with most parameters = 1, it stabilised at all -1 spins
    cout << "["<< acceptance << "/" << rejection << "]" << endl;
    double acceptance_ratio = acceptance / (rejection + acceptance);
    cout << "acceptance ratio = " << acceptance_ratio << endl;
    double percentage_acceptance = acceptance / rejection;
    cout << "percentage acceptance = " << percentage_acceptance << endl << endl;
}
void calculate_total_spin_f(const LaGenMatComplex& lattice, const int time_size, const int lattice_size){
    /* open the file */
    ofstream myfile;
    myfile.open("spin.txt", std::ios_base::app);
    /* calculate total spin */
    double total_spin = 0;
    for(int t = 0; t < time_size; t++){
        for(int l = 0; l < lattice_size; l++){
            total_spin += lattice(t,l).r;
            myfile << total_spin << endl;
        }
    }
    myfile << endl;
    double average_spin = total_spin / (time_size * lattice_size);
    myfile << average_spin << endl << endl;
    /* close the file */
    myfile.close();
}
void calculate_total_spin(const LaGenMatComplex& lattice, const int time_size, const int lattice_size, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* calculate total spin */
    double total_spin = 0;
    for(int t = 0; t < time_size; t++){
        for(int l = 0; l < lattice_size; l++){
            total_spin += lattice(t,l).r;
        }
    }
    myfile << total_spin << endl;
    double average_spin = total_spin / (time_size * lattice_size);
    myfile << average_spin << endl << endl;
    /* close the file */
    myfile.close();
}
/* -- Testing -- */
// - generic
void test_output_to_file(const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* print stuff */
    myfile << "Here, have a test :) .\n";
    /* close the file */
    myfile.close();
}
void test_output_functions(const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);

    /* initialise everything */
    int array_size = 5, matrix_size = 3;

    COMPLEX scalar;
    scalar.r = 1;
    scalar.i = 2;

    COMPLEX array[array_size];
    generate_slice(array_size, array);

    LaGenMatComplex matrix = LaGenMatComplex::rand(matrix_size, matrix_size, 0, 4);

    /* print everything */
    print_scalar(scalar, "scalar", file);
    myfile << endl;
    print_array(array, array_size, "array", file);
    myfile << endl;
    print_matrix(matrix, "matrix", file);
    myfile << endl;

    /* close the file */
    myfile.close();
}
void test_flip_spins(){
    /* initialise stuff */
    int lattice_size = 5, time_size = 8;
    int l = random_int(lattice_size-1), t = random_int(time_size-1);
    LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
    /* generate lattice */
    generate_lattice(lattice_size, time_size, lattice);
    print_matrix(lattice, "lattice");
    /* flip spins */
    flip_spin_v(lattice, t, l);
}
void test_flip_spins(const string file){
    /* initialise stuff */
    int lattice_size = 5, time_size = 8;
    int l = random_int(lattice_size-1), t = random_int(time_size-1);
    LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
    /* generate lattice */
    generate_lattice(lattice_size, time_size, lattice);
    print_matrix(lattice, "lattice", file);
    /* flip spins */
    flip_spin(lattice, t, l, file);
}
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
void test_matrix_equals_(){
    int matrix_size = 5;
    LaGenMatComplex matrix = LaGenMatComplex::rand(matrix_size, matrix_size, 0, 5);
    print_matrix(matrix, "initial matrix");
    matrix = 0;
    print_matrix(matrix, "matrix = 0");
    matrix = 1;
    print_matrix(matrix, "matrix = 1");
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
void test_matrix_multiple(){
    int matrix_size = 4;
    LaGenMatComplex matrix = LaGenMatComplex::rand(matrix_size, matrix_size, 0, 9);
    LaGenMatComplex result = LaGenMatComplex::zeros(matrix_size, matrix_size);
    matrix_multiple(matrix, matrix_size, 2.0, result);
    print_matrix(matrix, "initial matrix");
    print_matrix(result, "initial matrix * 2");

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
void test_simple_matrix_determinant(){
    /* initialise everything */
    LaGenMatComplex matrix = LaGenMatComplex::rand(2,2,0,5);
    print_matrix(matrix, "matrix");
    /* calculate determinant */
    simple_matrix_determinant(matrix);
    print_scalar(simple_matrix_determinant(matrix), "det(M)");
}
void test_determinant_coefficient(){
    /* initialise everything */
    LaGenMatComplex matrix = LaGenMatComplex::rand(4,4,0,5);
    print_matrix(matrix, "matrix");
    /* calculate coefficients */
    for(int i = 0; i < 4; i++){
        cout << determinant_coefficient(matrix, i) << endl;
    }
    cout << endl;
}
void test_cofactor_matrix(){
    /* initialise everything */
    LaGenMatComplex matrix = LaGenMatComplex::rand(4,4,0,5);
    LaGenMatComplex cofactor = LaGenMatComplex::zeros(3,3);
    print_matrix(matrix, "matrix");
    /* calculate reduced matrix */
    for(int i = 0; i < 4; i++){
        generate_cofactor_matrix(4, matrix, i, cofactor);
        print_matrix(cofactor, "cofactor matrix");
    }
}
void test_matrix_determinant(){
    /* initialise everything */
    double matrix_size = 4, scale = 8;
    LaGenMatComplex matrix = LaGenMatComplex::rand(4,4,0,5);
    LaGenMatComplex scaled_matrix = LaGenMatComplex::zeros(matrix_size, matrix_size);
    COMPLEX result;
    clear_scalar(result);
    /* scale the matrix */
    matrix_multiple(matrix, matrix_size, scale, scaled_matrix);
    print_matrix(matrix, "initial matrix");
    print_matrix(scaled_matrix, "scaled matrix");
    /* calculate determinant */
    print_scalar(matrix_determinant_v(4, matrix), "my determinant");
    matrix_determinant_e(4, matrix, result);
    print_scalar(result, "eigenvalue determinant");

    print_scalar(matrix_determinant_v(4, scaled_matrix), "my determinant");
    matrix_determinant_e(4, scaled_matrix, result);
    print_scalar(result, "eigenvalue determinant");
}
// - qmc
void test_initial_parameters(){
    double U = 1, beta = 10, lambda, delta_tau, mu = U / 2;
    int lattice_size = 5, time_size;
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size);
}
void test_initial_parameters(const string file){
    double U = 1, beta = 10, lambda, delta_tau, mu = U / 2;
    int lattice_size = 5, time_size;
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, file);
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
    H_generation(lattice_size, H);
    print_matrix(H);
    /* calculate eigenstuff */
    LaEigSolve(H, eigenvalues, eigenvectors);
    print_vector(eigenvalues, "eigenvalues");
    // eigenvalues are 2 cos(n pi / q), where q = the matrix size
}
void test_negH_exponential(){
    /* initialise everything */
    int lattice_size = 5;
    LaGenMatComplex H = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex negH = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex expH = LaGenMatComplex::zeros(lattice_size, lattice_size);
    /* calculate H */
    H_generation(lattice_size, H);
    print_matrix(H, "H");
    /* calculate -H */
    matrix_negative(lattice_size, H, negH);
    print_matrix(negH, "-H");
    /* calculate exponentials */
    matrix_exponential_v(H, lattice_size, expH);
    print_matrix(expH, "e^(-H)");
}
void test_V(){
    /* initialise everything */
    int lattice_size = 5, time_size;
    LaGenMatComplex V = LaGenMatComplex::zeros(lattice_size, lattice_size);
    COMPLEX slice[lattice_size];
    double U = 1, beta = 10, lambda, delta_tau, mu = U / 2;

    /* calculate initial parameters */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size);
    cout << endl;

    /* generate the lattice */
    generate_slice(lattice_size, slice);
    print_array(slice, lattice_size, "slice");

    /* calculate V */
    V_calculation(slice, lattice_size, U, lambda, 1, delta_tau, mu, V);
    print_matrix(V, "V");
}
void test_B_calculation(){
    /* initialise everything */
    int lattice_size = 5, time_size;
    double U = 1, beta = 10, lambda, delta_tau, mu = U / 2;
    LaGenMatComplex B = LaGenMatComplex::zeros(lattice_size, lattice_size);
    COMPLEX slice[lattice_size];
    /* generate initial conditions */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size);
    /* generate time slice */
    generate_slice(lattice_size, slice);
    /* calculate B */
    cout << "sigma = 1" << endl;
    B_calculation_v(slice, lattice_size, U, lambda, 1, delta_tau, mu, B);
    B = 0;
    cout << "sigma = -1" << endl;
    B_calculation_v(slice, lattice_size, U, lambda, -1, delta_tau, mu, B);
    /* print result */
    print_matrix(B,"B = e^-H e^-V");
}
void test_B_calculation(const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* initialise everything */
    int lattice_size = 5, time_size;
    double U = 1, beta = 10, lambda, delta_tau, mu = U / 2;
    LaGenMatComplex B = LaGenMatComplex::zeros(lattice_size, lattice_size);
    COMPLEX slice[lattice_size];
    /* generate initial conditions */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, file);
    /* generate time slice */
    generate_slice(lattice_size, slice);
    /* calculate B */
    myfile << "sigma = 1" << endl;
    B_calculation(slice, lattice_size, U, lambda, 1, delta_tau, mu, B, file);
    B = 0;
    myfile << "sigma = -1" << endl;
    B_calculation(slice, lattice_size, U, lambda, -1, delta_tau, mu, B, file);
    /* print result */
    print_matrix(B,"B = e^-H e^-V", file);
    /* close the file */
    myfile.close();
}
void test_O(){
    /* initialise everything */
    int lattice_size = 5, time_size = 0;
    double U = 1, beta = 10, lambda, delta_tau, mu = U / 2;
    LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
    LaGenMatComplex O = LaGenMatComplex::zeros(lattice_size, lattice_size);
    /* generate initial conditions */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size);
    /* generate lattice */
    generate_lattice(lattice_size, time_size, lattice);
    print_matrix(lattice, "lattice");
    /* calculate O */
    O_calculation_v(lattice, lattice_size, time_size, U, lambda, 1, delta_tau, mu, O);
    print_matrix(O, "O");
}
void test_O(const string file){
    /* initialise everything */
    int lattice_size = 5, time_size = 0;
    double U = 1, beta = 10, lambda, delta_tau, mu = U / 2;
    LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
    LaGenMatComplex O = LaGenMatComplex::zeros(lattice_size, lattice_size);
    /* generate initial conditions */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, file);
    /* generate lattice */
    generate_lattice(lattice_size, time_size, lattice);
    print_matrix(lattice, "lattice", file);
    /* calculate O */
    O_calculation(lattice, lattice_size, time_size, U, lambda, 1, delta_tau, mu, O, file);
    print_matrix(O, "O", file);
}
void test_weight(){
    /* initialise stuff */
    int lattice_size = 5, time_size;
    double U = 1, beta = 10, lambda, delta_tau, mu = U / 2;
    COMPLEX weight;
    weight.r = 0;
    weight.i = 0;
    /* generate initial conditions */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size);
    /* generate lattice */
    LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
    generate_lattice(lattice_size, time_size, lattice);
    print_matrix(lattice, "lattice");
    /* calculate the weight */
    weight_calculation(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, weight);
    print_scalar(weight, "weight");
}
void test_weight(const string file){
    /* initialise stuff */
    int lattice_size = 5, time_size;
    double U = 1, beta = 10, lambda, delta_tau, mu = U / 2;
    COMPLEX weight;
    weight.r = 0;
    weight.i = 0;
    /* generate initial conditions */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, file);
    /* generate lattice */
    LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
    generate_lattice(lattice_size, time_size, lattice);
    print_matrix(lattice, "lattice", file);
    /* calculate the weight */
    weight_calculation_v(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, weight, file);
}
void test_imaginary_weight(const string file){
    /* initialise stuff */
    int lattice_size = 5, time_size;
    double U = 1, beta = 10, lambda, delta_tau, mu = U / 2;
    COMPLEX weight;
    clear_scalar(weight);

    /* generate initial conditions */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, file);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, "weight_i.txt");
    LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);

    for(int i = 0; i < 10; i++){
        /* generate lattice */
        generate_lattice(lattice_size, time_size, lattice);

        /* calculate total spin */
        calculate_total_spin(lattice, time_size, lattice_size, file);

        /* calculate weight */
        clear_scalar(weight);
        weight_calculation_v(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, weight, file);
        cout << weight << endl;
    }
}
void test_sweep(){
    /* initialise everything */
    int lattice_size = 5, time_size, iterations = 1000;// = 10000;
    double U = .1, beta = 1, lambda, delta_tau, mu = U / 2;
    double acceptance = 0, rejection = 0;
    /* generate initial conditions */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size);
    /* generate lattice */
    LaGenMatComplex lattice = LaGenMatComplex::zeros(time_size, lattice_size);
    // print_matrix(lattice, "intialised lattice");
    generate_lattice(lattice_size, time_size, lattice);
    print_matrix(lattice, "lattice");
    /* sweep the lattice */
    sweep_lattice_v(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, iterations, acceptance, rejection);
}
void test_increasing_U(){
    /* initialise everything */
    int lattice_size = 5, time_size = 0, iterations = 120;
    double U, beta = 5.0, lambda = 1.0, delta_tau, mu = U / 2;
    double acceptance = 0.0, rejection = 0.0;
    /* test U = 0 to 10 */
    for(int i = 1; i <= 10; i++){
        /* generate initial conditions */
        U = i;
        initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
        print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size);
        /* generate a lattice of spins */
        LaGenMatComplex lattice = LaGenMatComplex::zeros(time_size, lattice_size);
        generate_lattice(lattice_size, time_size, lattice);
        /* sweep the lattice */
        sweep_lattice_v(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, iterations, acceptance, rejection);
    }
}

/* ------ TO TEST ------ */
void sweep_lattice(LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, const double mu, const int iterations, double& acceptance, double& rejection, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);

    /* initialise everything */
    COMPLEX weightBefore;
    COMPLEX weightAfter;
    clear_scalar(weightBefore);
    clear_scalar(weightAfter);
    double probability = 0;
    string result;
    int count = 0;
    acceptance = 0;
    rejection = 0;

    /* output headings */
    myfile.width(11);
    myfile << "weight";
    myfile << " lattice" << endl;

    /* sweep through the lattice */
    for(int i = 0; i < iterations; i++){
        for(int t = 0; t < time_size; t++){
            for(int l = 0; l < lattice_size; l++){
                /* calculate the weight before the flip */
                weight_calculation(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, weightBefore);

                /* propose the flip */
                flip_spin(lattice, t, l, file);

                /* calculate the weight after the flip */
                weight_calculation(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, weightAfter);

                /* calculate the ratio of weights */
                probability = weightAfter.r / weightBefore.r;

                /* accept or reject the flip */
                double prob = random_double();
                if(abs(probability) >= 1){
                    result = "accepted";
                    acceptance++;
                }else{
                    if(probability > prob){
                        result = "accepted";
                        acceptance++;
                    }else{
                        flip_spin(lattice, t, l, file);
                        result = "rejected";
                        rejection++;
                    }
                }
                /* comments */
                    //for negative values, we do some integration
                    //P\to\tilde{P} = |P| and  F\to \tilde
                    //you have to multiply each quan you measure bu the sign
                count++;
                if(count%1000 == 0){
                    myfile << " (" << count <<") " << "[" << acceptance << "/" << rejection << "] " << result << " - probability: " << probability;
                    myfile.width(15);
                    myfile << " - weightBefore: " << weightBefore << ", weightAfter: " << weightAfter << endl;
                }
                // if(result == "accepted"){
                //     print_matrix(lattice);
                // }else{
                //     myfile << endl;
                // }
            }
            /* Comments */
                //when you take measurements, there is noise
                //we're doing marcov chain
                //the simplest quan we measure is double occupancy \bra n_up n_down \ket
        }
    }
    //results
        // with most parameters = 1, it stabilised at all -1 spins
    myfile << "["<< acceptance << "/" << rejection << "]" << endl;
    double acceptance_ratio = acceptance / (rejection + acceptance);
    myfile << "acceptance ratio = " << acceptance_ratio << endl;
    double percentage_acceptance = acceptance / rejection;
    myfile << "percentage acceptance = " << percentage_acceptance << endl << endl;
    /* close the file */
    myfile.close();
}
void test_sweep(const string file){
    /* initialise everything */
    int lattice_size = 5, time_size, iterations = 1000;// = 10000;
    double U = .1, beta = 1, lambda, delta_tau, mu;
    double acceptance = 0, rejection = 0;
    /* generate initial conditions */
    initial_parameter_calculation(U, beta, lambda, delta_tau, mu, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, file);
    /* generate lattice */
    LaGenMatComplex lattice = LaGenMatComplex::zeros(time_size, lattice_size);
    // print_matrix(lattice, "intialised lattice", file);
    generate_lattice(lattice_size, time_size, lattice);
    print_matrix(lattice, "lattice", file);
    /* sweep the lattice */
    sweep_lattice_v(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, iterations, acceptance, rejection);
}
void calculate_greens_function(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, const double mu, COMPLEX& weight, const string file){
    // calculates the single particle Green’s function

    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);

    /* initialise everything */
    LaGenMatComplex OUP = LaGenMatComplex::zeros(lattice_size,lattice_size);
    LaGenMatComplex ODN = LaGenMatComplex::zeros(lattice_size,lattice_size);
    LaGenMatComplex invUP = LaGenMatComplex::zeros(lattice_size,lattice_size);
    LaGenMatComplex invDN = LaGenMatComplex::zeros(lattice_size,lattice_size);
    COMPLEX detOUP;
    COMPLEX detODN;

    /* calculate O */
    myfile << "sigma = 1" << endl;
    O_calculation(lattice, lattice_size, time_size, U, lambda, 1, delta_tau, mu, OUP);
    print_matrix(OUP, "O UP", file);
    myfile << "sigma = -1" << endl;
    O_calculation(lattice, lattice_size, time_size, U, lambda, -1, delta_tau, mu, ODN, file);
    print_matrix(ODN, "O DN", file);

    /* calculate det(O) */
    matrix_determinant_e(lattice_size, OUP, detOUP);
    matrix_determinant_e(lattice_size, ODN, detODN);
    print_scalar(detOUP, "det(O UP)", file);
    print_scalar(detODN, "det(O DN)", file);
    myfile << endl;

    /* calculate O^-1 */
    matrix_inverse(OUP, lattice_size, invUP);
    matrix_inverse(ODN, lattice_size, invDN);
    print_matrix(invUP, "(O DN)^-1", file);
    print_matrix(invDN, "(O DN)^-1", file);

    /* close the file */
    myfile.close();
}
void update_greens_function(){
    //
}
void test_increasing_mu(const string file){
    /* initialise everything */
    int lattice_size = 5, time_size, iterations = 40;
    double U = 1, beta = 10, lambda, delta_tau, mu;
    /* generate initial conditions */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, file);
    /* generate lattice */
    LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
    generate_lattice(lattice_size, time_size, lattice);
    print_matrix(lattice, "lattice", file);
    /* plot mu */
    for(double i = 0; i < iterations; i++){
        mu = i * U / 8;
        /* generate initial conditions */
        print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, file);
        /* sweep across the lattice */
        /* calculate average spin */
        calculate_total_spin_f(lattice, time_size, lattice_size);
    }
}

/* ------ Main QMC Program ------ */
int main(){
    test_sweep("test.txt");
}
