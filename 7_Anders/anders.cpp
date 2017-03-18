#include <iostream> //cout
#include <fstream>
#include <string>
#include "anders.h"
#include <gmd.h> 	//LaGenMatDouble
#include <laslv.h>  //LUFactorizeIP, LaLUInverseIP, etc.
#include <blas3pp.h>
#include <random>   //random_device, mt19937
#include <cstdlib>	//rand, srand
#include <math.h>
#include <ctime>

using namespace std;

/* ------ NOTES ------ */
// be careful with row ordering
	// true = row ordered + DOES NOT link to the array it was made with
	// false = column ordered + LINKS to the array it was made with
// check that everything is reset to 0 or I where needed

/* ------ WORKING ------ */
/* -- Output -- */
void print_space(const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* log stuff */
    myfile << endl;
    /* close the file */
    myfile.close();
}
void print_text(const string text, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* log stuff */
    myfile << text << endl;
    /* close the file */
    myfile.close();
}
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
void print_scalar_f(const double scalar, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* print the scalar */
    myfile << scalar << endl;
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
void print_matrix_f(const LaGenMatComplex& matrix, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* print the matrix */
	myfile << matrix << endl;
    /* close the file */
    myfile.close();
}
void print_initial_parameters(const double U, const double beta, const double lambda, const double delta_tau, const double mu, const int time_size, const int lattice_size, const int iterations){
	cout << "no of lattice points = " << lattice_size << endl;
	cout << "no of time slices = " << time_size << endl;
	cout << "U = " << U << endl;
	cout << "beta = " << beta << endl;
	cout << "lambda = " << lambda << endl;
	cout << "delta tau = " << delta_tau << endl;
    cout << "mu = " << mu << endl;
    cout << "iterations = " << iterations << endl;
}
void print_initial_parameters(const double U, const double beta, const double lambda, double delta_tau, const double mu, const int time_size, const int lattice_size, const int iterations, const string file){
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
    myfile << "iterations = " << iterations << endl;
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
// - rearrange
void vec_to_diag(const LaVectorComplex& vector, const int array_size, LaGenMatComplex& diag){
    diag = 0;
    for(int i = 0; i < array_size; i++){
        diag(i, i) = vector(i);
    }
}
void clear_scalar(COMPLEX& scalar){
    scalar.r = 0;
    scalar.i = 0;
}
void isolate_row(const LaGenMatComplex& matrix, const int matrix_width, const int row, COMPLEX array[]){
    for(int i = 0; i < matrix_width; i++){
        array[i] = matrix(row, i);
    }
}
// - clear
void clear_array(COMPLEX array[], const int array_size){
    for(int i = 0; i < array_size; i++){
        array[i].r = 0;
        array[i].i = 0;
    }
}
// - change
void flip_scalar(COMPLEX& spin){
    spin.r = -spin.r;
    spin.i = -spin.i;
}
void flip_spin(LaGenMatComplex& lattice, const int t, const int l){
    lattice(t,l).r = -lattice(t,l).r;
    lattice(t,l).i = -lattice(t,l).i;
}

// Generation
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
string generate_file_name(const double U, const double beta, const int iterations, const string test){
    string UU =  "U" + to_string(U);
    string BB = "_B" + to_string(beta);
    string i  = "_i" + to_string(iterations);
    string t  = "_"  + test;
    return UU.substr(0,6)+ BB.substr(0,7) + i + t + ".txt";
}

// Calculation
// - generic
// -- scalar operations
double check_size(const double scalar){
    return floor(log10(scalar));
}
void scalar_division(const COMPLEX& A, const int B, COMPLEX& result){
    result.r = A.r / (double) B;
    result.i = A.i / (double) B;
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
    int iterations = 6;
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
// -- matrix operations
// --- simple arithmetic
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
// --- inverse
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
// --- exponentials
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
void diagonal_matrix_exponential(const LaGenMatComplex& matrix, const int matrix_size, LaGenMatComplex& result){
    result = 0;
    for(int i = 0; i < matrix_size; i++){
        scalar_exponential(matrix(i,i), result(i,i));
    }
}
// --- determinants
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
    LaGenMatComplex scaled_matrix = LaGenMatComplex::zeros(matrix_size, matrix_size);
    LaGenMatComplex cofactorMatrix;
    COMPLEX determinant;
    COMPLEX coefficient;
    cofactorMatrix = 0;
    /* test size of elements */
    double scale = check_size(matrix(0,0).r);
    /* scale matrix */
    matrix_multiple(matrix, matrix_size, 1 / scale, scaled_matrix);
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
            /* calculate the cofactor */
            generate_cofactor_matrix(matrix_size, matrix, i, cofactorMatrix);
            /* finish calculation */
            scalar_sum(determinant, scalar_multiple(coefficient, matrix_determinant(cofactor_size, cofactorMatrix)));
        }
        scalar_product(determinant, pow (10, scale * matrix_size));
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

// QMC
// - initial parameters
void initial_parameter_calculation(const double U, const double beta, double& lambda, double& delta_tau, int& time_size){
    // mu = U / 2;                             // by definition
    lambda = acoshf(exp(sqrt(0.125*U)/2));  // by definition
    time_size = ceil(beta / lambda);        // by definition
    delta_tau = beta / time_size;           // = sqrt(0.125 / U) by convension
}
// - H
void H_calculation(const int lattice_size, LaGenMatComplex& H){
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
// - V
void V_calculation(const COMPLEX slice[], const int lattice_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& V){
    /* reset V */
    V = 0;

    /* calculate V */
    V(0,0).r += (mu - U / 2);
    for(int i = 0; i < lattice_size; i++){
        V(i,i).r += lambda * sigma * slice[i].r / delta_tau;
    }
}

// - B
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
    H_calculation(lattice_size, H);
    V_calculation(slice, lattice_size, U, lambda, sigma, delta_tau, mu, V);

    /* calculate H + V */
    sum = H.copy();
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
    H_calculation(lattice_size, H);
    V_calculation(slice, lattice_size, U, lambda, sigma, delta_tau, mu, V);
    print_matrix(sum, "H");
    print_matrix(V, "V");

    /* calculate H + V */
    sum = H.copy();
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
void B_calculation_f(const COMPLEX slice[], const int lattice_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& B, const string file){
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
    H_calculation(lattice_size, H);
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
// - O
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
void O_calculation_f(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& O, const string file){
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
// - weight
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
    O_calculation(lattice, lattice_size, time_size, U, lambda,  1, delta_tau, mu, OUP);
    O_calculation(lattice, lattice_size, time_size, U, lambda, -1, delta_tau, mu, ODN);
    /* calculate det(O) */
    detOUP = matrix_determinant(lattice_size, OUP);
    detODN = matrix_determinant(lattice_size, ODN);
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
    detOUP = matrix_determinant(lattice_size, OUP);
    detODN = matrix_determinant(lattice_size, ODN);
    print_scalar(detOUP, "det(O UP)");
    print_scalar(detODN, "det(O DN)");
    /* calculate weight */
    weight = scalar_multiple(detOUP, detODN);
    print_scalar(weight, "weight");
}
void weight_calculation_f(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, const double mu, COMPLEX& weight, const string file){
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
    myfile << endl;
    /* close the file */
    myfile.close();
}
void weight_calculation_O(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, const double mu, LaGenMatComplex& OUP, COMPLEX& weight){
    /* initialise everything */
    OUP = LaGenMatComplex::zeros(lattice_size, lattice_size);
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
}
// - sweep
int judge_acceptance(const double probability){
    double ran = random_double();
    if(abs(probability) >= 1){
        return 1; // accept
    }else{
        if(probability > ran){
            return 1; // accept
        }else{
            return 0; // reject
        }
    }
}
void sweep_lattice(LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double beta, const double lambda, const double delta_tau, const double mu, const int iterations){
    /* initialise everything */
    COMPLEX weightBefore, weightAfter;
    double probability = 0;
    string result;
    int count = 0, acceptance = 0, rejection = 0;
    int total_count = lattice_size * time_size * iterations;
    string rf = generate_file_name(U, beta, iterations, "results");
    string wf = generate_file_name(U, beta, iterations, "weights");
    string sf = generate_file_name(U, beta, iterations, "spins");
    string df = generate_file_name(U, beta, iterations, "occupancy");
    string nf = generate_file_name(U, beta, iterations, "n");
    LaGenMatComplex O = LaGenMatComplex::zeros(lattice_size, lattice_size);

    /* output initial conditions */
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, iterations, rf);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, iterations, df);
    print_matrix(lattice, "lattice", rf);

    /* sweep through the lattice */
    for(int i = 0; i < iterations; i++){
        for(int t = 0; t < time_size; t++){
            for(int l = 0; l < lattice_size; l++){
                clear_scalar(weightBefore);
                clear_scalar(weightAfter);
                count++;

                /* calculate the weight before the flip */
                if(count == 1){
                    weight_calculation_O(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, O, weightBefore);
                }else{
                    weightBefore = weightAfter;
                }

                /* propose the flip */
                flip_spin(lattice, t, l);

                /* calculate the weight after the flip */
                weight_calculation(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, weightAfter);

                /* calculate the ratio of weights */
                probability = weightAfter.r / weightBefore.r;

                /* accept or reject the flip */
                if(judge_acceptance(probability) == 1){
                    acceptance++;
                }else{
                    flip_spin(lattice, t, l);
                    weightAfter = weightBefore;
                    rejection++;
                }

                /* output results */
                if(count % (total_count / 300) == 0){
                    measure_result(count, acceptance, rejection, result, probability, rf);
                    measure_weight(count, probability, weightBefore, weightAfter, wf);
                    measure_spin(lattice, time_size, lattice_size, sf);
                    measure_double_occcupancy_ii(2, O, lattice_size, df);
                    measure_n(O, lattice_size, nf);
                }
            }
        }
    }
    print_space(rf);
    measure_acceptance(acceptance, rejection, total_count, rf);
}
// - measurements
// -- execution time
void measure_execution_time(const int iterations, const int start_s, const int stop_s, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* initialise everything */
    double execution_time = (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
    /* log stuff */
    myfile << "iterations = " << iterations << " - ";
    myfile << "execution time = " << execution_time << endl;
    cout << "iterations = " << iterations << " - ";
    cout << "execution time = " << execution_time << endl;
    /* close the file */
    myfile.close();
}
// -- acceptance
void measure_result(const int count, const int acceptance, const int rejection, const string result, const double probability, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* log results */
    myfile << " (" << count <<") ";
    myfile << "[" << acceptance << "/" << rejection << "] " << result;
    myfile << " - probability: " << probability;
    myfile << endl;
    /* close the file */
    myfile.close();
}
void measure_acceptance(const int acceptance, const int rejection, const int total_count, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* calculate stuff */
    double acceptance_ratio = (double) acceptance / (double) rejection;
    double percentage_acceptance = (double) acceptance / (double) total_count;
    /* log stuff */
    myfile << "["<< acceptance << "/" << rejection << "] - ";
    myfile << "acceptance ratio = " << acceptance_ratio << " - ";
    myfile << "percentage acceptance = " << percentage_acceptance << endl << endl;
    /* close the file */
    myfile.close();
}
// -- spin
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
    print_scalar(average_spin, "", "av_spin.txt");
    // myfile << average_spin << endl << endl;
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
void measure_spin(const LaGenMatComplex& lattice, const int time_size, const int lattice_size, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* initialise everything */
    double lattice_volume = (double) time_size * (double) lattice_size;
    /* calculate total spin */
    double total_spin = 0;
    for(int t = 0; t < time_size; t++){
        for(int l = 0; l < lattice_size; l++){
            total_spin += lattice(t,l).r;
        }
    }
    double average_spin = total_spin / lattice_volume;
    // cout << average_spin << endl;
    /* log stuff */
    myfile << average_spin << endl;
    /* close the file */
    myfile.close();
}
// -- weight
void measure_weight(const int count, const double probability, const COMPLEX weightBefore, const COMPLEX weightAfter, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* log results */
    myfile << " (" << count << ") - probability = " << probability << " - ";
    myfile << "weightBefore = " << weightBefore << " - ";
    myfile << "weightAfter = "  << weightAfter  << endl;
    /* close the file */
    myfile.close();
}
void measure_av_weight(){
    //!!!
}
// -- occupancy
void measure_double_occcupancy(const LaGenMatComplex& O, const int lattice_size, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* initialise stuff */
    LaGenMatComplex double_occcupancy = LaGenMatComplex::zeros(lattice_size, lattice_size);
    /* process stuff */
    matrix_inverse(O, lattice_size, double_occcupancy);
    /* log stuff */
    print_matrix_f(double_occcupancy, file);
    /* close the file */
    myfile.close();
}
double double_occupancy_ii(const int i, const LaGenMatComplex& O, const int lattice_size){
    /* initialise stuff */
    LaGenMatComplex double_occcupancy = LaGenMatComplex::zeros(lattice_size, lattice_size);
    /* calculate double occcupancy */
    matrix_inverse(O, lattice_size, double_occcupancy);
    return double_occcupancy(i,i).r;
}
void measure_double_occcupancy_ii(const int i, const LaGenMatComplex& O, const int lattice_size, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* initialise stuff */
    LaGenMatComplex double_occcupancy = LaGenMatComplex::zeros(lattice_size, lattice_size);
    /* process stuff */
    matrix_inverse(O, lattice_size, double_occcupancy);
    /* log stuff */
    // cout << double_occcupancy(i,j) << " - " << double_occcupancy(i,j).r << endl;
    print_scalar_f(double_occcupancy(i,i).r, file);
    /* close the file */
    myfile.close();
}
double n(const LaGenMatComplex& O, const int lattice_size){
    /* initialise stuff */
    double n = 0;
    LaGenMatComplex double_occcupancy = LaGenMatComplex::zeros(lattice_size, lattice_size);
    /* calculate double occcupancy */
    matrix_inverse(O, lattice_size, double_occcupancy);
    for(int i = 0; i < lattice_size; i++){
        n += double_occcupancy(i,i).r;
    }
    return n;
}
void measure_n(const LaGenMatComplex& O, const int lattice_size, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* initialise stuff */
    LaGenMatComplex double_occcupancy = LaGenMatComplex::zeros(lattice_size, lattice_size);
    double n = 0;
    /* process stuff */
    matrix_inverse(O, lattice_size, double_occcupancy);
    /* log stuff */
    for(int i = 0; i < lattice_size; i++){
        n += double_occcupancy(i,i).r;
    }
    print_scalar_f(n, file);
    /* close the file */
    myfile.close();
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
void test_scalar_exponential(){
    COMPLEX scalar, result;
    scalar.r = random_int(9);
    scalar.i = random_int(9);
    scalar_exponential(scalar, result);
    cout << "e^" << scalar << " = " << result << endl;
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
    LaGenMatComplex matrix = LaGenMatComplex::rand(4,4,0,4);
    LaGenMatComplex scaled_matrix = LaGenMatComplex::zeros(matrix_size, matrix_size);
    COMPLEX result;
    clear_scalar(result);
    /* calculate determinant */
    print_matrix(matrix, "initial matrix");
    print_scalar(matrix_determinant(4, matrix), "matrix determinant");
    /* scale the matrix */
    matrix_multiple(matrix, matrix_size, pow(10.0,scale), scaled_matrix);
    print_matrix(scaled_matrix, "scaled matrix");
    print_scalar(matrix_determinant(4, scaled_matrix), "my determinant");
}
// - qmc
void test_H(){
    /* start timing */
    int start_s = clock();

    /* initialise everything */
    string file = "test_H";
    int lattice_size = 5;
    LaGenMatComplex H;
    LaVectorComplex eigenvalues = LaVectorComplex(lattice_size);
    LaGenMatComplex eigenvectors = LaGenMatComplex::zeros(lattice_size, lattice_size);

    /* generate matrix */
    H_calculation(lattice_size, H);
    print_matrix(H);

    /* stop timing */
    int stop_s = clock();
    measure_execution_time(1, start_s, stop_s, file);

    /* check with eigenstuff */
    LaEigSolve(H, eigenvalues, eigenvectors);
    print_vector(eigenvalues, "eigenvalues");
    // eigenvalues are 2 cos(n pi / q), where q = the matrix size
}
void test_V(){
    /* start timing */
    int start_s = clock();

    /* initialise everything */
    string file = "test_V";
    int lattice_size = 5, time_size;
    LaGenMatComplex V = LaGenMatComplex::zeros(lattice_size, lattice_size);
    COMPLEX slice[lattice_size];
    double U = 1, beta = 10, lambda, delta_tau, mu = U;

    /* calculate initial parameters */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, 1);
    cout << endl;

    /* generate the lattice */
    generate_slice(lattice_size, slice);
    print_array(slice, lattice_size, "slice");

    /* calculate V */
    V_calculation(slice, lattice_size, U, lambda, 1, delta_tau, mu, V);
    print_matrix(V, "V");

    /* stop timing */
    int stop_s = clock();
    measure_execution_time(1, start_s, stop_s, file);
}
void test_B(){
    /* open the file */
    string file = "test_B.txt";
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* initialise everything */
    int lattice_size = 5, time_size;
    double U = 1, beta = 10, lambda, delta_tau, mu = U / 2;
    LaGenMatComplex B = LaGenMatComplex::zeros(lattice_size, lattice_size);
    COMPLEX slice[lattice_size];
    /* generate initial conditions */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, 1, file);
    /* generate time slice */
    generate_slice(lattice_size, slice);
    /* calculate B */
    myfile << "sigma = 1" << endl;
    B_calculation_f(slice, lattice_size, U, lambda, 1, delta_tau, mu, B, file);
    B = 0;
    myfile << "sigma = -1" << endl;
    B_calculation_f(slice, lattice_size, U, lambda, -1, delta_tau, mu, B, file);
    /* print result */
    print_matrix(B,"B = e^-H e^-V", file);
    /* close the file */
    myfile.close();
}
void test_O(){
    /* initialise everything */
    string file = "test_O.txt";
    int lattice_size = 5, time_size = 0;
    double U = 1, beta = 10, lambda, delta_tau, mu = U / 2;
    LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
    LaGenMatComplex O = LaGenMatComplex::zeros(lattice_size, lattice_size);
    /* generate initial conditions */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, 1, file);
    /* generate lattice */
    generate_lattice(lattice_size, time_size, lattice);
    print_matrix(lattice, "lattice", file);
    /* calculate O */
    O_calculation_f(lattice, lattice_size, time_size, U, lambda, 1, delta_tau, mu, O, file);
    print_matrix(O, "O", file);
}
void test_weight(){
    /* initialise stuff */
    string file = "test_weight.txt";
    int lattice_size = 5, time_size;
    double U = 1, beta = 10, lambda, delta_tau, mu = U / 2;
    COMPLEX weight;
    weight.r = 0;
    weight.i = 0;
    /* generate initial conditions */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, 1, file);
    /* generate lattice */
    LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
    generate_lattice(lattice_size, time_size, lattice);
    print_matrix(lattice, "lattice", file);
    /* calculate the weight */
    weight_calculation_v(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, weight, file);
}
void test_sweep(){
    /* initialise everything */
    string file = "test_sweep.txt";
    int lattice_size = 5, time_size, iterations = 1000;// = 10000;
    double U = .1, beta = 1, lambda, delta_tau, mu = U / 2;
    double acceptance = 0, rejection = 0;
    /* generate initial conditions */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, iterations);
    /* generate lattice */
    LaGenMatComplex lattice = LaGenMatComplex::zeros(time_size, lattice_size);
    // print_matrix(lattice, "intialised lattice");
    generate_lattice(lattice_size, time_size, lattice);
    print_matrix(lattice, "lattice");
    /* sweep the lattice */
    sweep_lattice(lattice, lattice_size, time_size, U, beta, lambda, delta_tau, mu, iterations);
}
void test_sweep_d(){
    /* initialise everything */
    int start_s = clock();
    int lattice_size = 5, time_size, iterations = 10000;
    double U = .1, beta = 1, lambda, delta_tau, mu = U/2;
    /* generate initial conditions */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    /* generate lattice */
    LaGenMatComplex lattice = LaGenMatComplex::zeros(time_size, lattice_size);
    generate_lattice(lattice_size, time_size, lattice);
    /* sweep the lattice */
    sweep_lattice(lattice, lattice_size, time_size, U, beta, lambda, delta_tau, mu, iterations);
    /* check time */
    int stop_s = clock();
    cout << "time: " << (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000 << endl;
}
void test_sweep_f_t(){
    /* initialise everything */
    int lattice_size = 5, time_size, iterations = 0;
    double U = .1, beta = 1, lambda, delta_tau, mu = U / 2;
    string file = generate_file_name(U, beta, iterations, "weight-time");
    /* do stuff */
    for(int i = 1; i <= 4; i++){
        /* generate initial conditions */
        iterations = pow(10.0, i);
        initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
        int start_s = clock();
        /* generate lattice */
        LaGenMatComplex lattice = LaGenMatComplex::zeros(time_size, lattice_size);
        generate_lattice(lattice_size, time_size, lattice);/* sweep the lattice */
        sweep_lattice(lattice, lattice_size, time_size, U, beta, lambda, delta_tau, mu, iterations);
        /* output results */
        int stop_s = clock();
        measure_execution_time(iterations, start_s, stop_s, file);
    }
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
        print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, iterations);
        /* generate a lattice of spins */
        LaGenMatComplex lattice = LaGenMatComplex::zeros(time_size, lattice_size);
        generate_lattice(lattice_size, time_size, lattice);
        /* sweep the lattice */
        sweep_lattice(lattice, lattice_size, time_size, U, beta, lambda, delta_tau, mu, iterations);
    }
}

/* ------ TO TEST ------ */
void measure_stuff(const int stuff, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* initialise stuff */
    /* log stuff */
    myfile << stuff << endl;
    /* close the file */
    myfile.close();
}
void test_increasing_mu(const string file){
    /* initialise everything */
    int lattice_size = 5, time_size, iterations = 40;
    double U = 1, beta = 10, lambda, delta_tau, mu;
    /* generate initial conditions */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, iterations, file);
    /* generate lattice */
    LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
    generate_lattice(lattice_size, time_size, lattice);
    print_matrix(lattice, "lattice", file);
    /* plot mu */
    for(double i = 0; i < iterations; i++){
        mu = i * U / 8;
        /* generate initial conditions */
        print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, iterations, file);
        /* sweep across the lattice */
        /* calculate average spin */
        calculate_total_spin_f(lattice, time_size, lattice_size);
    }
}

// Calculation
// - qmc

/* ------ TO IMPLEMENT ------ */
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
    O_calculation_f(lattice, lattice_size, time_size, U, lambda, -1, delta_tau, mu, ODN, file);
    print_matrix(ODN, "O DN", file);

    /* calculate det(O) */
    detOUP = matrix_determinant(lattice_size, OUP);
    detODN = matrix_determinant(lattice_size, ODN);
    print_scalar(detOUP, "det(O UP)", file);
    print_scalar(, "det(O DN)", file);
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
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, iterations, file);
    /* generate lattice */
    LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
    generate_lattice(lattice_size, time_size, lattice);
    print_matrix(lattice, "lattice", file);
    /* plot mu */
    for(double i = 0; i < iterations; i++){
        mu = i * U / 8;
        /* generate initial conditions */
        print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, iterations, file);
        /* sweep across the lattice */
        /* calculate average spin */
        calculate_total_spin_f(lattice, time_size, lattice_size);
    }
}

/* ------ Main QMC Program ------ */
int main(){
    test_matrix_determinant();
}
