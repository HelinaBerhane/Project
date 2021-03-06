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
// - to terminal
void print_scalar(const COMPLEX scalar){
    cout << scalar << endl;
}
void print_scalar(const COMPLEX scalar, const string name){
    cout << name << ": " << scalar << endl;
}
void print_scalar(const double scalar, const string name){
    cout << name << ": " << scalar << endl;
}
void print_array(const COMPLEX array[], int array_size, const string name){
	cout << name << ": ";
    for(int i = 0; i < array_size; i++){
        cout << array[i] << " ";
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
void print_scalar_f(const double scalar, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* print the scalar */
    myfile << scalar << endl;
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
// - to file
void print_space(const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* log stuff */
    myfile << " ";
    /* close the file */
    myfile.close();
}
void print_delimater(const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* log stuff */
    myfile << endl;
    /* close the file */
    myfile.close();
}
void print_line(const string file){
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
void print_double(const double number, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* log stuff */
    myfile << number << endl;
    /* close the file */
    myfile.close();
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
string generate_file_name(const double U, const double beta, const double mu, const int iterations, const string test){
    string UU =  "U" + to_string(U);
    string BB = "_B" + to_string(beta);
    string m  = "_m" + to_string(mu);
    string i  = "_i" + to_string(iterations);
    string t  = "_"  + test;
    return UU.substr(0,6) + BB.substr(0,7) + m.substr(0,7) + i + t + ".txt";
}

// Calculation
// - generic
// -- scalar operations
void scalar_division(const COMPLEX& A, const int B, COMPLEX& result){
    result.r = A.r / (double) B;
    result.i = A.i / (double) B;
}
void scalar_division(const COMPLEX& A, const COMPLEX& B, COMPLEX& result){
    la::complex<double> laA = la::complex<double>(A); //convert to la::complex<double>
    la::complex<double> laB = la::complex<double>(B);
    la::complex<double> laResult = la::complex<double>(result);
    laResult = laA / laB;
    result = laResult.toCOMPLEX();
}
COMPLEX scalar_division(const COMPLEX& A, const COMPLEX& B){
    COMPLEX result;
    //convert to la::complex<double>
    la::complex<double> laA = la::complex<double>(A);
    la::complex<double> laB = la::complex<double>(B);
    la::complex<double> laResult = la::complex<double>(result);
    laResult = laA / laB;
    result = laResult.toCOMPLEX();
    return result;
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
COMPLEX scalar_addition(COMPLEX& A, COMPLEX& B){
    COMPLEX result;
    result.r = A.r + B.r;
    result.i = A.i + B.i;
    return result;
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
void matrix_product(LaGenMatComplex& matrix, const double number, const int matrix_size){
    for(int i = 0; i < matrix_size; i++){
        for(int j = 0; j < matrix_size; j++){
            matrix(i, j).r *= number;
            matrix(i, j).i *= number;
        }
    }
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
// --- determinants !!!
void triangle_matrix_v(const LaGenMatComplex& matrix, const int matrix_size, LaGenMatComplex& triangle){
    print_matrix(matrix, "initial matrix");

    /* initialise everything */
    triangle = matrix.copy();

    // for every row but the first
    for(int i = 1; i < matrix_size; i++){

        cout << "i = " << i << endl << endl;

        for(int row = i; row < matrix_size; row++){

            if(triangle(i-1, i-1).r != 0 || triangle(i-1, i-1).i != 0){

                COMPLEX multiple = scalar_division(triangle(row, i-1), triangle(i-1, i-1));
                cout << "multiple = " << triangle(row, i-1) << " / " << triangle(i-1, i-1) << " = " << multiple << endl;

                cout << "initial row = ";
                for(int column = 0; column < matrix_size; column++){
                    cout << triangle(row, column);
                }
                cout << endl << endl;

                cout << "subtraction = " << endl;
                for(int column = 0; column < matrix_size; column++){
                    cout << multiple << " * " << triangle(i-1, column) << endl;
                }
                cout << endl;

                for(int column = 0; column < matrix_size; column++){
                    COMPLEX subtraction = scalar_multiple(triangle(i-1, column), multiple);
                    // cout << subtraction;
                    triangle(row, column).r -= subtraction.r;
                    triangle(row, column).i -= subtraction.i;
                }

                cout << "new row = ";
                for(int column = 0; column < matrix_size; column++){
                    cout << triangle(row, column);
                }
                cout << endl << endl;

            }else{
                cout << "n / a" << endl;
            }
        }
    }

    print_matrix(triangle, "triangle matrix");
}
void triangle_matrix(const LaGenMatComplex& matrix, const int matrix_size, LaGenMatComplex& triangle){
    /* initialise everything */
    triangle = matrix.copy();

    // for every row but the first
    for(int i = 1; i < matrix_size; i++){
        for(int row = i; row < matrix_size; row++){
            if(triangle(i-1, i-1).r != 0 || triangle(i-1, i-1).i != 0){
                COMPLEX multiple = scalar_division(triangle(row, i-1), triangle(i-1, i-1));
                for(int column = 0; column < matrix_size; column++){
                    COMPLEX subtraction = scalar_multiple(triangle(i-1, column), multiple);
                    triangle(row, column).r -= subtraction.r;
                    triangle(row, column).i -= subtraction.i;
                }
            }
        }
    }
}
COMPLEX matrix_determinant_v(const LaGenMatComplex& matrix, const int matrix_size){
    /* initialise everything */
    LaGenMatComplex triangle = LaGenMatComplex::zeros(matrix_size, matrix_size);
    COMPLEX product;
    product.r = 1;
    product.i = 0;
    /* calculate triangle matrix */
    triangle_matrix(matrix, matrix_size, triangle);
    print_matrix(triangle, "triangle matrix");
    /* calculate determinant */
    for(int i = 0; i < matrix_size; i++){
        scalar_product(product, triangle(i, i));
    }
    print_scalar(product, "determinant");
    return product;
}
COMPLEX matrix_determinant(const LaGenMatComplex& matrix, const int matrix_size){
    /* initialise everything */
    LaGenMatComplex triangle = LaGenMatComplex::zeros(matrix_size, matrix_size);
    COMPLEX product;
    product.r = 1;
    product.i = 0;
    /* calculate triangle matrix */
    triangle_matrix(matrix, matrix_size, triangle);
    /* calculate determinant */
    for(int i = 0; i < matrix_size; i++){
        scalar_product(product, triangle(i, i));
    }
    return product;
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
    for(int i = 0; i < lattice_size; i++){
        V(i,i).r += lambda * sigma * slice[i].r / delta_tau + (mu - U / 2);
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
    detOUP = matrix_determinant(OUP, lattice_size);
    detODN = matrix_determinant(ODN, lattice_size);
    /* calculate weight */
    weight = scalar_multiple(detOUP, detODN);
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
    detOUP = matrix_determinant(OUP, lattice_size);
    detODN = matrix_determinant(ODN, lattice_size);
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
void weight_calculation_O(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, const double mu, LaGenMatComplex& OUP, LaGenMatComplex& ODN, COMPLEX& weight){
    /* initialise everything */
    OUP = LaGenMatComplex::zeros(lattice_size, lattice_size);
    ODN = LaGenMatComplex::zeros(lattice_size,lattice_size);
    COMPLEX detOUP;
    COMPLEX detODN;
    clear_scalar(weight);
    clear_scalar(detOUP);
    clear_scalar(detODN);
    /* calculate O */
    O_calculation(lattice, lattice_size, time_size, U, lambda,  1, delta_tau, mu, OUP);
    O_calculation(lattice, lattice_size, time_size, U, lambda, -1, delta_tau, mu, ODN);
    /* calculate det(O) */
    detOUP = matrix_determinant(OUP, lattice_size);
    detODN = matrix_determinant(ODN, lattice_size);
    /* calculate weight */
    weight = scalar_multiple(detOUP, detODN);
}
// - sweep
int judge_acceptance(const double probability){
    double ran = random_double();
    if(probability >= 1){
        return 1; // accept
    }else{
        if(probability < 0){
            cout << "negative sign" << endl;
            return 0; // reject
        }else{
            if(probability > ran){
                return 1; // accept
            }else{
                return 0; // reject
            }
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
    string rf = generate_file_name(U, beta, mu, iterations, "results");
    string wf = generate_file_name(U, beta, mu, iterations, "weights");
    string sf = generate_file_name(U, beta, mu, iterations, "spins");
    string asf = generate_file_name(U, beta, mu, iterations, "average_spin");
    string df = generate_file_name(U, beta, mu, iterations, "occupancy");
    string nf = generate_file_name(U, beta, mu, iterations, "n");
    string cf = generate_file_name(U, beta, mu, iterations, "charge_density");
    string file = generate_file_name(U, beta, mu, iterations, "measurements");
    LaGenMatComplex OUP = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex ODN = LaGenMatComplex::zeros(lattice_size, lattice_size);
    /* output initial conditions */
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, iterations, rf);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, iterations, df);
    print_text("U", file);
    print_double(U, file);
    print_text("beta", file);
    print_double(beta, file);
    print_text("mu", file);
    print_double(mu, file);
    print_text("iterations", file);
    print_double(iterations, file);

    print_matrix(lattice, "lattice", rf);

    /* sweep through the lattice */
    double av_spin = 0, spin_sign = 0, av_charge = 0;
    double tot = 0;
    for(int i = 0; i < iterations; i++){
        for(int t = 0; t < time_size; t++){
            for(int l = 0; l < lattice_size; l++){
                count++;
                clear_scalar(weightBefore);
                clear_scalar(weightAfter);

                /* calculate the weight before the flip */
                // if(count == 1){
                    // cout << "recalculating" << count << endl;
                weight_calculation_O(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, OUP, ODN, weightBefore);
                // }else{
                //     weightBefore.r = weightAfter.r;
                //     weightBefore.i = weightAfter.i;
                //     clear_scalar(weightAfter);
                // }

                /* propose the flip */
                flip_spin(lattice, t, l);

                /* calculate the weight after the flip */
                weight_calculation(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, weightAfter);


                // if(count % (total_count / 500) == 0){
                //     measure_weight(count, probability, weightBefore, weightAfter, wf);
                // }

                /* calculate the ratio of weights */
                probability  = weightAfter.r / weightBefore.r;
                double sgn   = probability / abs(probability);
                probability *= sgn;
                /* accept or reject the flip */
                if(judge_acceptance(probability) == 1){
                    acceptance++;
                    result = "accepted";
                }else{
                    flip_spin(lattice, t, l);
                    // weightAfter.r = weightBefore.r;
                    // weightAfter.i = weightBefore.i;
                    // clear_scalar(weightBefore);
                    rejection++;
                }

                /* output results */
                if(count % (total_count / 500) == 0){
                    cout << ".";
                    measure_result(count, acceptance, rejection, result, probability, rf);
                    print_double(average_spin(lattice, time_size, lattice_size), sf);
                    measure_double_occcupancy_ii(2, OUP, lattice_size, df, sgn);
                    measure_n(OUP, lattice_size, nf, sgn);
                    if(count > (total_count / 50)){
                        av_spin += average_spin(lattice, time_size, lattice_size);
                        tot++;
                        if(av_spin < 0){
                            spin_sign++;
                        }
                        av_charge += charge_density(OUP, ODN, lattice_size, cf, sgn);
                    }
                }
            }
        }
    }
    av_spin /= tot;
    double av_spin_sign = spin_sign / tot;
    av_charge /= tot;
    cout << endl;
    cout << "average spin = " << av_spin << endl;
    print_text("U", asf);
    print_text("acceptance", file);
    print_text("average spin", file);
    print_text("average spin", asf);
    print_text("average spin sign", file);
    print_text("average spin sign", asf);
    print_text("average charge density", file);
    print_double(U, asf);
    print_double(av_spin, file);
    print_double(av_spin, asf);
    print_double(av_spin_sign, file);
    print_double(av_spin_sign, asf);
    print_double(av_charge, file);
    print_line(rf);
    measure_acceptance(acceptance, rejection, total_count, rf);
    measure_acceptance(acceptance, rejection, total_count, file);
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
    myfile << " (" << count <<") " << result << " - ";
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
double average_spin(const LaGenMatComplex& lattice, const int time_size, const int lattice_size){
    /* calculate total spin */
    double total_spin = 0;
    double lattice_volume = (double) time_size * (double) lattice_size;
    for(int t = 0; t < time_size; t++){
        for(int l = 0; l < lattice_size; l++){
            total_spin += lattice(t,l).r;
        }
    }
    return total_spin / lattice_volume;
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
void measure_double_occcupancy(const LaGenMatComplex& O, const int lattice_size, const string file, const double sign){
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
double double_occupancy_ii(const int i, const LaGenMatComplex& O, const int lattice_size, const double sign){
    /* initialise stuff */
    LaGenMatComplex double_occcupancy = LaGenMatComplex::zeros(lattice_size, lattice_size);
    /* calculate double occcupancy */
    matrix_inverse(O, lattice_size, double_occcupancy);
    return double_occcupancy(i,i).r * sign;
}
void measure_double_occcupancy_ii(const int i, const LaGenMatComplex& O, const int lattice_size, const string file, const double sign){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* initialise stuff */
    LaGenMatComplex double_occcupancy = LaGenMatComplex::zeros(lattice_size, lattice_size);
    /* process stuff */
    matrix_inverse(O, lattice_size, double_occcupancy);
    /* log stuff */
    // cout << double_occcupancy(i,j) << " - " << double_occcupancy(i,j).r << endl;
    print_scalar_f(double_occcupancy(i,i).r * sign, file);
    /* close the file */
    myfile.close();
}
double n(const LaGenMatComplex& O, const int lattice_size, const double sign){
    /* initialise stuff */
    double n = 0;
    LaGenMatComplex double_occcupancy = LaGenMatComplex::zeros(lattice_size, lattice_size);
    /* calculate double occcupancy */
    matrix_inverse(O, lattice_size, double_occcupancy);
    for(int i = 0; i < lattice_size; i++){
        n += double_occcupancy(i,i).r * sign;
    }
    return n;
}
double charge_density(const LaGenMatComplex& OUP, const LaGenMatComplex& ODN, const int lattice_size, const string file, const double sgn){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* calculate charge density */
    double density = n(OUP, lattice_size, sgn) + n(ODN, lattice_size, sgn);
    print_double(density, file);
    /* close the file */
    myfile.close();
    return density;
}
void measure_n(const LaGenMatComplex& O, const int lattice_size, const string file, const double sign){
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
        n += double_occcupancy(i,i).r * sign;
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
void test_triangle_matrix(){
    /* initialise everything */
    double matrix_size = 4, scale = 8;
    LaGenMatComplex matrix = LaGenMatComplex::zeros(matrix_size, matrix_size);
    LaGenMatComplex triangle = LaGenMatComplex::zeros(matrix_size, matrix_size);
    /* generate matrix */
    for(int i = 0; i < matrix_size; i++){
        for(int j = 0; j < matrix_size; j++){
            matrix(i,j).r = random_int(5) * pow(10, scale);
            matrix(i,j).i = 0;
        }
    }
    print_matrix(matrix, "initial matrix");
    /* calculate triangle matrix */

    triangle_matrix(matrix, matrix_size, triangle);
    print_matrix(triangle, "triangle matrix");
}
void test_matrix_determinant(){
    /* initialise everything */
    double matrix_size = 4, scale = 8;
    LaGenMatComplex matrix = LaGenMatComplex::zeros(matrix_size, matrix_size);
    /* generate matrix */
    for(int i = 0; i < matrix_size; i++){
        for(int j = 0; j < matrix_size; j++){
            matrix(i,j).r = random_int(5);
            matrix(i,j).i = 0;
        }
    }
    print_matrix(matrix, "initial matrix");
    /* calculate determinant */
    print_scalar(matrix_determinant(matrix, matrix_size), "determinant");
    /* generate matrix */
    for(int i = 0; i < matrix_size; i++){
        for(int j = 0; j < matrix_size; j++){
            matrix(i,j).r = random_int(5) * pow(10, scale);
            matrix(i,j).i = 0;
        }
    }
    print_matrix(matrix, "initial matrix");
    /* calculate determinant */
    print_scalar(matrix_determinant(matrix, matrix_size), "large determinant");
    /* generate matrix */
    for(int i = 0; i < matrix_size; i++){
        for(int j = 0; j < matrix_size; j++){
            matrix(i,j).r = random_double() * pow(10, scale);
            matrix(i,j).i = 0;
        }
    }
    print_matrix(matrix, "random matrix");
    /* calculate determinant */
    print_scalar(matrix_determinant(matrix, matrix_size), "random determinant");
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
    weight_calculation_f(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, weight, file);
}
void test_sweep(){
    /* initialise everything */
    int start_s = clock();
    int lattice_size = 2, time_size, iterations = 500;
    double U = 0.5, beta = 1, lambda, delta_tau, mu = U;
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
void test_sweep(const int lattice_size, const int iterations, const double U, const double beta, const double mu){
    /* initialise everything */
    int start_s = clock();
    int time_size;
    double lambda, delta_tau;
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
void test_increasing_U(){
    /* initialise everything */
    int lattice_size = 5, time_size, iterations = 10000;
    double U, beta = 10.0, lambda, delta_tau, mu;;
    LaGenMatComplex lattice;
    /* test U = 0 to 10 */
    for(int i = 1; i <= 5; i++){
        /* generate initial conditions */
        U = 2*i, mu = U / 2;
        int start_s = clock();
        initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
        print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, iterations);
        /* generate a lattice of spins */
        lattice = LaGenMatComplex::zeros(time_size, lattice_size);
        generate_lattice(lattice_size, time_size, lattice);
        /* sweep the lattice */
        sweep_lattice(lattice, lattice_size, time_size, U, beta, lambda, delta_tau, mu, iterations);
        /* check time */
        int stop_s = clock();
        cout << "time: " << (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000 << endl;
    }
}
void test_increasing_mu(){
    /* initialise everything */
    int lattice_size = 5, time_size, iterations = 5000;
    double U, beta, lambda, delta_tau, mu;
    /* plot mu */
    /* generate initial conditions */
    U = 8, beta = 1;
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, iterations);
    for(double i = 0; i < 32; i++){
        LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
        generate_lattice(lattice_size, time_size, lattice);
        mu = i * U / 8;
        cout << "mu = " << mu << endl;
        /* sweep across the lattice */
        sweep_lattice(lattice, lattice_size, time_size, U, beta, lambda, delta_tau, mu, iterations);
    }
    U = 4, beta = 5.76;
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, iterations);
    for(double i = 0; i < 32; i++){
        LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
        generate_lattice(lattice_size, time_size, lattice);
        mu = i * U / 8;
        cout << "mu = " << mu << endl;
        /* sweep across the lattice */
        sweep_lattice(lattice, lattice_size, time_size, U, beta, lambda, delta_tau, mu, iterations);
    }
    U = 4, beta = 2;
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, iterations);
    for(double i = 0; i < 16; i++){
        LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
        generate_lattice(lattice_size, time_size, lattice);
        mu = i * U / 8;
        cout << "mu = " << mu << endl;
        /* sweep across the lattice */
        sweep_lattice(lattice, lattice_size, time_size, U, beta, lambda, delta_tau, mu, iterations);
    }
}
void final_test(){
    /* initialise everything */
    int lattice_size, time_size, iterations = 5000;
    double U, beta = 1, lambda, delta_tau, mu;
    /* 1 */
    lattice_size = 2;
    U = 1;
    // test
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, iterations);
    for(double i = -1; i < 4; i++){
        LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
        generate_lattice(lattice_size, time_size, lattice);
        mu = i * U / 2;
        cout << "mu = " << mu << endl;
        /* sweep across the lattice */
        sweep_lattice(lattice, lattice_size, time_size, U, beta, lambda, delta_tau, mu, iterations);
    }
    /* 2 */
    U = 8;
    // test
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, iterations);
    for(double i = -1; i < 4; i++){
        LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
        generate_lattice(lattice_size, time_size, lattice);
        mu = i * U / 2;
        cout << "mu = " << mu << endl;
        /* sweep across the lattice */
        sweep_lattice(lattice, lattice_size, time_size, U, beta, lambda, delta_tau, mu, iterations);
    }
    /* 3 */
    lattice_size = 5;
    U = 1;
    // test
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, iterations);
    for(double i = 0; i < 5; i++){
        LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
        generate_lattice(lattice_size, time_size, lattice);
        mu = i * U / 2;
        cout << "mu = " << mu << endl;
        /* sweep across the lattice */
        sweep_lattice(lattice, lattice_size, time_size, U, beta, lambda, delta_tau, mu, iterations);
    }
    /* 4 */
    U = 8;
    // test
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, iterations);
    for(double i = 0; i < 5; i++){
        LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
        generate_lattice(lattice_size, time_size, lattice);
        mu = i * U / 2;
        cout << "mu = " << mu << endl;
        /* sweep across the lattice */
        sweep_lattice(lattice, lattice_size, time_size, U, beta, lambda, delta_tau, mu, iterations);
    }
    /* 5 */
    U = 1;
    beta = 8;
    // test
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, iterations);
    for(double i = 0; i < 5; i++){
        LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
        generate_lattice(lattice_size, time_size, lattice);
        mu = i * U / 2;
        cout << "mu = " << mu << endl;
        /* sweep across the lattice */
        sweep_lattice(lattice, lattice_size, time_size, U, beta, lambda, delta_tau, mu, iterations);
    }
}
void final_test_i(){
    /* initialise everything */
    int lattice_size, time_size, iterations = 5000;
    double U, beta = 1, lambda, delta_tau, mu;
    /* 1 */
    lattice_size = 2;
    U = 1;
    beta = 8;
    // test
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, iterations);
    for(double i = -4; i < 11; i++){
        LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
        generate_lattice(lattice_size, time_size, lattice);
        mu = i * U / 4;
        cout << "mu = " << mu << endl;
        /* sweep across the lattice */
        sweep_lattice(lattice, lattice_size, time_size, U, beta, lambda, delta_tau, mu, iterations);
    }
}


/* ------ Main QMC Program ------ */
int main(){
    test_sweep(5, 2000, 1, 8, -1);
    test_sweep(5, 2000, 8, 1, -4);
    test_sweep(5, 2000, 8, 1, -8);
    test_sweep(5, 2000, 1, 1, -0.5);
}
