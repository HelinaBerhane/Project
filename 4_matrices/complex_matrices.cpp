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
void generate_real_array(COMPLEX array[], const int array_length, const int max_rand){
    for(int i = 0; i < array_length; i++){
        array[i].r = basic_random_int(max_rand);	//1 to x
        array[i].i = 0;
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
void scalar_powers(const COMPLEX& number, const int power, COMPLEX& result){
    la::complex<double> laResult = la::complex<double>(number);
    la::complex<double> laNumber = la::complex<double>(number);
    for(int i = 1; i < power; i++){
        laResult *= laNumber;
    }
    result = laResult.toCOMPLEX();
}//working
void scalar_exponential_main(const COMPLEX& number, const int iterations, COMPLEX& result){
    COMPLEX division, total_division;
    result.r = 1;
    result.i = 0 ;
    for(int step = 1; step <= iterations; step++){   //sum (from 1 to n)
        total_division.r = 1;
        total_division.i = 0;
        for(int i = 1; i <= step; i++){        //    ( num^n / n!)
            scalar_division(number, i, division);
            scalar_product(total_division, division);
        }
        scalar_sum(result, total_division);
    }
}//probably working
//void scalar_exponential(const COMPLEX& number, const int iter, COMPLEX& result){
    //COMPLEX power;
    //COMPLEX division;
    //result.r = 0;
    //result.i = 0;
    //for(int i = 0; i < iter; i++){
    //    scalar_powers(number, i, power);
    //    scalar_division(power, factorial(i), division);
    //    scalar_addition(result, division, result);
    //}
//}
//COMPLEX rec_scalar_exp_step(const COMPLEX& number, const int step){
//    COMPLEX result, division, multiplication;
//	if(step <= 1){
//        result.r = 1;
//        return result;
//	}else{
//        scalar_division(number,step,division);
//        scalar_multiplication(division, rec_scalar_exp_step(step-1), multiplication);
//        return multiplication;
//	}
//}
//void recursive_scalar_exponential(const COMPLEX& number, const int iter, COMPLEX& result){
    //COMPLEX power;
    //COMPLEX division;
    //result.r = 0;
    //result.i = 0;
    //for(int i = 0; i < iter; i++){
    //    scalar_powers(number, i, power);
    //    scalar_division(power, factorial(i), division);
    //    scalar_addition(result, division, result);
    //}
//}
// QMC - [/1]

/* array manipulation [1/1] */
//void array_powers(COMPLEX array[], const int len, const int power){/**/
    /*
    for(int i = 0; i < len; i++){
        array[i] = complex_power(array[i], power, result);
    }
    */
//}                       //empty
void vector_exponential(const LaVectorComplex& vector, const int matrix_size, const int iterations, LaVectorComplex& result){
    for(int i = 0; i < matrix_size; i++){
        scalar_exponential_main(vector(i), iterations, result(i));
    }
    //print_vector(result, "vector exponential");
}//working

/* Matrix manipulation [15/15]*/
void matrix_negative(const int matrix_size, LaGenMatComplex& matrix){
    LaGenMatComplex result = LaGenMatComplex::zeros(matrix_size, matrix_size);
    for(int i = 0; i < matrix_size; i++){
        for(int j = 0; j < matrix_size; j++){
            result(i, j).r -= matrix(i, j).r;
            result(i, j).i -= matrix(i, j).i;
        }
    }
    matrix = result.copy();
}//working
void matrix_negative(const int matrix_size, const LaGenMatComplex& matrix, LaGenMatComplex& result){
    result = LaGenMatComplex::zeros(matrix_size, matrix_size);
    for(int i = 0; i < matrix_size; i++){
        for(int j = 0; j < matrix_size; j++){
            result(i, j).r -= matrix(i, j).r;
            result(i, j).i -= matrix(i, j).i;
        }
    }
}//working
void matrix_sum(const int matrix_size, LaGenMatComplex& sum, const LaGenMatComplex& matrix){//to test
    for(int i = 0; i < matrix_size; i++){
        for(int j = 0; j < matrix_size; j++){
            sum(i, j).r += matrix(i, j).r;
            sum(i, j).i += matrix(i, j).i;
        }
    }
}//should be working


void matrix_eigenvstuff(const LaGenMatComplex& matrix, LaVectorComplex& eigenvalues, LaGenMatComplex& eigenvectors){
    //LaEigSolve: http://lapackpp.sourceforge.net/html/laslv_8h.html#086357d17e9cdcaec69ab7db76998769
    LaEigSolve(matrix, eigenvalues, eigenvectors);
}//working
void recombine_diagonalised_matrices(const int matrix_size, LaGenMatComplex& eigenvectors, const LaVectorComplex& eigenvalues, LaGenMatComplex& result){
    /* initialise  everything */
    LaGenMatComplex eigenvalueMatrix = LaGenMatComplex::zeros(matrix_size, matrix_size);
    LaGenMatComplex transposeEigenvectors;
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
}//working
void matrix_inverse(LaGenMatComplex& matrix, int matrix_size){
    // LaLUInverseIP: http://lapackpp.sourceforge.net/html/laslv_8h.html#a042c82c5b818f54e7f000d068f14189
    LaVectorLongInt PIV = LaVectorLongInt(matrix_size);
    LUFactorizeIP(matrix, PIV);
    LaLUInverseIP(matrix, PIV);
}//working
void matrix_exponential(const LaGenMatComplex& matrix, const int matrix_size, const int iterations, LaGenMatComplex& result){
    /* initialise everything */
    LaVectorComplex eigenvalues = LaVectorComplex(matrix_size);
    LaGenMatComplex eigenvectors = LaGenMatComplex::zeros(matrix_size, matrix_size);
    LaGenMatComplex diagonalEigenExp = LaGenMatComplex::zeros(matrix_size, matrix_size);
    LaVectorComplex eigenExponential = LaVectorComplex(matrix_size);
    /* calculate eigenstuff */
    matrix_eigenvstuff(matrix, eigenvalues, eigenvectors);
    //print_matrix(eigenvectors, "eigenvectors");
    //print_vector(eigenvalues, "eigenvalues");
    /* calculate exponentials */
    vector_exponential(eigenvalues, matrix_size, iterations, eigenExponential);
    //print_vector(eigenExponential, "exponential eigenvalues");
    /* multiply them back together to get the matrix */
    recombine_diagonalised_matrices(matrix_size, eigenvectors, eigenExponential, result);
}//should be working
void diagonal_matrix_exponential(const LaGenMatComplex& matrix, const int matrix_size, const int iterations, LaGenMatComplex& result){
    result = LaGenMatComplex::zeros(matrix_size, matrix_size);
    for(int i = 0; i < matrix_size; i++){
        //cout << "matrix_ii "<< matrix(i,i) << endl;
        scalar_exponential_main(matrix(i,i), iterations, result(i,i));
        //cout << "e^matrix_ii "<< result(i,i) << endl;
    }
}//working
void matrix_transpose(const LaGenMatComplex& matrix, const int matrix_size, LaGenMatComplex& result){
    result = LaGenMatComplex::zeros(matrix_size, matrix_size);
    for(int i = 0; i < matrix_size; i++){
        for(int j = 0; j < matrix_size; j++){
            result(i, j) = matrix(j, i);
        }
    }
}//working
void matrix_product(LaGenMatComplex& product, const LaGenMatComplex& matrix){
    LaGenMatComplex result = matrix.copy();
    Blas_Mat_Mat_Mult(product, matrix, result);
    product = result.copy();
}//working
void five_matrix_multiplication(const LaGenMatComplex& matrixA, const LaGenMatComplex& matrixB, const LaGenMatComplex& matrixC, const LaGenMatComplex& matrixD, const LaGenMatComplex& matrixE, LaGenMatComplex& result){
    result = matrixA.copy();
    //print_matrix(result, "A");
    /* AB */
    //print_matrix(matrixB, "B");
    matrix_product(result, matrixB);
    //print_matrix(result, "AB");
    /* ABC */
    //print_matrix(matrixC, "C");
    matrix_product(result, matrixC);
    //print_matrix(result, "ABC");
    /* ABCD */
    //print_matrix(matrixD, "D");
    matrix_product(result, matrixD);
    //print_matrix(result, "ABCD");
    /* ABCDE */
    //print_matrix(matrixE, "E");
    matrix_product(result, matrixE);
    //print_matrix(result, "ABCDE");
}//working
COMPLEX simple_matrix_determinant(const LaGenMatComplex& matrix){
    /* initialise everything */
    COMPLEX A;
    COMPLEX B;
    /* multiply opposite corners */
    scalar_multiplication(matrix(0,0), matrix(1,1), A);
    scalar_multiplication(matrix(0,1), matrix(1,0), B);
    /* - B */
    B.r = -B.r;
    B.i = -B.i;
    /* calculate determinant */
    scalar_sum(A, B);
    return A;
}//working
COMPLEX determinant_coefficient(const LaGenMatComplex& matrix, const int element){
    COMPLEX coefficient;
    if(element % 2 == 1){
        // if odd
        coefficient.r = - matrix(0, element).r;
        coefficient.i = - matrix(0, element).i;
    }else{
        // if even
        coefficient.r = matrix(0, element).r;
        coefficient.i = matrix(0, element).i;
    }
    return coefficient;
}//working
COMPLEX my_matrix_determinant(const int matrix_size, const LaGenMatComplex& matrix){
    /* initialise everything */
    COMPLEX determinant;
    determinant.r = 0;
    determinant.i = 0;
    /* do stuff */
    if(matrix_size == 2){
        /* calculate the determinant */
        return simple_matrix_determinant(matrix);
    }else{
        for(int element = 0; element < matrix_size; element++){//for each element in the first row
            /* initialise everything */
            LaGenMatComplex cofactorMatrix;
            COMPLEX coefficient;
            int cofactor_size = matrix_size - 1;
            /* determine the coefficient */
            coefficient = determinant_coefficient(matrix, element); // = +- the element
            /* calculate the cofactor */
            cofactorMatrix = LaGenMatComplex::zeros(cofactor_size, cofactor_size);
            generate_cofactor_matrix(matrix_size, matrix, element, cofactorMatrix);
            //print_matrix(cofactorMatrix, "cofactorMatrix");
            /* finish calculation */
            scalar_sum(determinant, scalar_multiple(coefficient, my_matrix_determinant(cofactor_size, cofactorMatrix)));
        }
    }
    return determinant;
}//working
void matrix_determinant(const int matrix_size, const LaGenMatComplex& matrix, COMPLEX& result){
    /* initialise everything */
    LaVectorComplex eigenvalues = LaVectorComplex(matrix_size);
    LaGenMatComplex eigenvectors = LaGenMatComplex::zeros(matrix_size, matrix_size);
    result.r = 1;
    result.i = 0;
    /* calculate eigenvectors */
    matrix_eigenvstuff(matrix, eigenvalues, eigenvectors);
    /* calculate determinant */
    for(int i = 0; i < matrix_size; i++){
        scalar_product(result, eigenvalues(i));
    }
}//working

// QMC [7/8]
float lambda_calculation(const float U){
    return acoshf(exp(sqrt(0.125*U)/2));
}//working
float delta_tau_calculation(const float U){
    return sqrt(0.125 / U);
}//working
void initial_parameter_calculation(const float U, const float beta, float& lambda, float& delta_tau, int& time_size){
    // delta_tau = sqrt(0.125 / U);            // by convension
    lambda = acoshf(exp(sqrt(0.125*U)/2));   // by definition
    time_size = ceil(beta / lambda);      // by definition
    delta_tau = beta / time_size;
}
void print_initial_parameters(float U, float beta, float lambda, float delta_tau, int time_size, int lattice_size){
	cout << "no of lattice points = " << lattice_size << endl;
	cout << "no of time slices = " << time_size << endl;
	cout << "U = " << U << endl;
	cout << "beta = " << beta << endl;
	cout << "lambda = " << lambda << endl;
	cout << "delta tau = " << delta_tau << endl;
}


void V_calculation(const COMPLEX slice[], const int time_size, const float U, const float lambda, const float sigma, const float delta_tau, LaGenMatComplex& V){//should be working
    /* initialise everything */
    COMPLEX elements[time_size];
    float mu = 0;
    /* V = lambda sigma s_l + + mu - U / 2 */
    for(int i = 0; i < time_size; i++){
        scalar_multiplication(slice[i], lambda * sigma / delta_tau, elements[i]);
        elements[i].r = elements[i].r + mu - U / 2;
    }
    /* given a lattice */
    array_to_diag(elements, time_size, V);
}
void B_calculation(LaGenMatComplex& H, LaGenMatComplex& V, LaGenMatComplex& B, const int matrix_size, const int iterations){//should be working
    //B = exp(-H)exp(-V)
    /* initialise everything */
    LaGenMatComplex negH;
    LaGenMatComplex negV;
    LaGenMatComplex expH;
    LaGenMatComplex expV;
    /* negate matrices (not in place) */
    matrix_negative(matrix_size, H, negH);
    matrix_negative(matrix_size, V, negV);
    /* calculate exponentials */
    matrix_exponential(negH, matrix_size, iterations, expH);
    matrix_exponential(negV, matrix_size, iterations, expV);
    /* print exponential matrices */
    //print_matrix(expH, "e^-H");
    //print_matrix(expV, "e^-V");
    /* multiply exponentials */
    B = expH.copy();
    matrix_product(B, expV);
    /* print result */
    //print_matrix(B, "B");
}
void O_calculation(const int matrix_size, const LaGenMatComplex& BA, const LaGenMatComplex& BB, const LaGenMatComplex& BC, const LaGenMatComplex& BD, const LaGenMatComplex&BE, LaGenMatComplex& O){//should be working
    //O = 1 + B(m) B(m-1) B(...) B(1)
    /* initialise everything */
    LaGenMatComplex I = LaGenMatComplex::eye(matrix_size, matrix_size);
    //LaGenMatComplex multiplication;
    /* multiply exponentials */
    five_matrix_multiplication(BA, BB, BC, BD, BE, O);
    /* add I */
    matrix_sum(matrix_size, O, I);
}
void detO_calculation(const int matrix_size, const LaGenMatComplex& O, COMPLEX& detO){//to test
    /* calculate det O */
    detO = my_matrix_determinant(matrix_size, O);
}
void calculate_weight(const int matrix_size, const COMPLEX latticeUP[], const float U, const float lambda, const float delta_tau, COMPLEX& weight){//to test

    int lattice_size = matrix_size, time_size = matrix_size, iterations = 1000;
    COMPLEX latticeDOWN[matrix_size];
    LaGenMatComplex H;
    LaGenMatComplex VUP = LaGenMatComplex::zeros(matrix_size, matrix_size);
    LaGenMatComplex VDOWN = LaGenMatComplex::zeros(matrix_size, matrix_size);
    LaGenMatComplex proBUP = LaGenMatComplex::eye(matrix_size, matrix_size);
    LaGenMatComplex proBDOWN = LaGenMatComplex::eye(matrix_size, matrix_size);
    LaGenMatComplex I = LaGenMatComplex::eye(matrix_size, matrix_size);
    LaGenMatComplex BUP = LaGenMatComplex::zeros(matrix_size, matrix_size);
    LaGenMatComplex BDOWN = LaGenMatComplex::zeros(matrix_size, matrix_size);
    LaGenMatComplex OUP = LaGenMatComplex::eye(matrix_size, matrix_size);
    LaGenMatComplex ODOWN = LaGenMatComplex::eye(matrix_size, matrix_size);
    COMPLEX detOUP;
    COMPLEX detODOWN;

    /* generate lattices */
    for(int l = 0; l < matrix_size; l++){
        copy_negative_scalar(latticeUP[l], latticeDOWN[l]);
    }
    /* generate H */
    generate_H(matrix_size, H);

    /* generate V matrices */
    V_calculation(latticeUP, time_size, U, lambda, 1, delta_tau, VUP);
    V_calculation(latticeDOWN, time_size, U, lambda, 1, delta_tau, VDOWN);

    /* multiply B matrices */
    for(int t = time_size - 1; t >= 0 ; t--){
        /*   for each time slice   */

        /* calculate B(t) matrices */
        B_calculation(H, VUP, BUP, lattice_size, iterations);
        B_calculation(H, VDOWN, BDOWN, lattice_size, iterations);

        /* multiply the matrices */
        matrix_product(proBUP, BUP);
        matrix_product(proBDOWN, BUP);
    }

    /* calculate O matrices */
    matrix_sum(matrix_size, OUP, I);
    matrix_sum(matrix_size, ODOWN, I);

    /* calculate det(O)s */
    matrix_determinant(matrix_size, OUP, detOUP);
    matrix_determinant(matrix_size, ODOWN, detODOWN);

    /* calculate the weight */
    scalar_multiplication(detOUP, detODOWN, weight);

}



void sweep_lattice(const int matrix_size, LaGenMatComplex& lattice, const float U, const int iterations){//in progress
    /* Plan */

        /* Input */
            // matrix_size      - int
            // lattice          - LaGenMatComplex&
            // U                - float
            // iterations       - int

        /* Processing */
            // Calculate initial parameters
                // Calculate lambda
                // Calculate delta_tau
            // for each iteration
                // for each time slice
                    // isolate the spins in an array
                    // for each lattice point
                        // calculate the probability of the spin flipping
                        // decide whether it flips or not
                        // record the flip in the original matrix
                        // record the measurements

        /* Output */
            // a pritout of the lattice over time?
            // probabiliy of flipping at each stage
            // average spin
            // ... ?

    /* initialise everything */
        // int lattice_size = matrix_size, time_size = matrix_size;
    COMPLEX weightBefore, weightAfter, slice[matrix_size];
    float probability, prob;
    float lambda = lambda_calculation(U), delta_tau = delta_tau_calculation(U);

    string result;
    int count = 0;

    /* set up output headings */
    cout.width(11);
    cout << "weight";
    cout << " lattice" << endl;

    /* sweep through the lattice */
    for(int i = 0; i < iterations; i++){
        for(int time_slice = 0; time_slice < matrix_size; time_slice++){

            /* isolate the time slice as a COMPLEX array */
            isolate_row(lattice, matrix_size, time_slice, slice);

            for(int lattice_site = 0; lattice_site < matrix_size; lattice_site++){

                /* calculate the weight before the flip */
                calculate_weight(matrix_size, slice, U, lambda, delta_tau, weightBefore);

                /* propose the flip */
                flip_scalar(slice[lattice_site]);

                /* calculate the weight after the flip */
                calculate_weight(matrix_size, slice, U, lambda, delta_tau, weightAfter);

                /* calculate the ratio of weights */
                probability = weightAfter.r / weightBefore.r;

                /* accept or reject the flip */
                    // cout << "prob: " << probability << endl;
                if(abs(probability) >= 1){
                    flip_scalar(lattice(time_slice, lattice_site)); //accept
                    result = "accepted";
                }else{
                    prob = random_probability();
                    if(probability > prob){
                        flip_scalar(lattice(time_slice, lattice_site)); //accept
                        result = "accepted";
                    }else{
                        // cout << " rejected" << endl;
                        result = "rejected";
                    }
                }
                /* comments */
                    //for negative values, we do some integration
                    //P\to\tilde{P} = |P| and  F\to \tilde
                    //you have to multiply each quan you measure bu the sign
                count++;
                cout << " (" << count <<") " << result << " - " << probability;
                cout.width(15);
                cout << " - weightBefore: " << weightBefore << ", weightAfter: " << weightAfter << endl;
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
}


/* Testing [20/21] */
void test_random_int(){
    int max_rand = 9, iterations = 1000000000, test;
    int zero = 0, one = 0, two = 0, three = 0, four = 0;
    int five = 0, six = 0, seven = 0, eight = 0, nine = 0;
    for(int i = 0; i < iterations; i++){
        test = basic_random_int(max_rand);
        switch(test){
            case 0: zero++;
                break;
            case 1: one++;
                break;
            case 2: two++;
                break;
            case 3: three++;
                break;
            case 4: four++;
                break;
            case 5: five++;
                break;
            case 6: six++;
                break;
            case 7: seven++;
                break;
            case 8: eight++;
                break;
            case 9: nine++;
                break;
        }
    }
    cout << "0: " << zero << endl;
    cout << "1: " << one << endl;
    cout << "2: " << two << endl;
    cout << "3: " << three << endl;
    cout << "4: " << four << endl;
    cout << "5: " << five << endl;
    cout << "6: " << six << endl;
    cout << "7: " << seven << endl;
    cout << "8: " << eight << endl;
    cout << "9: " << nine << endl;
}//working
void test_random_float(){
    int count = 10, min = 0, max = 10;
    for (int i = 0; i < count; i++) {
        cout << random_float(min, max) << endl;
    }
}//working
void test_negative_scalar(){
    COMPLEX scalar;
    COMPLEX negativeScalar;
    for(int i = 0; i < 10; i ++){
        generate_scalar(scalar, 9);
        copy_negative_scalar(scalar, negativeScalar);
        print_scalar(negativeScalar, "negative scalar");
    }
}//working
void test_scalar_manipulation(const int max_rand){
    COMPLEX compA;
    generate_scalar(compA, max_rand);
    COMPLEX compB;
    generate_scalar(compB, max_rand);
    int realB = 0;
    generate_scalar(realB, max_rand);
    COMPLEX result;

    cout << "compA = " << compA << endl;
    cout << "compB = " << compB << endl;
    cout << "realB = " << realB << endl;
    cout << endl;

    for(int i = 1; i < 5; i++){                                                 //factorials
        cout << "factorial(" << i << "): " << factorial(i) << endl;
    }

    scalar_addition(compA, compB, result);                                      //addition/subtraction
    cout << "scalar addition: " << result << endl << endl;

    scalar_multiplication(compA, realB, result);                                //multiplication
    cout << "scalar multiplication by scalar: " << result << endl << endl;

    scalar_multiplication(compA, compB, result);
    cout << "scalar multiplication by complex: " << result << endl << endl;

    scalar_division(compA, realB, result);         //division
    cout << "scalar division by scalar: " << result << endl << endl;

    scalar_division(compA, compB, result);
    cout << "scalar division by complex: " << result << endl << endl;

    for(int i = 1; i < 5; i++){
        scalar_powers(compA, i, result);
        cout << "scalar powers - A^" << i << " = " << result << endl;
    }
    cout << endl;

    COMPLEX sum;
    sum.r = 0;
    sum.i = 0;
    for(int i = 0; i < 5; i++){
        cout << "sum(" << i << ") = " << sum << endl;
        scalar_sum(sum, compA);
    }
    cout << "sum = " << sum << endl;

}//working
void test_eigenvalues(const int matrix_size, const int max_rand){
    /* initialise everything */
    LaGenMatComplex matrix;
    LaVectorComplex eigenvalues = LaVectorComplex(matrix_size);
    LaGenMatComplex eigenvectors = LaGenMatComplex::zeros(matrix_size, matrix_size);
    LaGenMatComplex result;
    /* generate matrix */
    generate_matrix(matrix_size, max_rand, matrix);
    print_matrix(matrix, "initial matrix");
    /* calculate eigenstuff */
    matrix_eigenvstuff(matrix, eigenvalues, eigenvectors);
    /* multiply them back together to get the matrix */
    recombine_diagonalised_matrices(matrix_size, eigenvectors, eigenvalues, result);
    /* wolfram test */
        // 2x2 real:
        // 3x3 complex: {{1+7i, 1+3i, 5+7i},{7i, 6+i, 5+4i},{5+7i, 5+4i, 6}}
}//working
void test_inverse(const LaGenMatComplex& initialMatrix, const int matrix_size){
    LaGenMatComplex inverseMatrix;
    inverseMatrix = initialMatrix.copy();
    matrix_inverse(inverseMatrix, matrix_size);
    print_matrix(inverseMatrix, "inverse matrix");
}//working
void test_scalar_sum(const int max_rand, const int iterations){
    COMPLEX number, step;
    generate_scalar(number, max_rand);
    step = number;
    for(int i = 0; i < iterations; i++){
        cout << step << endl;
        scalar_sum(number, step);
    }
    cout << number << endl;
}//working
void test_scalar_exponential(const int max_rand, const int iterations){
    COMPLEX number, result;
    generate_scalar(number, max_rand);
    cout << endl << "scalar exponential test no.: " << number << endl << endl;
    scalar_exponential_main(number, iterations, result);
    cout << "e^" << number << " = " << result << endl;
}//working
void test_scalar_exponential(COMPLEX& number, const int iterations, COMPLEX& result){
    scalar_exponential_main(number, iterations, result);
    cout << "e^" << number << " = " << result << endl;
}//working
void test_matrix_subtraction(const int matrix_size, const int max_rand){
    LaGenMatComplex matrix;
    LaGenMatComplex result = LaGenMatComplex::zeros(matrix_size, matrix_size);
    generate_matrix(matrix_size, max_rand, matrix);
    print_matrix(matrix, "Matrix");
    matrix_negative(matrix_size, matrix, result);
    print_matrix(result, "- Matrix");
    matrix_negative(matrix_size, matrix);
    print_matrix(matrix, "- Matrix (in place)");
}
void test_matrix_multiplication(const int matrix_size, const int max_rand){
    int matrix_volume = matrix_size * matrix_size;
    /* generate matrix A */
    COMPLEX elements[matrix_volume];
    generate_array(elements, matrix_volume, max_rand);
	LaGenMatComplex matrixA = LaGenMatComplex(elements, matrix_size, matrix_size, false );
    print_matrix(matrixA, "Matrix A");
    /* generate matrix B */
    generate_array(elements, matrix_volume, max_rand);
	LaGenMatComplex matrixB = LaGenMatComplex(elements, matrix_size, matrix_size, false );
    print_matrix(matrixB, "Matrix B");
    /* generate matrix A^T */
    LaGenMatComplex transposeA;
    matrix_transpose(matrixA, matrix_size, transposeA);
    print_matrix(transposeA, "transpose A");
    /* generate matrix B^T */
    LaGenMatComplex transposeB;
    matrix_transpose(matrixB, matrix_size, transposeB);
    print_matrix(transposeB, "transpose B");
    /* initialise result */
    LaGenMatComplex result = LaGenMatComplex::zeros(matrix_size, matrix_size);
    /* A * B */
    Blas_Mat_Mat_Mult(matrixA, matrixB, result);
    print_matrix(result, "Matrix A * Matrix B");
    /* A^T * B */
    Blas_Mat_Mat_Mult(transposeA, matrixB, result);
    print_matrix(result, "Matrix A^T * Matrix B");
    /* A * B^T */
    Blas_Mat_Mat_Mult(matrixA, transposeB, result);
    print_matrix(result, "Matrix A * Matrix B^T");
    /* A^T * B^T */
    Blas_Mat_Mat_Mult(transposeA, transposeB, result);
    print_matrix(result, "Matrix A^T * Matrix B^T");
}//working
void test_matrix_product(const int matrix_size, const int max_rand){
    /* initialise everything */
    LaGenMatComplex matrixA;
    LaGenMatComplex matrixB;
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




// void test_matrix_arrays(){
//
//     /* initialise everything */
//     int n = 3, matrix_size = 5, max_rand = 9;
//     LaGenMatComplex* matrices[1];
//     *matrices[0] = LaGenMatComplex::eye(2, 2);
//
//     LaGenMatComplex mgds;
//
//     /* generate everything */
//     // generate_matrix(matrix_size, max_rand, matrices[0]);
//     // *matrices[0] = LaGenMatComplex::eye(matrix_size, matrix_size);
//
//     mgds = *matrices[0];
//     cout << matrices[0] << endl;
//     cout << mgds(1,1) << endl;
//         // this is an array of pointers to matrices
//     // LaGenMatComplex product = LaGenMatComplex::eye(2, 2);
//     //
//     // for(int i = 0; i < n; i++){
//     //
//     //     /* generate everything */
//     //     generate_matrix(matrix_size, max_rand, *matrices[n]);
//     //
//     //     /* print everything */
//     //     cout << "(" << n << ")" << endl;
//     //     print_matrix(*matrices[n]);
//     // }
// }



void test_five_matrix_multiplication(const int matrix_size, const int max_rand){
    /* initialise everything */
    LaGenMatComplex matrixA;
    LaGenMatComplex matrixB;
    LaGenMatComplex matrixC;
    LaGenMatComplex matrixD;
    LaGenMatComplex matrixE;
    LaGenMatComplex result;
    /* generate everything */
    generate_matrix(matrix_size, max_rand, matrixA);
    generate_matrix(matrix_size, max_rand, matrixB);
    generate_matrix(matrix_size, max_rand, matrixC);
    generate_matrix(matrix_size, max_rand, matrixD);
    generate_matrix(matrix_size, max_rand, matrixE);
    /* ABCDE */
    five_matrix_multiplication(matrixA, matrixB, matrixC, matrixD, matrixE, result);
    print_matrix(result, "ABCDE");
    //{{1+7i, 5+7i},{7i, 1+3i}}*{{6+i, 5+7i},{5+4i, 5+4i}}*{{6, 8+8i},{7+i, 6+6i}}*{{8+8i, 1+i},{8+4i, 5}}*{{3, 1+7i},{5+3i, 4+7i}}
}//working
void test_matrix_exponential(const int matrix_size, const int max_rand, const int iterations){
    /* initialise everything */
    LaGenMatComplex matrix;
    LaGenMatComplex result;
    result = LaGenMatComplex::zeros(matrix_size, matrix_size);
    /* generate matrix */
    generate_matrix(matrix_size, max_rand, matrix);
    print_matrix(matrix, "initial matrix");
    /* calculate exponential */
    matrix_exponential(matrix, matrix_size, iterations, result);
    print_matrix(result, "e^(matrix)");
}//working
void test_idenpotent_exponential(const int iterations){
    /* generate the matrix */
    int numbers [] = {2, -2, -4, -1, 3, 4, 1, -2, -3};
    COMPLEX elements[9];
    for(int i = 0; i < 9; i++){
        elements[i].r = numbers[i];
        elements[i].i = 0;
    }
    LaGenMatComplex matrix = LaGenMatComplex(elements, 3, 3, false );
    LaGenMatComplex result = LaGenMatComplex::zeros(3, 3);
    //print_matrix(matrix, "initial matrix");
    /* calculate the exponential */
    for(int j = 1; j <= 5; j++){
        matrix_exponential(matrix, 3, j, result);
        cout << j << " iterations:" << endl;
        print_matrix(result);
    }
    matrix_exponential(matrix, 3, iterations, result);
    print_matrix(result, "idenpotent exponential");
}//working
void test_diagonal_exponential(const int iterations){
    LaGenMatComplex I = LaGenMatComplex::eye(3, 3);
    LaGenMatComplex result = LaGenMatComplex::zeros(3, 3);
    diagonal_matrix_exponential(I, 3, iterations, result);
    print_matrix(I, "I");
    print_matrix(result);
}
void test_simple_matrix_determinant(const int max_rand){
    /* initialise everything */
    LaGenMatComplex matrix;
    /* generate matrix */
    generate_matrix(2, max_rand, matrix);
    print_matrix(matrix, "matrix");
    /* calculate determinant */
    simple_matrix_determinant(matrix);
    print_scalar(simple_matrix_determinant(matrix), "det(M)");
}//working
void test_determinant_coefficient(){
    /* initialise everything */
    LaGenMatComplex matrix;
    int matrix_size = 4, max_rand = 9;
    /* generate matrix */
    generate_matrix(matrix_size, max_rand, matrix);
    print_matrix(matrix, "matrix");
    /* calculate coefficients */
    for(int element = 0; element < matrix_size; element++){
        cout << determinant_coefficient(matrix, element) << " ";
    }
    cout << endl;
}//working
void test_reduced_matrix(){
    /* initialise everything */
    int matrix_size = 4, max_rand = 9;
    LaGenMatComplex matrix;
    LaGenMatComplex cofactorMatrix = LaGenMatComplex::zeros(matrix_size - 1, matrix_size - 1);
    /* generate matrix */
    generate_matrix(matrix_size, max_rand, matrix);
    print_matrix(matrix, "matrix");
    /* calculate reduced matrix */
    for(int element = 0; element < matrix_size; element++){
        generate_cofactor_matrix(matrix_size, matrix, element, cofactorMatrix);
        print_matrix(cofactorMatrix, "cofactorMatrix");
    }
}//working
void test_matrix_determinant(){//in progress
    /* initialise everything */
    int matrix_size = 4, max_rand = 9;
    LaGenMatComplex matrix;
    COMPLEX result;
    /* generate matrix */
    generate_matrix(matrix_size, max_rand, matrix);
    print_matrix(matrix, "initial matrix");
    /* calculate determinant */
    print_scalar(my_matrix_determinant(matrix_size, matrix), "my determinant");
    matrix_determinant(matrix_size, matrix, result);
    print_scalar(result, "eigenvalue determinant");
}//working
void test_isolate_row(){

    /* initialise everything  */
    int matrix_size = 4, max_rand = 9;
    LaGenMatComplex matrix;
    COMPLEX array[matrix_size];

    /* generate matrix */
    generate_matrix(matrix_size, max_rand, matrix);
    print_matrix(matrix, "matrix");

    /* isolate and print rows */
    for(int row = 0; row < matrix_size; row++){
        isolate_row(matrix, matrix_size, row, array);
        cout << " Row (" << row << "):";
        print_array(array, matrix_size);
    }
}//working

// QMC [8/10]
void test_random_probability(){
    int count = 10;
    for (int i = 0; i < count; i++) {
        cout << random_probability() << endl;
    }
}//working
void test_lattice_generation(){

    /* initialise everything */
    int matrix_size = 5;
    LaGenMatComplex lattice;

    /* generate the lattice */
    generate_lattice_matrix(matrix_size, lattice);

    /* print result */
    print_matrix(lattice, "random spins");

}//working
void test_parameter_calculation(){
    float U = 1, lambda = lambda_calculation(U), delta_tau = delta_tau_calculation(U);
    cout << "U = " << U << endl;
    cout << "l = " << lambda << endl;
    cout << "t = " << delta_tau << endl;
}//working
void test_H(const int matrix_size){
    /* initialise everything */
    LaGenMatComplex H;
    LaVectorComplex eigenvalues = LaVectorComplex(matrix_size);
    LaGenMatComplex eigenvectors = LaGenMatComplex::zeros(matrix_size, matrix_size);
    /* generate matrix */
    generate_H(matrix_size, H);
    print_matrix(H);
    /* calculate eigenstuff */
    matrix_eigenvstuff(H, eigenvalues, eigenvectors);
    print_vector(eigenvalues, "eigenvalues");
    // eigenvalues are 2 cos(n pi / q), where q = the matrix size
}//working
void test_V_generation(){//should work

    /* initialise everything */
        // float delta tau = ?, lamba = 1, s_i(l) = \pm 1, mu = 0?
    int matrix_size = 5;
    LaGenMatComplex V = LaGenMatComplex::zeros(matrix_size, matrix_size);
    COMPLEX time_slice[matrix_size];
    float U = 1, lambda = lambda_calculation(U), delta_tau = delta_tau_calculation(U);

    /* generate the lattice */
    generate_lattice_array(matrix_size, time_slice);
    V_calculation(time_slice, matrix_size, U, lambda, 1, delta_tau, V);

    /* print result */
    print_matrix(V);
}
void test_B_generation(){//should work
    /* initialise everything */
    int time_size = 5, max_rand = 9, iterations = 1000;
    LaGenMatComplex H = LaGenMatComplex::eye(time_size, time_size);
    LaGenMatComplex V = LaGenMatComplex::eye(time_size, time_size);
    LaGenMatComplex B = LaGenMatComplex::zeros(time_size, time_size);
    /* generate matrices */
    for(int i = 0; i < time_size; i++){
        H(i,i).r = basic_random_int(max_rand);
        V(i,i).r = basic_random_int(max_rand);
    }
    /* print matrices */
    print_matrix(H, "H");
    print_matrix(V, "V");
    /* calculate B */
    B_calculation(H, V, B, time_size, iterations);
    /* print result */
    print_matrix(B,"B = e^-H e^-V");
}
void test_O_generation(const int time_size, const int iterations){//should work
    /* initialise everything */
    COMPLEX elements[time_size];
    LaGenMatComplex H;
    LaGenMatComplex V = LaGenMatComplex::zeros(time_size, time_size);
    LaGenMatComplex BA = LaGenMatComplex::zeros(time_size, time_size);
    LaGenMatComplex BB = LaGenMatComplex::zeros(time_size, time_size);
    LaGenMatComplex BC = LaGenMatComplex::zeros(time_size, time_size);
    LaGenMatComplex BD = LaGenMatComplex::zeros(time_size, time_size);
    LaGenMatComplex BE = LaGenMatComplex::zeros(time_size, time_size);
    LaGenMatComplex O = LaGenMatComplex::zeros(time_size, time_size);
    float U = 1, lambda = lambda_calculation(U), delta_tau = delta_tau_calculation(U);

    /* generate matrices */
    generate_H(time_size, H);
    for(int i = 0; i < time_size; i++){
        /* generate matrices */
        generate_lattice_array(time_size, elements);
        V_calculation(elements, time_size, U, lambda, 1, delta_tau, V);
        /* calculate B */
        if(i == 0){
            B_calculation(H, V, BA, time_size, iterations);
        }else if(i == 1){
            B_calculation(H, V, BB, time_size, iterations);
        }else if(i == 2){
            B_calculation(H, V, BC, time_size, iterations);
        }else if(i == 3){
            B_calculation(H, V, BD, time_size, iterations);
        }else if(i == 4){
            B_calculation(H, V, BE, time_size, iterations);
        }
    }
    O_calculation(time_size, BA, BB, BC, BD, BE, O);
    /* print result */
    print_matrix(BA, "BA");
    print_matrix(BB, "BB");
    print_matrix(BC, "BC");
    print_matrix(BD, "BD");
    print_matrix(BE, "BE");
    print_matrix(O, "O");
}
void test_detO(){//in progress
    //
}
void test_weight(){//working

    /* Plan */
        /* Input */
            // matrix_size  - int
            // slice        - COMPLEX[]
            // U            - float

        /* Processing */
            // calculate initial parameters
                // lambda    - float
                // delta_tau - float
            // calculate the weight

        /* Output */
            // the weight

    /* initialise everything */
    int matrix_size = 5, matrix_volume = matrix_size * matrix_size;
    float U = 1, lambda = lambda_calculation(U), delta_tau = delta_tau_calculation(U);
    COMPLEX lattice[matrix_volume];
    COMPLEX weight;

    /* generate the lattice */
    generate_lattice_array(matrix_volume, lattice);
    print_array(lattice, matrix_size, "lattice");

    /* calculate the weight */
    calculate_weight(matrix_size, lattice, U, lambda, delta_tau, weight);

    /* output the weight */
    print_scalar(weight, "weight");
}
void test_sweep(){
    /* Plan */
        /* [x] Input */
            // [x] matrix_size  - int
            // [x] iterations   - int
            // [x] lattice      - LaGenMatComplex
            // [x] U            - float

        /* [x] Processing */
            // [x] generate a lattice of spins
            // [x] sweep the lattice

        /* [ ] Output */
            // [ ] average spins
            // [ ] acceptance probabilities

    /* initialise everything */
    int matrix_size = 5, iterations = 3;
    LaGenMatComplex lattice;
    float U = 1;
    float lambda = lambda_calculation(U), delta_tau = delta_tau_calculation(U);

    /* print initial parameters */
    cout.width(11);
    cout << "U = " << U << endl;
    cout.width(11);
    cout << "lambda = " << lambda << endl;
    cout.width(11);
    cout << "delta tau = " << delta_tau<< endl;

    /* generate a lattice of spins */
    generate_lattice_matrix(matrix_size, lattice);

    /* sweep the lattice */
    sweep_lattice(matrix_size, lattice, U, iterations);
}
void test_increasing_U(){//in progress

    /* Plan */

        /* [x] Input */
            // [x] matrix_size  - int
            // [x] iterations   - int
            // [x] lattice      - LaGenMatComplex
            // [x] U            - float

        /* [x] Processing */
            // [x] for n increasing values of U
                // [x] generate a lattice of spins
                // [x] sweep the lattice

        /* [ ] Output */
            // [ ] acceptance probabilities

    /* initialise everything */
    int matrix_size = 5, iterations = 1;
    LaGenMatComplex lattice;
    float U, lambda, delta_tau;

    /* test U = 0 to 1 */
    for(int i = 0; i <= 10; i++){

        /* calculate initial parameters */
        U = 2 + 0.8*i;
        lambda = lambda_calculation(U);
        delta_tau = delta_tau_calculation(U);

        /* print initial parameters */
        cout.width(11);
        cout << "U = " << U << endl;
        cout.width(11);
        cout << "lambda = " << lambda << endl;
        cout.width(11);
        cout << "delta tau = " << delta_tau<< endl;

        /* generate a lattice of spins */
        generate_lattice_matrix(matrix_size, lattice);

        /* sweep the lattice */
        sweep_lattice(matrix_size, lattice, U, iterations);
    }
}

/* Generalisation */
/* - Working - */
void n_matrix_product(const COMPLEX storage[], const int matrix_size, const int n, LaGenMatComplex& result){
    /* initialise everything */
    LaGenMatComplex matrix;
    int matrix_volume = matrix_size * matrix_size;
    /* reset variables */
    result = LaGenMatComplex::eye(matrix_size, matrix_size);
    //for each matrix
    for(int m = 0; m < n; m++){
        // reset variables
        matrix = LaGenMatComplex::eye(matrix_size, matrix_size);
        // convert the storage to a matrix
        for(int r = 0; r < matrix_size; r++){
            for(int c = 0; c < matrix_size; c++){
                int e = r * matrix_size + c;
                int i = m * matrix_volume + e;
                matrix(r, c).r = storage[i].r;
                matrix(r, c).i = storage[i].i;
            }
        }
        // multiply with the result
        matrix_product(result, matrix);
        // test
        // print_matrix(matrix, "current matrix");
        // print_matrix(result, "current product");
    }
}
void test_n_matrix_product(){

    /* initialise everything */
    int n = 4, matrix_size = 3, max_rand = 5;
    int storage_size = matrix_size * matrix_size * n;
    COMPLEX storage[storage_size];
    LaGenMatComplex result = LaGenMatComplex::eye(matrix_size, matrix_size);

    /* generate matrices (skip to storage) */
    generate_real_array(storage, storage_size, max_rand);
    print_array(storage, storage_size, "storage");

    /* multiply everything */
    n_matrix_product(storage, matrix_size, n, result);

    print_matrix(result, "result");
}

/* - Testing - */
void general_weight(const int lattice_size, const int time_size, const LaGenMatComplex& lattice, const float U, const float lambda, const float delta_tau, COMPLEX& weight){
    /* Plan */
        /* Input */
            // a lattice        - LaGenMatComplex
            // matrix_size      - int
            // no of matrices   - int
        /* Processing */
            // for sigma = +-1
                // for each time_slice
                    // calculate the V matrix
                    // calculate the B matrix
                    // add the elements of the B matrix to an array
                    // multiply all B matrices together
                // add 1
                // calculate det O
            // calculate det O up * det O down
    /* initialise everything */

    /* initialise everything */
    int iterations = 100;
    int lattice_volume = lattice_size * lattice_size;
    int storage_size = lattice_volume * time_size;
    float sigma;
    LaGenMatComplex H;
    LaGenMatComplex V;
    LaGenMatComplex B;
    LaGenMatComplex O;
    LaGenMatComplex I = LaGenMatComplex::eye(lattice_size, lattice_size);
    COMPLEX product;
    COMPLEX detO;
    COMPLEX slice;
    COMPLEX storage[storage_size];
    for(int i = 0; i < storage_size; i++){
        storage[i].r = 0;
        storage[i].i = 0;
    }
    print_array(storage, storage_size, "storage");

    /* generate H */
    generate_H(lattice_size, H);

    // for sigma = +-1
    for(int s = 0; s < 1; s++){
        // reset all variables
        O = LaGenMatComplex::zeros(lattice_size, lattice_size);
        sigma = s * 2 - 1;
        cout << "sigma = " << sigma << endl;
        // for each time_slice
        for(int t = 0; t < time_size; t++){    // check the order of multiplication!!!
            cout << "current time slice = " << t << endl;
            // reset all variables
            V = LaGenMatComplex::zeros(lattice_size, lattice_size);
            B = LaGenMatComplex::zeros(lattice_size, lattice_size);
            // isolate the time slice
            isolate_row(lattice, lattice_size, t, slice);
            // calculate all variables
            V_calculation(slice, time_size, U, lambda, sigma, delta_tau, V);
            B_calculation(H, V, B, lattice_size, iterations);
            // store the elements of the B matrix in an array
            print_matrix(B, "B");
            for(int r = 0; r < lattice_size; r++){
                for(int c = 0; c < lattice_size; c++){
                    int i = (t * lattice_volume) + (r * lattice_size) + c;
                    storage[i].r = B(r, c).r;
                    storage[i].i = B(r, c).i;
                    // print_scalar(storage[i]);
                }
            }
            cout << endl;
            print_array(storage, storage_size, "storage");
            /* multiply the matrices */
            n_matrix_product(storage, lattice_size, time_size, O);
        }
        print_array(storage, storage_size, "storage");
        // add 1
        matrix_sum(lattice_size, O, I);
        // calculate det O
        matrix_determinant(lattice_size, O, detO);
        // calculate det O up * det O down
        scalar_product(product, detO);
    }
}
void general_sweep(const int lattice_size, LaGenMatComplex& lattice, const float U, const float beta, const int iterations){
    /* Plan */

        /* Input */
            // lattice_size      - int
            // lattice          - LaGenMatComplex&
            // U                - float
            // iterations       - int

        /* Processing */
            // Calculate initial parameters
                // Calculate lambda
                // Calculate delta_tau
            // for each iteration
                // for each time slice
                    // isolate the spins in an array
                    // for each lattice point
                        // calculate the probability of the spin flipping
                        // decide whether it flips or not
                        // record the flip in the original matrix
                        // record the measurements

        /* Output */
            // a pritout of the lattice over time?
            // probabiliy of flipping at each stage
            // average spin
            // ... ?

    /* initialise everything */
    COMPLEX weightBefore, weightAfter, slice[lattice_size];
    float probability, prob, lambda, delta_tau;
    string result;
    int time_size, count = 0;

    /* calculate initial parameters */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, time_size, lattice_size);

    /* set up output headings */
    cout.width(11);
    cout << "weight";
    cout << " lattice" << endl;

    /* sweep through the lattice */
    for(int i = 0; i < iterations; i++){
        for(int t = 0; t < time_size; t++){

            for(int lattice_site = 0; lattice_site < lattice_size; lattice_site++){

                /* calculate the weight before the flip */
                general_weight(lattice_size, time_size, lattice, U, lambda, delta_tau, weightBefore);
                /* propose the flip */
                flip_scalar(slice[lattice_site]);

                /* calculate the weight after the flip */
                general_weight(lattice_size, time_size, lattice, U, lambda, delta_tau, weightAfter);

                /* calculate the ratio of weights */
                probability = weightAfter.r / weightBefore.r;

                /* accept or reject the flip */
                if(abs(probability) >= 1){
                    flip_scalar(lattice(t, lattice_site)); //accept
                    result = "accepted";
                }else{
                    prob = random_probability();
                    if(probability > prob){
                        flip_scalar(lattice(t, lattice_site)); //accept
                        result = "accepted";
                    }else{
                        // cout << " rejected" << endl;
                        result = "rejected";
                    }
                }
                /* comments */
                    //for negative values, we do some integration
                    //P\to\tilde{P} = |P| and  F\to \tilde
                    //you have to multiply each quan you measure bu the sign
                count++;
                cout << " (" << count <<") " << result << " - " << probability;
                cout.width(15);
                cout << " - weightBefore: " << weightBefore << ", weightAfter: " << weightAfter << endl;
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
}
void test_general_sweep(){
    /* Plan */
        /* [x] Input */
            // [x] lattice_size  - int
            // [x] iterations   - int
            // [x] lattice      - LaGenMatComplex
            // [x] U            - float
            // [x] beta         - float

        /* [x] Processing */
            // [x] generate a lattice of spins
            // [x] sweep the lattice

        /* [ ] Output */
            // [ ] average spins
            // [ ] acceptance probabilities

    /* initialise everything */
    int lattice_size = 5, iterations = 3;
    float U = 1, beta = 10;
    LaGenMatComplex lattice;

    /* generate the lattice */
    generate_lattice_matrix(lattice_size, lattice);

    /* sweep the lattice */
    general_sweep(lattice_size, lattice, U, beta, iterations);
}


/* --- Main QMC Program --- */
int main(){

    cout << "---- TESTING GENERALISED SWEEP ----" << endl;
    test_general_sweep();
    /* notes */

}
