#include <iostream> //cout
#include <cstdlib>	//rand, sran
#include <string>
#include "complex_matrices.h"
#include <gmc.h> 	//LaGenMatComplex
#include <laslv.h>  //LUFactorizeIP, LaLUInverseIP, etc.
#include <blas3pp.h>
#include <random>   //random_device, mt19937

using namespace std;

/* Total [13/20] - QMC [3/3] */

/* Randomisation [1/1]*/
int ran(int max_rand){
    return rand() % max_rand;
}//working

/* Printing [7/7] */
void print_scalar(const COMPLEX scalar){
    cout << scalar << endl;
}//working
void print_scalar(const COMPLEX scalar, const string name){
    cout << name << ":" << scalar << endl;
}//working
void print_array(const COMPLEX array[], int len){
    for(int i = 0; i < len; i++){
        cout << array[i] << endl;
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
}
void print_matrix(const LaGenMatComplex& matrix){
	cout << matrix << endl;
}//working
void print_matrix(const LaGenMatComplex& matrix, const string name){
	cout << name << ":" << endl << matrix << endl;
}//working

/* Number generation [4/4] - QMC [3/3]*/
void generate_scalar(COMPLEX& A, const int x){
    A.r = ran(1, x);	//1 to x
    A.i = ran(1, x);
}//working
void generate_scalar(int number, const int x){
    number = ran(1, x);	//1 to x
}//working
void generate_array(COMPLEX array[], const int len, const int x){
    for(int i = 0; i < len; i++){
        array[i].r = ran(1, x);	//1 to x
        array[i].i = ran(1, x);
	}
}//working
void generate_matrix(const int matrix_size, const int max_rand, LaGenMatComplex& matrix){
    int matrix_volume = matrix_size*matrix_size;
    COMPLEX elements[matrix_volume];
    generate_array(elements, matrix_volume, max_rand);
    matrix = LaGenMatComplex(elements, matrix_size, matrix_size, false);
}//working
// QMC
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
void generate_H(const int matrix_size, LaGenMatComplex& hamiltonian){
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
                elements[n].r = 1;
            }else{
                elements[n].r = 0;
            }
            elements[n].i = 0;
        }
        //cout << endl;
    }
    //cout << endl;
    hamiltonian = LaGenMatComplex(elements, matrix_size, matrix_size, false );
    /* print result */
}//working
void generate_lattice_array(const int time_slices, COMPLEX& elements[]){
    COMPLEX elements[time_slices];
    for(int i = 0; i < time_slices; i++){
        elements[i].r = generate_spins();
        elements[i].i = 0;
    }
}

/* Matrix conversion [3/3] */
void vec_to_array(const LaVectorComplex& vector, const int len, COMPLEX array[ ]){
    for(int i = 0; i < len; i++){
        array[i] = vector(i);
    }
}//working
void array_to_diag(COMPLEX array[], const int len, LaGenMatComplex& diag){
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

/* Scalar manipulation [9/9] */
int factorial(int x){
	if(x <= 1){
        return 1;
	}else{
        return x * factorial(x - 1);
	}
}//working
void scalar_addition(const COMPLEX& A, const COMPLEX& B , COMPLEX& result){
    result.r = A.r + B.r;
    result.i = A.i + B.i;
}//working
void scalar_addition(COMPLEX& result, const COMPLEX addition){//probably working
    result.r += addition.r;
    result.i += addition.i;
}
void scalar_multiplication(const COMPLEX& A, const int B, COMPLEX& result){//to test
    result.r = A.r * B;
    result.i = A.i * B;
}
void scalar_multiplication(const COMPLEX& A, const COMPLEX& B, COMPLEX& result){
    la::complex<double> laA = la::complex<double>(A); //convert to la::complex<double>
    la::complex<double> laB = la::complex<double>(B);
    la::complex<double> laResult = la::complex<double>(result);
    laResult = laA * laB;
    result = laResult.toCOMPLEX();
}//working
void scalar_product(COMPLEX& total, const COMPLEX& number){
    COMPLEX part;
    part.r = (total.r * number.r) - (total.i * number.i);
    part.i = (total.r * number.i) + (total.i * number.r);
    total = part;
}//working
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
        scalar_addition(result, total_division);
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
//        scalar_multiplication(division, rec_scalar_exp_step(step-1),  multiplication);
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
}

/* Matrix manipulation [7/7] - [0/3]*/
void matrix_sum(const int matrix_size, const LaVectorComplex& sum, const LaVectorComplex& matrix){//to test
    for(int i = 0; i < matrix_size){
        for(int j = 0; j < matrix_size; j++){
            scalar_addition(sum, matrix);
        }
    }
}
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
// QMC
void V_matrix_calculation(const int matrix_size, LaGenMatComplex& V){//in progress
    //V = ??
    V = LaGenMatComplex::eye(matrix_size, matrix_size);
}
void B_matrix_calculation(LaGenMatComplex& H, LaGenMatComplex& V, LaGenMatComplex& B, const int matrix_size, const int iterations){//in progress
    //B = exp(-H)exp(-V)
    /* initialise everything */
    LaGenMatComplex expH;
    LaGenMatComplex expV;
    /* calculate exponentials */
    matrix_exponential(H, matrix_size, iterations, Hexp);
    matrix_exponential(V, matrix_size, iterations, Vexp);
    /* print exponential matrices */
    print_matrix(Hexp, "e^H");
    print_matrix(Vexp, "e^V");
    /* multiply exponentials */
    B = Hexp.copy;
    matrix_product(B, Vexp);
    /* print result */
    print_matrix(B, "B");
}
void O_matrix_calculation(const LaGenMatComplex& BA, const LaGenMatComplex& BB, const LaGenMatComplex& BC, const LaGenMatComplex& BD, const LaGenMatComplex&BE, const LaGenMatComplex& O){//in progress
    //O = 1 + B(m) B(m-1) B(...) B(1)
    /* initialise everything */
    LaGenMatComplex I = LaGenMatComplex::eye(matrix_size, matrix_size);
    //LaGenMatComplex multiplication;
    /* multiply exponentials */
    five_matrix_multiplication(BA, BB, BC, BD, BE, O);
    /* add I */
    matrix_sum(O, I);
}
void partition_function(){//in progress
    //
}

/* Testing [11/11] - QMC [2/5]*/
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
        scalar_addition(sum, compA);
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
        scalar_addition(number, step);
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
    for(int j = 1; j <= iterations; j++){
        matrix_exponential(matrix, 3, j, result);
        cout << j << " iterations:" << endl;
        print_matrix(result);
    }
}//working
// QMC
void test_lattice_generation(const int matrix_size, const int time_slices){
    LaGenMatComplex lattice;
    generate_lattice(matrix_size, lattice);
    print_matrix(lattice);
    COMPLEX string[time_slices];
    generate_lattice_array(time_slices, elements[]);
    print_array(elements, "string");
}//working
void test_hamiltonian(const int matrix_size){
    /* initialise everything */
    LaGenMatComplex hamiltonian;
    LaVectorComplex eigenvalues = LaVectorComplex(matrix_size);
    LaGenMatComplex eigenvectors = LaGenMatComplex::zeros(matrix_size, matrix_size);
    /* generate matrix */
    generate_H(matrix_size, hamiltonian);
    print_matrix(hamiltonian);
    /* calculate eigenstuff */
    matrix_eigenvstuff(hamiltonian, eigenvalues, eigenvectors);
    print_vector(eigenvalues, "eigenvalues");
    // eigenvalues are 2 cos(n pi / q), where q = the matrix size
}//working
void test_V_generation(){//in progress
    //delta test_scalar_sum
    //lamba
    //sigma
    //s_i(l)
    //mu
}
void test_B_generation(){//in progress
    /* initialise everything */
    int delta_tau = 1;
    LaGenMatComplex hamiltonian;
    /* generate matrices */
    generate_H(matrix_size, hamiltonian);
}
void test_QMC(const int matrix_size, const int time_slices){
    test_lattice_generation(matrix_size);
    test_hamiltonian(time_slices);
}

/* --- Main QMC Program --- */
int main(){
    /* initialise everything */
    int matrix_size = 5; time_slices = 5;//tau
    /* generate matrices */
    test_lattice_generation(const int matrix_size, const int time_slices);
    //test_QMC(matrix_size, time_slices);
}
