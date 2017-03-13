/* ------ TO CONVERT ------ */

// [x] void ...(...)        -> void ..._f(...)
// [x] void ...(...)        -> void ...(..., const string file)

// [x] print_scalar         -> print_scalar_f
// [x] print_scalar_f(...)  -> print_scalar_f(...,"")
// [x] print_matrix         -> print_matrix_f
// [x] print_matrix_f(...)  -> print_matrix_f(...,"")

// [ ] cout << ...          -> myfile << ...
// [ ] storage              -> array


// - Example
void test_output_to_file(){
    string file = test;
    /* open the file */
    ofstream myfile;
    myfile.open(file);
    /* print stuff */
    myfile << "Here, have a test :) .\n";
    /* close the file */
    myfile.close();
}

/* -- Output -- */

void print_array_f(const COMPLEX array[], int array_size, const string name, const string file){
	cout << name << ": ";
    for(int i = 0; i < array_size; i++){
        cout << array[i] << " ";
    }
    cout << endl;
}
void print_vector_f(const LaVectorComplex& vector, const string name, const string file){
    cout << name << ":" << endl << vector << endl;
}

void print_initial_parameters_f(double U, double beta, double lambda, double delta_tau, double mu, int time_size, int lattice_size, const string file){
	cout << "no of lattice points = " << lattice_size << endl;
	cout << "no of time slices = " << time_size << endl;
	cout << "U = " << U << endl;
	cout << "beta = " << beta << endl;
	cout << "lambda = " << lambda << endl;
	cout << "delta tau = " << delta_tau << endl;
    cout << "mu = " << mu << endl;
}

/* -- Processing -- */
// Manipulation
void flip_spin_f(LaGenMatComplex& lattice, const int t, const int l){
    cout << "flipped ("<<t<<", "<<l<<"): " << lattice(t,l);
    lattice(t,l).r = -lattice(t,l).r;
    lattice(t,l).i = -lattice(t,l).i;
    cout << " -> " << lattice(t,l) << endl;
}
// Calculation
// - qmc
void B_calculation_f(const COMPLEX slice[], const int lattice_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& B){
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
    print_matrix_f(sum, "H");
    print_matrix_f(V, "V");

    /* calculate H + V */
    matrix_sum(lattice_size, sum, V);
    print_matrix_f(sum, "H + V");

    /* calculate delta_tau * (H + V) */
    matrix_multiple(sum, lattice_size, delta_tau, product);
    print_matrix_f(product, "delta_tau * (H + V)");

    /* calculate - delta_tau * (H + V) */
    matrix_negative(lattice_size, product, negative);
    print_matrix_f(negative, "- delta_tau * (H + V)");

    /* calculate exp(- delta_tau * (H + V)) */
    matrix_exponential(negative, lattice_size, B);
    print_matrix_f(B, "B = exp(- delta_tau * (H + V))");
}
void O_calculation_f(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& O){
    /* initialise everything */
    LaGenMatComplex B = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex I = LaGenMatComplex::eye(lattice_size, lattice_size);
    COMPLEX slice[lattice_size];
    O = LaGenMatComplex::eye(lattice_size, lattice_size);
    print_matrix_f(O, "product");
    /* calculate B matrices */
    for(int x = 0; x < time_size; x++){
        clear_storage(slice, lattice_size);
        int t = time_size - x - 1;
        cout << "t = " << t << ": ";
        isolate_row(lattice, lattice_size, t, slice);
        // print_array(slice, lattice_size, "slice");
        B_calculation_v(slice, lattice_size, U, lambda, sigma, delta_tau, mu, B);
        // print_matrix_f(B, "B");
        matrix_product(O, B);
        print_matrix_f(O, "product");
    }
    /* add I */
    matrix_sum(lattice_size, O, I);
}
void weight_calculation_f(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, const double mu, COMPLEX& weight){
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
    print_matrix_f(OUP, "O UP");
    print_matrix_f(ODN, "O DN");
    /* calculate det(O) */
    matrix_determinant_e(lattice_size, OUP, detOUP);
    matrix_determinant_e(lattice_size, ODN, detODN);
    print_scalar_f(detOUP, "det(O UP)");
    print_scalar_f(detODN, "det(O DN)");
    /* calculate weight */
    weight = scalar_multiple(detOUP, detODN);
    print_scalar_f(weight, "weight");
}
void sweep_lattice_f(LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, const double mu, const int iterations, double& acceptance, double& rejection){

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
                //     print_matrix_f(lattice);
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

/* -- Testing -- */
// - generic
void test_output_to_file(){
    ofstream myfile;
    myfile.open ("test.txt");
    myfile << "Here, have a test :) .\n";
    myfile.close();
}
void test_flip_spins(){
    /* initialise stuff */
    int lattice_size = 5, time_size = 8;
    int l = random_int(lattice_size-1), t = random_int(time_size-1);
    LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
    /* generate lattice */
    generate_lattice(lattice_size, time_size, lattice);
    print_matrix_f(lattice, "lattice");
    /* flip spins */
    flip_spin_v(lattice, t, l);
}
void test_inverse(){
    int matrix_size = 3;
    LaGenMatComplex initialMatrix = LaGenMatComplex::rand(matrix_size, matrix_size, 0, 9);
    LaGenMatComplex inverseMatrix = LaGenMatComplex::zeros(matrix_size, matrix_size);
    matrix_inverse(initialMatrix, matrix_size, inverseMatrix);
    print_matrix_f(initialMatrix, "inverse matrix");
    print_matrix_f(inverseMatrix, "inverse matrix");
}
void test_matrix_product(){
    /* initialise everything */
    int matrix_size = 5;
    LaGenMatComplex matrixA = LaGenMatComplex::rand(matrix_size, matrix_size, 0, 9);
    LaGenMatComplex matrixB = LaGenMatComplex::rand(matrix_size, matrix_size, 0, 9);
    /* print everything */
    print_matrix_f(matrixA, "Matrix A");
    print_matrix_f(matrixB, "Matrix B");
    /* matrix product */
    matrix_product(matrixA, matrixB);
    print_matrix_f(matrixB, "result");
}
void test_recombine_diagonalised_matrices(){
    /* initialise everything */
    int matrix_size = 5;
    LaGenMatComplex matrix = LaGenMatComplex::rand(matrix_size, matrix_size, 0, 9);
    LaVectorComplex eigenvalues = LaVectorComplex(matrix_size);
    LaGenMatComplex eigenvectors = LaGenMatComplex::zeros(matrix_size, matrix_size);
    LaGenMatComplex result = LaGenMatComplex::zeros(matrix_size, matrix_size);
    /* generate matrix */
    print_matrix_f(matrix, "initial matrix");
    /* calculate eigenstuff */
    LaEigSolve(matrix, eigenvalues, eigenvectors);
    /* multiply them back together */
    recombine_diagonalised_matrices(matrix_size, eigenvectors, eigenvalues, result);
    print_matrix_f(result, "final matrix");
}
void test_matrix_equals_(){
    int matrix_size = 5;
    LaGenMatComplex matrix = LaGenMatComplex::rand(matrix_size, matrix_size, 0, 5);
    print_matrix_f(matrix, "initial matrix");
    matrix = 0;
    print_matrix_f(matrix, "matrix = 0");
    matrix = 1;
    print_matrix_f(matrix, "matrix = 1");
}
void test_matrix_negative(){
    int matrix_size = 3;
    LaGenMatComplex matrix = LaGenMatComplex::rand(matrix_size, matrix_size, 0, 9);
    LaGenMatComplex result = LaGenMatComplex::zeros(matrix_size, matrix_size);
    print_matrix_f(matrix, "Matrix");
    matrix_negative(matrix_size, matrix, result);
    print_matrix_f(result, "- Matrix");
    matrix_negative(matrix_size, matrix);
    print_matrix_f(matrix, "- Matrix (in place)");
}
void test_matrix_multiple(){
    int matrix_size = 4;
    LaGenMatComplex matrix = LaGenMatComplex::rand(matrix_size, matrix_size, 0, 9);
    LaGenMatComplex result = LaGenMatComplex::zeros(matrix_size, matrix_size);
    matrix_multiple(matrix, matrix_size, 2.0, result);
    print_matrix_f(matrix, "initial matrix");
    print_matrix_f(result, "initial matrix * 2");

}
void test_matrix_exponential(){
    int matrix_size = 5;
    /* initialise everything */
    LaGenMatComplex matrix = LaGenMatComplex::rand(matrix_size, matrix_size, 0, 9);
    LaGenMatComplex result = LaGenMatComplex::zeros(matrix_size, matrix_size);
    print_matrix_f(matrix, "initial matrix");
    /* calculate exponential */
    matrix_exponential(matrix, matrix_size, result);
    print_matrix_f(result, "e^(matrix)");
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
    print_matrix_f(test, "test");
    print_matrix_f(result, "result");
}
void test_store_matrix(){
    /* initialise everything */
    int matrix_size = 5, storage_size = matrix_size * matrix_size * 3;
    COMPLEX storage[storage_size];
    clear_storage(storage, storage_size);
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
void test_matrix_determinant(){
    /* initialise everything */
    int matrix_size = 4;
    LaGenMatComplex matrix = LaGenMatComplex::rand(matrix_size, matrix_size, 0, 9);
    print_matrix_f(matrix, "initial matrix");
    COMPLEX result;
    // clear_scalar(result);
    /* calculate determinant */
    // result = matrix_determinant(matrix_size, matrix);
    // print_scalar_f(result, "determinant");
    clear_scalar(result);
    matrix_determinant_e(matrix_size, matrix, result);
    print_scalar_f(result, "determinant (from eigenstuff)");
}
// - qmc
void test_initial_parameters(){
    double U = 1, beta = 10, lambda, delta_tau, mu;
    int lattice_size = 5, time_size;
    initial_parameter_calculation(U, beta, lambda, delta_tau, mu, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size);
}
void test_generate_lattice(){
    int lattice_size = 5, time_size = 17;
    LaGenMatComplex lattice;
    generate_lattice(lattice_size, time_size, lattice);
    print_matrix_f(lattice, "lattice");
}
void test_H(){
    /* initialise everything */
    int lattice_size = 5;
    LaGenMatComplex H;
    LaVectorComplex eigenvalues = LaVectorComplex(lattice_size);
    LaGenMatComplex eigenvectors = LaGenMatComplex::zeros(lattice_size, lattice_size);
    /* generate matrix */
    H_generation(lattice_size, H);
    print_matrix_f(H, "H");
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
    print_matrix_f(H, "H");
    /* calculate -H */
    matrix_negative(lattice_size, H, negH);
    print_matrix_f(negH, "-H");
    /* calculate exponentials */
    matrix_exponential_v(H, lattice_size, expH);
    print_matrix_f(expH, "e^(-H)");
}
void test_V(){
    /* initialise everything */
    int lattice_size = 5, time_size;
    LaGenMatComplex V = LaGenMatComplex::zeros(lattice_size, lattice_size);
    COMPLEX slice[lattice_size];
    double U = 1, beta = 10, lambda, delta_tau, mu;

    /* calculate initial parameters */
    initial_parameter_calculation(U, beta, lambda, delta_tau, mu, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size);
    cout << endl;

    /* generate the lattice */
    generate_slice(lattice_size, slice);
    print_array(slice, lattice_size, "slice");

    /* calculate V */
    V_calculation(slice, lattice_size, U, lambda, 1, delta_tau, mu, V);
    print_matrix_f(V, "V");
}
void test_B_calculation(){
    /* initialise everything */
    int lattice_size = 5, time_size;
    double U = 1, beta = 10, lambda, delta_tau, mu;
    LaGenMatComplex B = LaGenMatComplex::zeros(lattice_size, lattice_size);
    COMPLEX slice[lattice_size];
    /* generate initial conditions */
    initial_parameter_calculation(U, beta, lambda, delta_tau, mu, time_size);
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
    print_matrix_f(B,"B = e^-H e^-V");
}
void test_O(){
    /* initialise everything */
    int lattice_size = 5, time_size = 0;
    double U = 1, beta = 10, lambda, delta_tau, mu;
    LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
    LaGenMatComplex O = LaGenMatComplex::zeros(lattice_size, lattice_size);
    /* generate initial conditions */
    initial_parameter_calculation(U, beta, lambda, delta_tau, mu, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size);
    /* generate lattice */
    generate_lattice(lattice_size, time_size, lattice);
    print_matrix_f(lattice, "lattice");
    /* calculate O */
    O_calculation_v(lattice, lattice_size, time_size, U, lambda, 1, delta_tau, mu, O);
    print_matrix_f(O, "O");
}
void test_weight(){
    /* initialise stuff */
    int lattice_size = 5, time_size;
    double U = 1, beta = 10, lambda, delta_tau, mu;
    COMPLEX weight;
    weight.r = 0;
    weight.i = 0;
    /* generate initial conditions */
    initial_parameter_calculation(U, beta, lambda, delta_tau, mu, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size);
    /* generate lattice */
    LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
    generate_lattice(lattice_size, time_size, lattice);
    print_matrix_f(lattice, "lattice");
    /* calculate the weight */
    weight_calculation(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, weight);
    print_scalar_f(weight, "weight");
}
void test_sweep(){
    /* initialise everything */
    int lattice_size = 5, time_size, iterations = 1000;// = 10000;
    double U = .1, beta = 1, lambda, delta_tau, mu;
    double acceptance = 0, rejection = 0;
    /* generate initial conditions */
    initial_parameter_calculation(U, beta, lambda, delta_tau, mu, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size);
    /* generate lattice */
    LaGenMatComplex lattice = LaGenMatComplex::zeros(time_size, lattice_size);
    // print_matrix_f(lattice, "intialised lattice");
    generate_lattice(lattice_size, time_size, lattice);
    print_matrix_f(lattice, "lattice");
    /* sweep the lattice */
    sweep_lattice_v(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, iterations, acceptance, rejection);
}
void test_increasing_U(){
    /* initialise everything */
    int lattice_size = 5, time_size = 0, iterations = 120;
    double U, beta = 5.0, lambda = 1.0, delta_tau, mu;
    double acceptance = 0.0, rejection = 0.0;
    /* test U = 0 to 10 */
    for(int i = 1; i <= 10; i++){
        /* generate initial conditions */
        U = i;
        initial_parameter_calculation(U, beta, lambda, delta_tau, mu, time_size);
        print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size);
        /* generate a lattice of spins */
        LaGenMatComplex lattice = LaGenMatComplex::zeros(time_size, lattice_size);
        generate_lattice(lattice_size, time_size, lattice);
        /* sweep the lattice */
        sweep_lattice_v(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, iterations, acceptance, rejection);
    }
}
