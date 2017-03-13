/* ------ TO CONVERT ------ */

// [x] void ...(...)        -> void ...(...)
// [x] void ...(...)        -> void ...(..., const string file)

// [x] print_scalar         -> print_scalar_f
// [x] print_scalar(...)  -> print_scalar(..., " ", file)
// [x] print_array          -> print_array_f
// [x] print_array(...)   -> print_array(..., " ", file)
// [x] print_matrix         -> print_matrix_f
// [x] print_matrix(...)  -> print_matrix(..., " ", file)

// [x] print_initial_parameters        -> print_initial_parameters_f
// [x] print_initial_parameters(...) -> print_initial_parameters(...,file)

// [x] cout << ... -> myfile << ...
// [x] storage     -> array

// [x] void test_...()      -> void test_(file)
// [x] open the files
// [x] close the files

// [x] flip_spin(...) -> flip_spin(..., file)


// - Example
void test_output_to_file(const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* print stuff */
    myfile << "Here, have a test :) .\n";
    /* close the file */
    myfile.close();
}

/* -- Processing -- */
// Calculation
// - qmc

void O_calculation(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& O, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* initialise everything */
    LaGenMatComplex B = LaGenMatComplex::zeros(lattice_size, lattice_size);
    LaGenMatComplex I = LaGenMatComplex::eye(lattice_size, lattice_size);
    COMPLEX slice[lattice_size];
    O = LaGenMatComplex::eye(lattice_size, lattice_size);
    print_matrix(O, "product", file);
    /* calculate B matrices */
    for(int x = 0; x < time_size; x++){
        clear_array(slice, lattice_size);
        int t = time_size - x - 1;
        myfile << "t = " << t << ": ";
        isolate_row(lattice, lattice_size, t, slice, file);
        // print_array(slice, lattice_size, "slice", file);
        B_calculation_v(slice, lattice_size, U, lambda, sigma, delta_tau, mu, B);
        // print_matrix(B, "B", file);
        matrix_product(O, B);
        print_matrix(O, "product", file);
    }
    /* add I */
    matrix_sum(lattice_size, O, I);
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
    myfile << "sigma = 1" << endl;
    O_calculation_v(lattice, lattice_size, time_size, U, lambda, 1, delta_tau, mu, OUP);
    myfile << "sigma = -1" << endl;
    O_calculation_v(lattice, lattice_size, time_size, U, lambda, -1, delta_tau, mu, ODN);
    print_matrix(OUP, "O UP", file);
    print_matrix(ODN, "O DN", file);
    /* calculate det(O) */
    matrix_determinant_e(lattice_size, OUP, detOUP);
    matrix_determinant_e(lattice_size, ODN, detODN);
    print_scalar(detOUP, "det(O UP)", file);
    print_scalar(detODN, "det(O DN)", file);
    /* calculate weight */
    weight = scalar_multiple(detOUP, detODN);
    print_scalar(weight, "weight", file);
    /* close the file */
    myfile.close();
}
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

/* -- Testing -- */
// - generic
void test_O(const string file){
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
    print_matrix(lattice, "lattice");
    /* calculate O */
    O_calculation_v(lattice, lattice_size, time_size, U, lambda, 1, delta_tau, mu, O);
    print_matrix(O, "O");
}
void test_weight(const string file){
    string file = "test.txt";
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* initialise stuff */
    int lattice_size = 5, time_size;
    double U = 1, beta = 10, lambda, delta_tau, mu;
    COMPLEX weight;
    weight.r = 0;
    weight.i = 0;
    /* generate initial conditions */
    initial_parameter_calculation(U, beta, lambda, delta_tau, mu, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, mu, time_size, lattice_size, file);
    /* generate lattice */
    LaGenMatComplex lattice = LaGenMatComplex::zeros(lattice_size, time_size);
    generate_lattice(lattice_size, time_size, lattice);
    print_matrix(lattice, "lattice", file);
    /* calculate the weight */
    weight_calculation(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, weight);
    print_scalar(weight, "weight", file);
    /* close the file */
    myfile.close();
}
void test_sweep(const string file){
    string file = "test.txt";
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
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
    /* close the file */
    myfile.close();
}
void test_increasing_U(const string file){
    string file = "test.txt";
    /* open the file */
    ofstream myfile;
    myfile.open(file, std::ios_base::app);
    /* initialise everything */
    int lattice_size = 5, time_size = 0, iterations = 120;
    double U, beta = 5.0, lambda = 1.0, delta_tau, mu;
    double acceptance = 0.0, rejection = 0.0;
    /* test U = 0 to 10 */
    for(int i = 1; i <= 10; i++){
        /* generate initial conditions */
        U = i;
        initial_parameter_calculation(U, beta, lambda, delta_tau, mu, time_size);
        print_initial_parameters_f(U, beta, lambda, delta_tau, mu, time_size, lattice_size, file);
        /* generate a lattice of spins */
        LaGenMatComplex lattice = LaGenMatComplex::zeros(time_size, lattice_size);
        generate_lattice(lattice_size, time_size, lattice);
        /* sweep the lattice */
        sweep_lattice_v(lattice, lattice_size, time_size, U, lambda, delta_tau, mu, iterations, acceptance, rejection);
    }
    /* close the file */
    myfile.close();
}
