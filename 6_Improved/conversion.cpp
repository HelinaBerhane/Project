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



/* -- Testing -- */
// - generic



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
