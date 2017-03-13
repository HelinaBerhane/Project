/* ------ TO CONVERT ------ */

// [x] void ...(...)        -> void ..._f(...)
// [x] void ...(...)        -> void ...(..., const string file)

// [x] print_scalar         -> print_scalar_f
// [x] print_scalar_f(...)  -> print_scalar_f(...,"")
// [x] print_matrix         -> print_matrix_f
// [x] print_matrix_f(...)  -> print_matrix_f(...,"")

// [ ] cout << ...          -> myfile << ...
// [ ] storage              -> array

// [ ] void test_...()      -> void test_(file)
// [ ] open the files
// [ ] close the files


// - Example
void test_output_to_file(const string file){
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
    /* open the file */
    ofstream myfile;
    myfile.open(file);
    /* print stuff */
	myfile << name << ": ";
    for(int i = 0; i < array_size; i++){
        myfile << array[i] << " ";
    }
    myfile << endl;
    /* close the file */
    myfile.close();
}
void print_vector_f(const LaVectorComplex& vector, const string name, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file);
    /* print stuff */
    myfile << name << ":" << endl << vector << endl;
    /* close the file */
    myfile.close();
}

void print_initial_parameters_f(double U, double beta, double lambda, double delta_tau, double mu, int time_size, int lattice_size, const string file){
    /* open the file */
    ofstream myfile;
    myfile.open(file);
    /* print stuff */
	myfile << "no of lattice points = " << lattice_size << endl;
	myfile << "no of time slices = " << time_size << endl;
	myfile << "U = " << U << endl;
	myfile << "beta = " << beta << endl;
	myfile << "lambda = " << lambda << endl;
	myfile << "delta tau = " << delta_tau << endl;
    myfile << "mu = " << mu << endl;
    /* close the file */
    myfile.close();
}

/* -- Processing -- */
// Manipulation
void flip_spin_f(LaGenMatComplex& lattice, const int t, const int l){
    /* open the file */
    ofstream myfile;
    myfile.open(file);
    /* print stuff */
    myfile << "flipped ("<<t<<", "<<l<<"): " << lattice(t,l);
    lattice(t,l).r = -lattice(t,l).r;
    lattice(t,l).i = -lattice(t,l).i;
    myfile << " -> " << lattice(t,l) << endl;
    /* close the file */
    myfile.close();
}
// Calculation
// - qmc
void B_calculation_f(const COMPLEX slice[], const int lattice_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& B){
    /* open the file */
    ofstream myfile;
    myfile.open(file);
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
    /* close the file */
    myfile.close();
}
void O_calculation_f(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double sigma, const double delta_tau, const double mu, LaGenMatComplex& O){
    /* open the file */
    ofstream myfile;
    myfile.open(file);
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
        myfile << "t = " << t << ": ";
        isolate_row(lattice, lattice_size, t, slice);
        // print_array(slice, lattice_size, "slice");
        B_calculation_v(slice, lattice_size, U, lambda, sigma, delta_tau, mu, B);
        // print_matrix_f(B, "B");
        matrix_product(O, B);
        print_matrix_f(O, "product");
    }
    /* add I */
    matrix_sum(lattice_size, O, I);
    /* close the file */
    myfile.close();
}
void weight_calculation_f(const LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, const double mu, COMPLEX& weight){
    /* open the file */
    ofstream myfile;
    myfile.open(file);
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
    /* close the file */
    myfile.close();
}
void sweep_lattice_f(LaGenMatComplex& lattice, const int lattice_size, const int time_size, const double U, const double lambda, const double delta_tau, const double mu, const int iterations, double& acceptance, double& rejection){
    /* open the file */
    ofstream myfile;
    myfile.open(file);

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
                    myfile << " (" << count <<") " << "[" << acceptance << "/" << rejection << "] " << result << " - probability: " << probability;
                    myfile.width(15);
                    myfile << " - weightBefore: " << weightBefore << ", weightAfter: " << weightAfter << endl;
                }
                // if(result == "accepted"){
                //     print_matrix_f(lattice);
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
void test_weight(){
    string file = "test.txt";
    /* open the file */
    ofstream myfile;
    myfile.open(file);
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
    /* close the file */
    myfile.close();
}
void test_sweep(){
    string file = "test.txt";
    /* open the file */
    ofstream myfile;
    myfile.open(file);
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
    /* close the file */
    myfile.close();
}
void test_increasing_U(){
    string file = "test.txt";
    /* open the file */
    ofstream myfile;
    myfile.open(file);
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
    /* close the file */
    myfile.close();
}
