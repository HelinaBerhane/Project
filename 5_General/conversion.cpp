/* Randomisation */
int basic_random_int(int max_rand){
    return rand() % (max_rand+1);
}
float basic_random_float(){
    return basic_random_int(1000)/1000;
}
float random_float(float min, float max){
    random_device rd;
    mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}
float random_probability(){
    random_device rd;
    mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    return dis(gen);
}

/* Printing [7/7] */
void print_array(const COMPLEX array[], int len){
    for(int i = 0; i < len; i++){
        cout.width(7);
        cout << array[i] << " ";
    }
    cout << endl;
}//working
void print_vector(const LaVectorComplex& vector, const string name){
    cout << name << ":" << endl << vector << endl;
}//working


// Manipulation
void copy_array(const int len, const COMPLEX array[], COMPLEX copy[]){//in progress
    for(int i = 0; i < len; i++){
        //
    }
}

/* Scalar manipulation [14/14] */
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
void scalar_multiplication(const COMPLEX& A, const int B, COMPLEX& result){//to test
    result.r = A.r * B;
    result.i = A.i * B;
}//working
void scalar_multiplication_f(const COMPLEX& A, const float B, COMPLEX& result){//to test
    result.r = A.r * B;
    result.i = A.i * B;
}//working

COMPLEX scalar_multiple(COMPLEX& A, const COMPLEX& B){
    COMPLEX part;
    part.r = (A.r * B.r) - (A.i * B.i);
    part.i = (A.r * B.i) + (A.i * B.r);
    return part;
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
}

void test_random_probability(){
    int count = 10;
    for (int i = 0; i < count; i++) {
        cout << random_probability() << endl;
    }
}//working
void test_parameter_calculation(){
    float U = 1, lambda = lambda_calculation(U), delta_tau = delta_tau_calculation(U);
    cout << "U = " << U << endl;
    cout << "l = " << lambda << endl;
    cout << "t = " << delta_tau << endl;
}//working

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
    generate_lattice_matrix(matrix_size, matrix_size, lattice);

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
        generate_lattice_matrix(matrix_size, matrix_size, lattice);

        /* sweep the lattice */
        sweep_lattice(matrix_size, lattice, U, iterations);
    }
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
    int iterations = 4;
    int lattice_volume = lattice_size * time_size;
    int matrix_volume = lattice_size * lattice_size;
    int storage_size = matrix_volume * time_size;
    float sigma = 0.0, beta = 10.0;
    LaGenMatComplex H;
    LaGenMatComplex V;
    LaGenMatComplex B;
    LaGenMatComplex O;
    LaGenMatComplex I = LaGenMatComplex::eye(lattice_size, lattice_size);
    weight.r = 1;
    weight.i = 0;
    COMPLEX detO;
    COMPLEX slice[lattice_size];
    COMPLEX storage[storage_size];
    for(int i = 0; i < storage_size; i++){
        storage[i].r = 0;
        storage[i].i = 0;
    }
    // print_array(storage, storage_size, "storage");

    /* generate H */
    generate_H(lattice_size, H);

    // for sigma = +-1
    for(int s = 0; s < 2; s++){
        // reset all variables
        O = LaGenMatComplex::eye(lattice_size, lattice_size);
        sigma = s * 2.0 - 1.0;
        // cout << "sigma = " << sigma << endl;
        // for each time_slice
        // cout << "no of time slices = " << time_size << endl;

        for(int t = 0; t < time_size; t++){    // check the order of multiplication!!!
            // cout << "current time slice = " << t << endl;
            // reset all variables
            V = LaGenMatComplex::zeros(lattice_size, lattice_size);
            B = LaGenMatComplex::zeros(lattice_size, lattice_size);

            // isolate the time slice
            isolate_row(lattice, lattice_size, t, slice);
            // print_matrix(lattice, "lattice");
            // cout << endl << "current time slice = " << t << endl;
            // print_array(slice, lattice_size, "lattice");

            // calculate all variables
            // print_array(slice, lattice_size, "slice");
            V_calculation(slice, lattice_size, U, lambda, sigma, delta_tau, V);
            // print_matrix(V, "V");
            B_calculation(slice, lattice_size, U, lambda, sigma, delta_tau, B);

            // print_matrix(B, "B");

            // store the elements of the B matrix in an array
            for(int r = 0; r < lattice_size; r++){
                for(int c = 0; c < lattice_size; c++){
                    int i = (t * matrix_volume) + (r * lattice_size) + c;
                    // cout << "i = " << i << endl;
                    storage[i].r = B(r, c).r;
                    storage[i].i = B(r, c).i;
                    // print_scalar(storage[i]);
                }
            }
            print_array(storage, storage_size, "storage");
            // cout << endl;
            // print_array(storage, storage_size, "storage");
            /* multiply the matrices */
            print_matrix(O, "B_1 B_2 B_...");
        }
        n_matrix_product(storage, lattice_size, time_size, O);
        // print_array(storage, storage_size, "storage");
        // add 1
        matrix_sum(lattice_size, O, I);
        print_matrix(O, "O");
        // calculate det O
        matrix_determinant(lattice_size, O, detO);
        print_scalar(detO, "detO");
        // calculate det O up * det O down
        scalar_product(weight, detO);
        print_scalar(weight, "weight");
    }
}
void general_sweep(const int lattice_size, LaGenMatComplex& lattice, const float U, const float beta, const float lambda, const float delta_tau, const int time_size, const int iterations){
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
    float probability, prob;
    string result;
    int count = 0;

    /* set up output headings */
    cout.width(11);
    cout << "weight";
    cout << " lattice" << endl;

    /* sweep through the lattice */
    for(int i = 0; i < iterations; i++){
        for(int t = 0; t < time_size; t++){

            for(int l = 0; l < lattice_size; l++){

                /* calculate the weight before the flip */
                general_weight(lattice_size, time_size, lattice, U, lambda, delta_tau, weightBefore);
                /* propose the flip */
                // print_scalar(lattice(t,l), "before");

                flip_scalar(lattice(t,l));
                // print_scalar(lattice(t,l), "after");

                /* calculate the weight after the flip */
                general_weight(lattice_size, time_size, lattice, U, lambda, delta_tau, weightAfter);

                /* calculate the ratio of weights */
                probability = weightAfter.r / weightBefore.r;

                /* accept or reject the flip */
                if(abs(probability) >= 1){
                    result = "accepted";
                }else{
                    prob = random_probability();
                    if(probability > prob){
                        result = "accepted";
                    }else{
                        // cout << " rejected" << endl;
                        flip_scalar(lattice(t, l));
                        result = "rejected";
                    }
                }

                // print_scalar(lattice(t,l), "after");
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
    int lattice_size = 5, iterations = 3, time_size;
    float U = 5, beta = 10, lambda, delta_tau;
    LaGenMatComplex lattice;

    /* calculate initial parameters */
    initial_parameter_calculation(U, beta, lambda, delta_tau, time_size);
    print_initial_parameters(U, beta, lambda, delta_tau, time_size, lattice_size);
    int lattice_volume = lattice_size * time_size;

    /* generate the lattice */
    generate_lattice_matrix(lattice_size, time_size, lattice);

    /* sweep the lattice */
    general_sweep(lattice_size, lattice, U, beta, lambda, delta_tau, time_size, iterations);
}
