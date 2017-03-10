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

void print_vector(const LaVectorComplex& vector, const string name){
    cout << name << ":" << endl << vector << endl;
}//working

/* Generation [5/5]*/

void generate_array(COMPLEX array[], const int array_length, const int max_rand){
    for(int i = 0; i < array_length; i++){
        array[i].r = basic_random_int(max_rand);	//1 to x
        array[i].i = basic_random_int(max_rand);
	}
}//working


void generate_general_matrix(const int matrix_width, const int matrix_length, const int max_rand, LaGenMatComplex& matrix){
    int matrix_volume = matrix_width * matrix_length;
    COMPLEX elements[matrix_volume];
    generate_array(elements, matrix_volume, max_rand);
    matrix = LaGenMatComplex(elements, matrix_length, matrix_width, true);
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

void copy_array(const int len, const COMPLEX array[], COMPLEX copy[]){//in progress
    for(int i = 0; i < len; i++){
        //
    }
}
void isolate_row(const LaGenMatComplex& matrix, const int matrix_width, const int row, COMPLEX array[]){
    for(int i = 0; i < matrix_width; i++){
        array[i] = matrix(row, i);
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

void scalar_multiplication(const COMPLEX& A, const COMPLEX& B, COMPLEX& result){
    la::complex<double> laA = la::complex<double>(A); //convert to la::complex<double>
    la::complex<double> laB = la::complex<double>(B);
    la::complex<double> laResult = la::complex<double>(result);
    laResult = laA * laB;
    result = laResult.toCOMPLEX();
}//working

COMPLEX scalar_multiple(COMPLEX& A, const COMPLEX& B){
    COMPLEX part;
    part.r = (A.r * B.r) - (A.i * B.i);
    part.i = (A.r * B.i) + (A.i * B.r);
    return part;
}

void matrix_sum(const int matrix_size, LaGenMatComplex& sum, const LaGenMatComplex& matrix){//to test
    for(int i = 0; i < matrix_size; i++){
        for(int j = 0; j < matrix_size; j++){
            sum(i, j).r += matrix(i, j).r;
            sum(i, j).i += matrix(i, j).i;
        }
    }
}//should be working

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
    V_calculation(latticeUP, lattice_size, U, lambda, 1, delta_tau, VUP);
    V_calculation(latticeDOWN, lattice_size, U, lambda, 1, delta_tau, VDOWN);

    /* multiply B matrices */
    for(int t = time_size - 1; t >= 0 ; t--){
        /*   for each time slice   */

        /* calculate B(t) matrices */
        B_calculation(latticeUP, lattice_size, U, lambda, 1, delta_tau, BUP);
        B_calculation(latticeDOWN, lattice_size, U, lambda, 1, delta_tau, BDOWN);

        /* multiply the matrices */
        matrix_product(proBUP, BUP);
        matrix_product(proBDOWN, BDOWN);
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
    int matrix_width = 3, matrix_length = 5, max_rand = 9;

    LaGenMatComplex matrix;
    COMPLEX array[matrix_width];

    /* generate matrix */
    generate_general_matrix(matrix_width, matrix_length, max_rand, matrix);
    print_matrix(matrix, "matrix");

    /* isolate and print rows */
    for(int row = 0; row < matrix_length; row++){
        isolate_row(matrix, matrix_width, row, array);
        cout << " Row (" << row << "):";
        print_array(array, matrix_width);
    }
}//working

// QMC [8/10]
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

/* Generalisation */
/* - Working - */


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


/* --- Main QMC Program --- */
int main(){

        cout << "---- TESTING MATRIX MANIPULATION ----" << endl;
    test_matrix_equals_();
    /* notes */

}
