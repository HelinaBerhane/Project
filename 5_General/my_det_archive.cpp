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
}
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
}
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
}
COMPLEX matrix_determinant(const int matrix_size, const LaGenMatComplex& matrix){
    /* initialise everything */
    COMPLEX determinant;
    LaGenMatComplex cofactorMatrix;
    COMPLEX coefficient;
    cofactorMatrix = 0;
    /* do stuff */
    if(matrix_size == 2){
        return simple_matrix_determinant(matrix);
    }else{
        //for each element in the first row
        for(int element = 0; element < matrix_size; element++){
            /* initialise everything */
            int cofactor_size = matrix_size - 1;
            clear_scalar(determinant);
            clear_scalar(coefficient);
            cofactorMatrix = LaGenMatComplex::zeros(cofactor_size, cofactor_size);
            /* determine the coefficient */
            coefficient = determinant_coefficient(matrix, element);
                // = +- the element
            /* calculate the cofactor */
            generate_cofactor_matrix(matrix_size, matrix, element, cofactorMatrix);
            //print_matrix(cofactorMatrix, "cofactorMatrix");
            /* finish calculation */
            scalar_sum(determinant, scalar_multiple(coefficient, matrix_determinant(cofactor_size, cofactorMatrix)));
        }
    }
    return determinant;
}
