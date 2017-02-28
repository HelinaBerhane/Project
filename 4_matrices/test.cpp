void multiply( LaGenMatComplex* working, LaGenMatComplex** matrices ) {
    Blas_Mat_Mat_Mult(working, matrices[0], working);
    multiply( working, matrices + 1 );
}

int main() {
    LaGenMatComplex* matrices[n] = new LaGenMatComplex[n];
    LaGenMatComplex multiply = new LaGenMatComplex;
    //make this new matrix in to a identity matrix

    multiply( multiply, matrices );
}
