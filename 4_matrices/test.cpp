void multiply( LaGenMatComplex* product, LaGenMatComplex** matrices, int n) {
    if( n == 0 ) {
        return;
    }
    Blas_Mat_Mat_Mult(product, matrices[0], product);
    multiply(product, matrices + 1, n - 1 );
}

int main() {
    LaGenMatComplex* matrices[n] = new LaGenMatComplex[n];
    LaGenMatComplex multiply = new LaGenMatComplex;
    //make this new matrix in to a identity matrix

    multiply( multiply, matrices, 5 );
}




void multiply_matrices(LaGenMatComplex* product, LaGenMatComplex** matrices, int n){
    if(n <= 0){
        return;
    }
    Blas_Mat_Mat_Mult(product, matrices[0], product);
    multiply(product, matrices + 1, n - 1 );
}
