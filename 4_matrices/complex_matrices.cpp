#define LA_COMPLEX_SUPPORT

#include <iostream> //cout
#include <cstdlib>	//rand, sran
#include <gmc.h> 	//LaGenMatComplex
#include <string>
#include <laslv.h>  //LUFactorizeIP, LaLUInverseIP, etc.
#include <lavli.h>  //LaVectorLongInt
using namespace std;

/* Testing */
void print_matrix(LaGenMatComplex* matrix, const string name){					//working
	cout << name << ":" << endl << *matrix << endl;
}

/* Matrix generation */
void generate_elements(COMPLEX* elements, int m, int x){	//working
	srand(time(NULL));				//seed
    for(int i = 0; i < m*m; i++){
        elements[i].r = rand() % x;	//1 to x
        elements[i].i = rand() % x;
	}
}

/* Matrix manipulation */
void matrix_eigenvalues(){
    
}                                                   //empty                       
void matrix_eignevectors(LaGenMatComplex* matrix){
	/* see laslv.h */
}                           //empty
void matrix_powers(){
	/* documentation {
		void Blas_Mat_Mat_Mult(
			const LaGenMatComplex &  A,
			const LaGenMatComplex &  B,
			LaGenMatComplex &  C,
			bool  hermit_A,
			bool  hermit_B = false,
			LaComplex  alpha = 1.0,
			LaComplex  beta = 0.0
			)
		
		Perform the matrix-matrix operation C := alpha*A*B + beta*C
			 A and B are used in either non-hermitian or hermitian form depending on the function arguments.
				e.g. matrix transpose and complex conjugate
		
		Parameters:
			A		 - ???
			B		 - ???
			
			hermit_A - If true, use hermitian A
						i.e. A* (sometimes denoted conj(A')),
						the matrix transpose and complex conjugate, instead of A.
					 - If false, use A directly in non-hermitian form
			hermit_B - If true, use hermitian B
						i.e. B* (sometimes denoted conj(B')),
						the matrix transpose and complex conjugate, instead of B.
					 - If false, use B directly in non-hermitian form
		Internally this uses zgemm .
 ) 

	} */
}                                                        //empty
void matrix_inverse(LaGenMatComplex* matrix, int len){
    // http://lapackpp.sourceforge.net/html/laslv_8h.html#a042c82c5b818f54e7f000d068f14189
    LaVectorLongInt PIV = LaVectorLongInt(len);
    LUFactorizeIP(matrix, &PIV);
    LaLUInverseIP(matrix, &PIV);
}                       //untested

/* Main Program */
int main(){
	
    /* initialise stuff */
     matrix;
    
	/* generate the matrix */           //working
	int m = 10, x = 10;                                                     //m = matrix size, x = max
	COMPLEX* elements = new COMPLEX[m*m];                                   //initialise elements
	generate_elements(elements, m, x);                                      //generate random elements
	LaGenMatComplex* matrix = new LaGenMatComplex( elements, m, m, false ); //copy the elements to a matrix
    
    /* save and print the matrix */     //working
    LaGenMatComplex* initialMatrix;
    initialMatrix = matrix->copy();
    print_matrix(initialMatrix, "initial matrix");   
    
	/* calculate the eigenvalue matrix */
//    matrix_eignevalues(matrix);
    
    /* save and print the eigenvalue matrix */
//    LaGenMatComplex* eigenvalueMatrix;
//    eigenvalueMatrix = matrix->copy();
//    print_matrix(eigenvalueMatrix, "eigenvalue matrix");   
//    matrix = initialMatrix->copy();                      //reset matrix
    
	/* calculate the eigenvector matrix */
    
    /* save and print the eigenvector matrix */
//    LaGenMatComplex* eigenvectorMatrix;
//    eigenvectorMatrix = matrix->copy();
//    print_matrix(eigenvectorMatrix, "eigenvector matrix"); 
//    matrix = initialMatrix->copy();                      //reset matrix
    
	/* calculate the transpose eigenvector matrix */
    
    /* save and print the transpose eigenvector matrix */
//    LaGenMatComplex* tEigenvectorMatrix;
//    tEigenvectorMatrix = matrix->copy();
//    print_matrix(tEigenvectorMatrix, "transpose eigenvector matrix");
//    matrix = initialMatrix->copy();                      //reset matrix
    
	/* calculate the power matrix */
    
    /* save and print the power matrix */
    
	/* calculate the exponential matrix */
    
    /* save and print the exponential matrix */
    
    /* invert the matrix */
    matrix_inverse(matrix, m);
    
    /* save and print the inverse matrix */
    LaGenMatComplex* inverseMatrix;
    inverseMatrix = matrix->copy();
    print_matrix(inverseMatrix, "inverse matrix");   
    matrix = initialMatrix->copy();                      //reset matrix
    
	/* delete things on the heap: */
    delete [] elements;
	delete matrix;
	
}
