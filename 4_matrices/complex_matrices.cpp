#define LA_COMPLEX_SUPPORT

#include <iostream> //cout
#include <cstdlib>	//rand, sran
#include <gmc.h> 	//LaGenMatComplex
//#include <f2c.h>    //complex type
//#include <complex>

using namespace std;

/* Testing */
/*
void print_matrix(complex* elements, int m) {				//out of date
	for( int i = 0; i = m*m; i++ ) {
		cout << elements[i] << " ";
		if( i%m == 0 ) {
			cout << endl;
		}
	}
}
*/
void print_states(LaGenMatComplex* matrix){							//not sure wat does
	cout << matrix->info();
}
void ideal_print_matrix(LaGenMatComplex* matrix){					//untested
	 /* documentation{
		ostream documentation{
			std::ostream & operator <<  ( std::ostream & , const LaGenMatComplex & ) [friend]
			
			Print the matrix to the given output stream.
			If the matrix info flag is set, then this prints only the matrix info
				see LaGenMatDouble::info(). 
			Otherwise all matrix elements are printed.
		}
		setPrintFormat documentation{
			static void LaPreferences::setPrintFormat ( pFormat  p, bool newlines = true ) [static] 
			Set the output display format. The default is LaPreferences::NORMAL.
			
			Parameters:
				p		 - The preferred output format
				newlines - Toggles multiline matrix output 
					true = place a newline after each matrix row
					false = use only the appropriate MATLAB/MAPLE delimiter
					false = use only the appropriate MATLAB/MAPLE delimiter
			
			The following is how to modify the output display format to be compatible with your favourite math program:
				Place " #include LA_PREFS_H " in the include statements, somewhere after "lafnames.h".
				
				At the beginning of your code, call e.g. " LaPreferences::setPrintFormat( LaPreferences::MATLAB, true); "
		}
	} */
	cout << *matrix << endl;
}

/* Number generation */
void generate_elements(COMPLEX* elements, int m, int x){	//untested but should work
	srand(time(NULL));				//seeds
    for(int i = 0; i < m*m; i++){
        elements[i].r = rand() % x;	//1-x
        elements[i].i = rand() % x;
	}
}

/* Matrix manipulation */

//void matrix_powers(){
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
//}



//void solve_eignevalues(LaGenMatComplex* matrix){
	/* see laslv.h */
//}


/* Main Program */
int main(){  											//untested
	
	/* generate the matrix */
	int m = 10, x = 10; 	//m = matrix size, x = max element size
	COMPLEX* elements = new COMPLEX[m*m];
	LaGenMatComplex diag;
	generate_elements(elements, m, x);
	LaGenMatComplex* matrix = new LaGenMatComplex( elements, m, m, false );
	diag = matrix->diag();
	
	/* test */
	//print_matrix( elements, m);
	//print_states(matrix);
	ideal_print_matrix(matrix);
    ideal_print_matrix(&diag);
	
	/* diagonalise the matrix */
	/* when the program is done, we need to delete the things on the heap: */
	delete [] elements;
	delete matrix;
	
}
