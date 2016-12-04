#include <cstdlib>  //rand, srand
#include <iostream> //cout
#include <ctime>    //time(NULL)

using namespace std;

/*
void do_array_things( int array[10] ) {
    cout << array[0] << endl;
}
*/

int main(){
    srand(time(NULL));
    int lattice [10][10];
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            int ran = (rand() % 2)*2 - 1;
            lattice[i][j] = ran;
            cout.width(3);
            cout << lattice [i][j];
        }
        cout << endl;
    }
    /*
    int array[10];
    do_array_things( array );
    */
}

/*

make separate functions for each task and take them out of main

fix the randomisation using http://stackoverflow.com/questions/19665818/generate-random-numbers-using-c11-random-library

*/
