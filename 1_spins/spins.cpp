#include <cstdlib>  //rand, srand
#include <iostream> //cout
#include <ctime>    //time(NULL)

using namespace std;

void random_seed(){
    srand(time(NULL));
}

void generate_lattice(int array[10][10]){
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            int ran = (rand() % 2)*2 - 1;
            array[i][j] = ran;
        }
    }
}

void print_lattice(int array[10][10]){
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            cout.width(3);
            cout << array[i][j];
        }
        cout << endl;
    }
}

int main(){
    int lattice [10][10];
    random_seed();
    generate_lattice(lattice);
    print_lattice(lattice);
}

/*
fix the randomisation using http://stackoverflow.com/questions/19665818/generate-random-numbers-using-c11-random-library
*/
