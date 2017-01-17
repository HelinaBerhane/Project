#include <random>   //random_device, mt19937
#include <iostream> //cout
#include <ctime>    //time(NULL)

using namespace std;

int rand(int a, int b){
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dist(a, b);
    return dist(gen);
}

void generate_lattice(int array[10][10]){
    int test = 0;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            array[i][j] = 2*rand(0, 1)-1;
            test += array[i][j];
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
    generate_lattice(lattice);
    print_lattice(lattice);

}

/*
fix the randomisation using http://stackoverflow.com/questions/19665818/generate-random-numbers-using-c11-random-library
*/
