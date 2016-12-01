#include <iostream>
using namespace std;

int main(){
    int lattice [10][10];
    for(int i = 1; i < 11; i = i + 1){
        for(int j = 1; j < 11; j = j + 1){
            cout << "i = " << i << ", j = " << j << endl;
        }
    }
}
