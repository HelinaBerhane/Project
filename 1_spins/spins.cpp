#include <iostream>
#include <string>
using namespace std;

int main(){
    string lattice [10][10];
    for(int i = 1; i < 11; i = i + 1){
        for(int j = 1; j < 11; j = j + 1){
            /*
            string a = "i";
            string b = to_string(i);
            string c = "j";
            string d = to_string(j);
            cout << a << b << c << d << endl;
            */
            lattice[i-1][j-1] = "i" + to_string(i) + "j" + to_string(j);
            cout << lattice [i-1][j-1] << endl;
        }
    }
    /*    
    int testint = 1;
    string test = "test";
    string teststr = test + to_string(testint);
    cout <<  teststr << endl;
    */
}
