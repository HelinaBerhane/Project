#include <random>   //random_device, mt19937
#include <cstdlib>	//rand, sran
#include <iostream> //cout
#include <ctime>    //time(NULL)
#include <iostream>
// #include <iomanip>
#include <string>
// #include <map>
#include <math.h>


using namespace std;

// random_device rd;
// mt19937 gen(rd());

/* output */
void print_lattice(int array[10][10]){
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            cout.width(3);
            cout << array[i][j];
        }
        cout << endl;
    }
}

/* random number generation */
// int basic_random_int(int max_rand){
//     return rand() % (max_rand+1);
// }
// float basic_random_float(float steps, float max_rand){
//     return basic_random_int(max_rand*steps)/(steps);
// }
// int random_int(int min, int max){
//     uniform_int_distribution<> dist(min, max);
//     return dist(gen);
// }
// int random_float(float a, float b, string hi){
//     // random_device rd;
//     // mt19937 gen(rd());
//     // uniform_int_distribution<> dist(a, b);
//     // return dist(gen);
//     random_device rd;
//     mt19937 e2(rd());
//     uniform_real_distribution<> dist(a, b);
//     if(hi == "y"){
//         map<int, int> hist;
//         for (int n = 0; n < 10000; ++n) {
//             ++hist[floor(dist(e2))];
//         }
//         for (auto p : hist) {
//             cout << fixed << setprecision(1) << setw(2)
//                       << p.first << ' ' << string(p.second/200, '*') << '\n';
//         }
//         return 0;
//     }else{
//         return dist(e2);
//     }
// }
float random_float(float min, float max){
    random_device rd;                                       //maybe find a way to
    mt19937 gen(rd());                                      //not seed ebvery time?
    std::uniform_real_distribution<> dis(min, max);         //but cba
    return dis(gen);
    // random_device rd;                                   //seed
    // mt19937 mt(rd());                                   //generate numbers (0 to 2^32)
    // uniform_real_distribution<float> dist(min, max+1);  //set distribution
    // return dist(mt);                                    //run the no
    //info
        //see https://channel9.msdn.com/Events/GoingNative/2013/rand-Considered-Harmful
        //speed - 499 Mb/s = 6.5 cycles per byte (see 16:20)
        //has a period of 2^19937, so will pretty much never repeat (see 16:50)
        //seedable, so you can replicate a run (could be useful)
        //reproducable (17:45)
        //you can't have multiple threads acting on the same generator at the same tim(26:0)
}

/* lattice generation */
// void generate_lattice(int array[10][10]){
//     int test = 0;
//     for(int i = 0; i < 10; i++){
//         for(int j = 0; j < 10; j++){
//             array[i][j] = 2*random_int(0, 1)-1;
//             test += array[i][j];
//         }
//     }
// }

/* calculations */
float lambda(const float U){
    return acoshf(exp(sqrt(0.125*U)/2));
}

/* testing */
// void test_basic_random_int(){
//     int max_rand = 9, iterations = 1000;
//     for(int i = 0; i < iterations; i++){
//         cout << basic_random_int(max_rand) << endl;
//     }
// }
// void test_basic_random_float(){
//     int iterations = 1000;
//     float steps = 10000, max_rand = 9;
//     for(int i = 0; i < iterations; i++){
//         cout << basic_random_float(steps, max_rand) << endl;
//     }
// }
void test_random_float(){
    int count = 100, min = 0, max = 10;
    // random_device rd;
    // mt19937 gen(rd());
    // uniform_real_distribution<> dist(min, max);
    for (int i = 0; i < count; i++) {
        // cout << dist(gen) << endl;
        cout << random_float(min, max) << endl;
    }
}

int main(){
    // int lattice [10][10];
    // generate_lattice(lattice);
    // print_lattice(lattice);
    cout << lambda(1) << endl;
}

/*
fix the randomisation using http://stackoverflow.com/questions/19665818/generate-random-numbers-using-c11-random-library
*/
