#include <iostream>     // cout
#include <stdlib.h>     // srand, rand
#include <time.h>       // time
using namespace std;

void random_seed(){
    srand(time(NULL));
}
int ran(int x){
    return rand() % x;
}
int fact(int x){
	if(x<=1){
        return 1;
	}else{
        return x*fact(x-1);
	}
}
long double powers(long long int x, long long int y){
    if(y<=1){
        return x;
    }else{
        y--;
        return x * powers(x, y);
    }
}

/* Redundant */
long double straight_exponential(long long int number, long long int iterations, int n){
    if(iterations == 1){
        return powers(number, n) / fact(n);
    }else{
        return powers(number, n) / fact(n) + straight_exponential(number, iterations-1, n+1);
    }
    /*
    float partial_result = 0;
    float result = 0;
    // over each iteration, calculate \sum_i{matrix^n/n!}
    for(float i = 0; i > 10; i++){
        partial_result = matrix;
        result = result + partial_result;
        cout << partial_result << endl;
        //multiply the resulting

        // To Do:
           // find out how long each operation takes
           // or find some way to decide how accurate you want to be
    }
    cout << result << endl;

    int y = 3;
    int powers_test = powers();
    return powers_test;
    */
}
long double exp_step(const long double number, const double step){
  double result;
  double A;
  double B;
  if(step <= 0){
    cout << "step (0): 1 " << endl;
    return number;
  }else{
      A = exp_step(number, step - 1);
      B = number / step;
      //cout << "step (" << step << "): ";
      //cout << number << " / " << step << " = " << B << endl;
      result = B * A;
      //cout << B << " * " << A << " = " << result << endl;
      return result;
  }
}
long double recursive_exponential(const long double number, int iterations){
  double result = 0;
  long double step;
  cout << result;
  for(int i = 0; i < iterations; i ++){
    step = exp_step(number, i);
    cout << result << " + " << step;
    result += step;
    cout << " = " << result << endl;
  }
  return result;
}

/* Working Version */
double exponent(const double number, const int iterations){
    double A, B = 1, C = 1;
    for(int n = 1; n <= iterations; n++){
        for(int i = 1; i <= n; i++){
            A = number / i;
            B *= A;
        }
        C += B;
        B = 1;
    }
    return C;
}

int main(){
    int iterations = 1000, n = 1;
    double number = 5;
    cout.precision(50);
    cout << "better: " << exponent(number, iterations) << endl;
}
