#include <iostream>
#include "zeta0.cpp"
using namespace std;
int main(int argc, char* argv[])
{
    //check if we got one argument
    if(argc < 1){
        printf("Please type one number n");
        return 1;
    }
    /*else if(argc >= 2){
        printf("Please type only one number n, %f", argc);
        return 2;
    } */

    int n = strtol(argv[1],NULL, 0);

    printf("n set as %i\n",n);
    double pi = zeta0(n);
    
    printf("Pi is aproximatly %f with %i itterations", pi, n);
}