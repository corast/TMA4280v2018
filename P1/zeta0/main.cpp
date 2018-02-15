#include <cstdio> //printf function
#include <string> // stoi function (string to int)

#include "zeta0.h"

int main(int argc, char* argv[])
{
    //check if we got one argument
    if(argc != 2){
        std::printf("Please type one number n as argument to this program\n");
        return 1;
    }

    int n = std::stoi(argv[1],NULL,0);

    double pi = zeta(n);
    printf("Pi is aproximatly %f with n as %i \n", pi, n);
}