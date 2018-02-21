/*
    Unit test, compare value of n = 3 to an pre computed value.
*/
#include <cmath>
#include "zeta0.h"
#include <cassert>
#include <string>

void unitTest();

int main(int argc, char* argv[]){
    unitTest();
    return 0;
}

void unitTest(){

    double pi = zeta(3);
    double calculatedPi = 2.857738;

    printf("Pi calculated is %f, pi computed is %f",calculatedPi,pi);

    //check if these are the same
    int tmp = pi * 6;
    double pi_2 = tmp/6.0;//3.¤¤¤¤¤¤
    double difference = std::abs(pi_2 - calculatedPi);
    printf("Difference is %f\n",difference);
    assert(difference < 0.001);
}