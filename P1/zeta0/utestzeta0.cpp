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

    //printf("Pi calculated is %f, pi computed is %f\n",calculatedPi,pi);

    //check if these are the same
    double difference = std::abs(pi - calculatedPi);
    assert(difference < 0.00001);
    printf("unit test PASSED\n");
}