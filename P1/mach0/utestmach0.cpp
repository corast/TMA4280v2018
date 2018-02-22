/*
    Unit test, compare value of n = 3 to an pre computed value.
*/
#include <cmath>
#include "mach0.h"
#include <cassert>
#include <string>

void unitTest();

int main(int argc, char* argv[]){
    unitTest();
    return 0;
}

void unitTest(){

    double pi = mach(3);
    double computedPi = 3.141621;
    //check if these are the same
    //printf("pi = %f, cpi = %f",pi,computedPi);
    double difference = std::abs(pi - computedPi);
    assert(difference < 0.00001);
    printf("unit test PASSED\n");    
}