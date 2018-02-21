/*
    Unit test, compare value of n = 3 to an pre computed value.
*/
#include <cmath>
#include "mach0.h"
#include <cassert>

void unitTest();

int main(int argc, char* argv[]){
    unitTest();
    return 0;
}

void unitTest(){

    double pi = mach(3);
    double computedPi = 2.857738;
    //check if these are the same
    int tmp = pi * 6;
    double pi_2 = tmp/6.0;//3.¤¤¤¤¤¤
    double difference = std::abs(pi_2 - computedPi);
    assert(difference < 0.001);
}