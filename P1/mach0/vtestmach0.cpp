#include <cmath> //for pow function
#include "mach0.h" //for zeta function
//#include <cassert>
#include <string> //for printf
#include <fstream> //for file handing: creation, opening, writing etc.
#include <iomanip>

int validTest();
//double can relibly hold 15 decimal places.
//const double PI = 3.141592653589793;
const double PI = 4.0 * atan(1.0);
int main(int argc, char* argv[]){
    
    if(!validTest()){
        printf("Validation PASSED\n");
    }
    return 0;
}

int validTest(){
    //create file to write error for each k value.
    std::ofstream myfile;
    myfile.open("error.txt");
    int n = 0;
    double pi;
    double error_old;
    double error;

    //calculate error for each k, and vertify that it converges to 0 towards k=24.
    for(int k = 1; k<=24;k++){
        n = 1 << k; //binary bitshit to the right, 2 bits.
        pi = mach(n);
        error = std::abs(PI-pi);
        
        myfile << "For n = " << std::setw(10) << n << "\t\t Error: " << error << std::endl;
        if(error > error_old && k != 1){
            //means we are not convering to a better value than before, we can stop test
            printf("not convering to pi!\n");
            printf("error:%f , error_old:%f\n",error, error_old);
            return 1;
        }
        error_old = error;
    }
    myfile.close();
    return 0;
}
