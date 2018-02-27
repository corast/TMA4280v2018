#include <cmath> //for pow function
#include "zeta0.h" //for zeta function
//#include <cassert>
#include <string> //for printf
#include <fstream> //for file handing: creation, opening, writing etc.
#include <iomanip>

int validTest();
//double can relibly hold 15 decimal places.
const double PI = 3.141592653589793;

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
        n = pow(2,k);
        pi = zeta(n); //calculate pi
        error = std::abs(PI-pi); //calculate error between calculate and costant pi.
        std::string line;
        //line = "For n = " + std::to_string(n) + "\t\tError: " + std::to_string(error) + "\n";

        //output to file.
        myfile << "For n = " << std::setw(10) << n << "\t\t Error: " << error << std::endl;
        if(error > error_old && k != 1){ //check if our new error, is better than the old value, or else we are not converging.
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