#include <iostream>

#include <Eigen/Dense>

int main(int argc, char* argv[]){
    std::cout << "Hello World" << std::endl;
    Eigen::Matrix<double, 3, 1> xyz(1.0, 2.0, 3.0);

    std::cout << xyz << std::endl;
    return true;    
}