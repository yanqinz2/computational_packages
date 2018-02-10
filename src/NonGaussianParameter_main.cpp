#include <iostream>

#ifndef NONGAUSSIANPARAMETER_HPP
#define NONGAUSSIANPARAMETER_HPP
#include "NonGaussianParameter.hpp"
#endif

int main(int argc, char* argv[])
{ 
    std::cout << "\nComputing non-Gaussian parameter...\n" << std::endl;

    NonGaussianParameter non_gaussian_parameter_;

    non_gaussian_parameter_.read_commands(argc, argv);
    non_gaussian_parameter_.read_input_file();
    non_gaussian_parameter_.initialize_rest_members();
    non_gaussian_parameter_.compute_ngp();
    non_gaussian_parameter_.write_output_file();

    std::cout << "Successfully computed non-Gaussian parameter.\n" << std::endl;

    return 0;
}