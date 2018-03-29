#include <iostream>

#ifndef SELFVANHOVECORRELATIONFUNCTION_HPP
#define SELFVANHOVECORRELATIONFUNCTION_HPP
#include "SelfVanHoveCorrelationFunction.hpp"
#endif

int main(int argc, char* argv[])
{
    std::cout << "\nComputing self van Hove correlation function...\n" << std::endl;

    SelfVanHoveCorrelationFunction self_van_hove_correlation_function_;

    self_van_hove_correlation_function_.read_commands(argc, argv);
    self_van_hove_correlation_function_.read_input_file();
    self_van_hove_correlation_function_.initialize_rest_members();
    self_van_hove_correlation_function_.compute_gs_rt();
    self_van_hove_correlation_function_.write_output_file();

    std::cout << "Successfully computed self van Hove correlation function.\n" << std::endl;

    return 0;
}
