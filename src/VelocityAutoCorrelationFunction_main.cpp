#include <iostream>

#ifndef VELOCITYAUTOCORRELATIONFUNCTION_HPP
#define VELOCITYAUTOCORRELATIONFUNCTION_HPP
#include "VelocityAutoCorrelationFunction.hpp"
#endif

int main(int argc, char* argv[])
{
    std::cout << "\nComputing velocity auto correlation function...\n" << std::endl;

    VelocityAutoCorrelationFunction velocity_auto_correlation_function_;

    velocity_auto_correlation_function_.read_commands(argc, argv);
    velocity_auto_correlation_function_.read_input_file();
    velocity_auto_correlation_function_.initialize_rest_members();
    velocity_auto_correlation_function_.compute_vacf();
    velocity_auto_correlation_function_.write_output_file();

    std::cout << "Successfully computed velocity auto correlation function.\n" << std::endl;

    return 0;
}
