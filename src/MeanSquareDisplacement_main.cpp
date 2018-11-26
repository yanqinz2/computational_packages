#include <iostream>

#ifndef MEANSQUAREDISPLACEMENT_HPP
#define MEANSQUAREDISPLACEMENT_HPP
#include "MeanSquareDisplacement.hpp"
#endif

int main(int argc, char* argv[])
{
    std::cout << "\nComputing mean square displacement...\n" << std::endl;

    MeanSquareDisplacement mean_square_displacement_;

    mean_square_displacement_.read_commands(argc, argv);
    mean_square_displacement_.read_input_file();
    mean_square_displacement_.initialize_rest_members();
    mean_square_displacement_.compute_msd();
    mean_square_displacement_.write_output_file();

    std::cout << "Successfully computed mean square displacement.\n" << std::endl;

    return 0;
}
