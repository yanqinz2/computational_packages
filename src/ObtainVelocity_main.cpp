#include <iostream>

#ifndef OBTAINVELOCITY_HPP
#define OBTAINVELOCITY_HPP
#include "ObtainVelocity.hpp"
#endif

int main(int argc, char* argv[])
{
    std::cout << "\nComputing velocity...\n" << std::endl;

    ObtainVelocity obtain_velocity_;

    obtain_velocity_.read_commands(argc, argv);
    obtain_velocity_.initialize_rest_members();
    obtain_velocity_.compute_and_write_output_file();

    std::cout << "Successfully computed the velocity.\n" << std::endl;

    return 0;
}