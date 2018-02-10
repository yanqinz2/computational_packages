#include <iostream>

#ifndef UNWRAPTRAJECTORY_HPP
#define UNWRAPTRAJECTORY_HPP
#include "UnwrapTrajectory.hpp"
#endif

int main(int argc, char* argv[])
{
    std::cout << "\nUnwrapping trajectory...\n" << std::endl;

    UnwrapTrajectory unwrap_trajectory_;

    unwrap_trajectory_.read_commands(argc, argv);
    unwrap_trajectory_.initialize_rest_members();
    unwrap_trajectory_.unwrap_trajectory();

    std::cout << "Successfully unwrapped the trajectory.\n" << std::endl;

    return 0;
}