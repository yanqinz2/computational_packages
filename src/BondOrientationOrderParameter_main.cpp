#include <iostream>

#ifndef BONDORIENTATIONORDERPRAMETER_HPP
#define BONDORIENTATIONORDERPRAMETER_HPP
#include "BondOrientationOrderParameter.hpp"
#endif

int main(int argc, char* argv[])
{
    std::cout << "\nComputing bond orientation order parameter of each layer...\n" << std::endl;

    BondOrientationOrderParameter bond_orientation_order_parameter_;

    bond_orientation_order_parameter_.read_commands(argc, argv);
    bond_orientation_order_parameter_.read_input_file();
    bond_orientation_order_parameter_.initialize_rest_members();
    bond_orientation_order_parameter_.compute_bop();
    bond_orientation_order_parameter_.write_output_file();

    std::cout << "Successfully computed bond orientation order parameter of each layer.\n" << std::endl;

    return 0;
}
