#include <iostream>

#ifndef ATOMNUMBERDENSITY_HPP
#define ATOMNUMBERDENSITY_HPP
#include "AtomNumberDensity.hpp"
#endif

int main(int argc, char* argv[])
{
    std::cout << "\nComputing atom number density...\n" << std::endl;

    AtomNumberDensity atom_number_density_;

    atom_number_density_.read_commands(argc, argv);
    atom_number_density_.read_input_file();
    atom_number_density_.initialize_rest_members();
    atom_number_density_.compute_atom_number_density();
    atom_number_density_.write_output_file();

    std::cout << "Successfully computed atom number density.\n" << std::endl;

    return 0;
}
