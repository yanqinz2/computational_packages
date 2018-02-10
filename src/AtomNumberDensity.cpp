#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifndef ATOMNUMBERDENSITY_HPP
#define ATOMNUMBERDENSITY_HPP
#include "AtomNumberDensity.hpp"
#endif

#ifndef OMP
#define OMP
#include "omp.h"
#endif

AtomNumberDensity::AtomNumberDensity():
        print_headline_(false),
        trajectory_file_name_(""),
        output_file_name_(""),
        start_frame_(1),
        end_frame_(1),
        dimension_(3),
        number_of_layers_(1),
        output_precision_(15)
{
}

AtomNumberDensity::~AtomNumberDensity()
{
}

void AtomNumberDensity::read_commands(int argc, char* argv[])
{
    for (int input = 1; input < argc; ++input) {
        if (strcmp(argv[input], "-in") == 0) {
            input_file_name_ = argv[++input];
            continue;
        }

        std::cerr << "ERROR: Unrecognized flag '" << argv[input] << "' from command inputs." << std::endl;
        exit(1);
    }
}

void AtomNumberDensity::read_input_file()
{
    std::ifstream input_file_(input_file_name_);

    if (!input_file_) {
        std::cerr << "ERROR: Input file not found!\n" << std::endl;
        exit(1);
    }

    std::string input_word_;

    // check for comment
    while (input_file_ >> input_word_) {
        if (input_word_[0] == '#') {
            getline(input_file_, input_word_);
            continue;
        }

        if (input_word_ == "print_headline") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            if (input_word_ == "true" || input_word_ == "yes") {
                print_headline_ = true;
            }
            else if(input_word_ == "false" || input_word_ == "no") {
                print_headline_ = false;
            }
            else {
                print_headline_ = stoi(input_word_);
            }
            continue;
        }
        if (input_word_ == "trajectory_file_path") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            trajectory_file_name_ = input_word_;
            continue;
        }
        if (input_word_ == "output_file_path") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            output_file_name_ = input_word_;
            continue;
        }
        if (input_word_ == "start_frame") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            start_frame_ = stod(input_word_);
            continue;
        }
        if (input_word_ == "end_frame") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            end_frame_ = stod(input_word_);
            continue;
        }
        if (input_word_ == "number_of_layers") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            number_of_layers_ = stod(input_word_);
            continue;
        }
        if (input_word_ == "output_precision") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            output_precision_ = stod(input_word_);
            continue;
        }
    }

    input_file_.close();
}

void AtomNumberDensity::initialize_rest_members()
{
    std::cout << "Initializing...\n" << std::endl;

    if (start_frame_ == 0) {
        std::cerr << "ERROR: Start frame should be positive!\n" << std::endl;
        exit(1);
    }
    if (end_frame_ < start_frame_) {
        std::cerr << "ERROR: End frame should be equal or greater than start frame!\n" << std::endl;
        exit(1);
    }
    if (number_of_layers_ == 0) {
        std::cerr << "ERROR: Number of layers should be positive!\n" << std::endl;
        exit(1);
    }

    number_of_frames_ = end_frame_ - start_frame_ + 1;

    left_point_of_box_.resize(dimension_);
    right_point_of_box_.resize(dimension_);
    box_length_.resize(dimension_);

    std::ifstream trajectory_file_(trajectory_file_name_);

    if (!trajectory_file_) {
        std::cerr << "ERROR: Trajectory file not found!\n" << std::endl;
        exit(1);
    }

    std::string line_stream_;
    unsigned int atom_type_temp_;

    for (int i_line = 0; i_line < 3; ++i_line) {
        getline(trajectory_file_, line_stream_); // Skip the first line
    }
    trajectory_file_ >> number_of_atoms_;
    getline(trajectory_file_, line_stream_);
    for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
        getline(trajectory_file_, line_stream_);
        trajectory_file_ >> left_point_of_box_[i_dimension] >> right_point_of_box_[i_dimension];
        box_length_[i_dimension] = right_point_of_box_[i_dimension] - left_point_of_box_[i_dimension];
    }
    for (int i = 0; i < number_of_atoms_ + 1; ++i) {
        getline(trajectory_file_, line_stream_);
    }
    trajectory_file_ >> line_stream_ >> number_of_types_of_atoms_;

    trajectory_file_.seekg(0, trajectory_file_.beg);
    for (int i = 0; i < dimension_ + 6; ++i) {
        getline(trajectory_file_, line_stream_);
    }
    number_of_atoms_of_each_type_.resize(number_of_types_of_atoms_);
    for (int i_type = 0; i_type < number_of_types_of_atoms_; ++i_type) {
        number_of_atoms_of_each_type_[i_type] = 0;
    }
    for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
        trajectory_file_ >> line_stream_ >> atom_type_temp_;
        ++number_of_atoms_of_each_type_[atom_type_temp_ - 1];
        getline(trajectory_file_, line_stream_);
    }

    trajectory_file_.close();

    width_of_layers_ = right_point_of_box_[2] / number_of_layers_;

    atom_type_.resize(number_of_atoms_);

    trajectory_.resize(number_of_atoms_);
    for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom){
        trajectory_[i_atom].resize(dimension_);
    }

    number_density_of_each_layer_.resize(number_of_layers_);
    for (int i_layer = 0; i_layer < number_of_layers_; ++i_layer) {
        number_density_of_each_layer_[i_layer].resize(number_of_types_of_atoms_ + 1);
    }
    for (int i_layer = 0; i_layer < number_of_layers_; ++i_layer){
        for (int i_type = 0; i_type < number_of_types_of_atoms_ + 1; ++i_type) {
            number_density_of_each_layer_[i_layer][i_type] = 0.0;
        }
    }
}

void AtomNumberDensity::compute_atom_number_density()
{
    std::cout << "Computing...\n" <<std::endl;

    std::ifstream trajectory_file_(trajectory_file_name_);
    std::string line_stream_;

    // jump over needless frames
    for (int i_frame = 1; i_frame < start_frame_; ++i_frame) {
        for (int i_line = 0; i_line < number_of_atoms_ + dimension_ + 6; ++i_line) {
            getline(trajectory_file_, line_stream_);
        }
    }

    for (int i_frame = 0; i_frame < number_of_frames_; ++i_frame) {
        // reading the trajectory
        for (int i_line = 0; i_line < dimension_ + 6; ++i_line) {
            getline(trajectory_file_, line_stream_);
        }
        for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
            trajectory_file_ >> line_stream_ >> atom_type_[i_atom];
            for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                trajectory_file_ >> trajectory_[i_atom][i_dimension];
            }
            getline(trajectory_file_, line_stream_);
        }

#pragma omp parallel for
        for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
            --atom_type_[i_atom];
            for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                trajectory_[i_atom][i_dimension] =
                        trajectory_[i_atom][i_dimension] * box_length_[i_dimension] + left_point_of_box_[i_dimension];
            }
        }

        for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
            int temp_atom_layer_ = int(fabs(trajectory_[i_atom][2]) / width_of_layers_);
            if (temp_atom_layer_ >= number_of_layers_) {
                ++number_density_of_each_layer_[number_of_layers_ - 1][atom_type_[i_atom]];
                ++number_density_of_each_layer_[number_of_layers_ - 1][number_of_types_of_atoms_];           
            }
            else {
                ++number_density_of_each_layer_[temp_atom_layer_][atom_type_[i_atom]];
                ++number_density_of_each_layer_[temp_atom_layer_][number_of_types_of_atoms_];
            }
        }

        std::cout << "\r" << float(i_frame + 1) / number_of_frames_ * 100 << "% has been completed..." << std::flush;
    }

    std::cout << "\n" << std::endl;

    trajectory_file_.close();

    double temp_normalization_ = 2 * number_of_frames_ * width_of_layers_ / box_length_[0];
#pragma omp parallel for
    for (int i_layer = 0; i_layer < number_of_layers_; ++i_layer) {
        for (int i_type = 0; i_type < number_of_types_of_atoms_; ++i_type) {
            number_density_of_each_layer_[i_layer][i_type] /= temp_normalization_ * number_of_atoms_of_each_type_[i_type];
        }
        number_density_of_each_layer_[i_layer][number_of_types_of_atoms_] /= temp_normalization_ * number_of_atoms_;
    }
}

void AtomNumberDensity::write_output_file()
{
    std::cout << "Writing output file...\n" << std::endl;

    std::ofstream number_density_output_file_(output_file_name_);
    number_density_output_file_.precision(output_precision_);

    if (print_headline_) {
        number_density_output_file_ << "z     ";
        for (int i_type = 0; i_type < number_of_types_of_atoms_; ++i_type) {
            number_density_output_file_ << "atom_" << i_type + 1 << "     ";
        }
        number_density_output_file_ << "total_density\n";
    }

    for (int i_layer = 0; i_layer < number_of_layers_; ++i_layer) {
        number_density_output_file_ << (i_layer + 0.5) * width_of_layers_ << " ";
        for (int i_type = 0; i_type < number_of_types_of_atoms_; ++i_type) {
            number_density_output_file_ << number_density_of_each_layer_[i_layer][i_type] << " ";
        }
        number_density_output_file_ << number_density_of_each_layer_[i_layer][number_of_types_of_atoms_] << "\n";
    }

    number_density_output_file_.close();
}