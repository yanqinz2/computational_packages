#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifndef VELOCITYAUTOCORRELATIONFUNCTION_HPP
#define VELOCITYAUTOCORRELATIONFUNCTION_HPP
#include "VelocityAutoCorrelationFunction.hpp"
#endif

#ifndef OMP
#define OMP
#include "omp.h"
#endif

VelocityAutoCorrelationFunction::VelocityAutoCorrelationFunction():
        print_headline_(false),
        trajectory_file_name_(""),
        output_file_name_(""),
        dimension_(3),
        calculation_dimension_(123),
        start_frame_(1),
        interval_(1),
        number_of_timepoints_(1),
        number_of_frames_to_be_averaged_(1),
        number_of_layers_(1),
        output_precision_(15),
        frame_interval_time_(0.005),
        layer_left_point_(0.0),
        layer_right_point_(0.0)
{
}

VelocityAutoCorrelationFunction::~VelocityAutoCorrelationFunction()
{
}

void VelocityAutoCorrelationFunction::read_commands(int argc, char* argv[])
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

void VelocityAutoCorrelationFunction::read_input_file()
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
        if (input_word_ == "calculation_dimension") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            calculation_dimension_ = stod(input_word_);
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
        if (input_word_ == "interval") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            interval_ = stod(input_word_);
            continue;
        }
        if (input_word_ == "number_of_timepoints") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            number_of_timepoints_ = stod(input_word_);
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
        if (input_word_ == "number_of_frames_to_be_averaged") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            number_of_frames_to_be_averaged_ = stod(input_word_);
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
        if (input_word_ == "frame_interval_time") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            frame_interval_time_ = stod(input_word_);
            continue;
        }
        if (input_word_ == "left_point_of_the_layer") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            layer_left_point_ = stod(input_word_);
            continue;
        }
        if (input_word_ == "right_point_of_the_layer") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            layer_right_point_ = stod(input_word_);
            continue;
        }
    }
    input_file_.close();
}

void VelocityAutoCorrelationFunction::initialize_rest_members()
{
    std::cout << "Initializing...\n" << std::endl;

    // check for legality
    if (calculation_dimension_ != 1 && calculation_dimension_ != 2 && calculation_dimension_ != 3 &&
        calculation_dimension_ != 12 && calculation_dimension_ != 13 && calculation_dimension_ != 23 && calculation_dimension_ != 123) {
        std::cerr << "ERROR: Calculation dimension can only be 1, 2, 3, 12, 13, 23 or 123!" << std::endl;
    }
    if (start_frame_ == 0) {
        std::cerr << "ERROR: Start frame should be positive!\n" << std::endl;
        exit(1);
    }
    if (interval_ == 0) {
        std::cerr << "ERROR: Interval should be positive!\n" << std::endl;
        exit(1);
    }
    if (number_of_timepoints_ == 0) {
        std::cerr << "ERROR: Number of timepoints should be positive!\n" << std::endl;
        exit(1);
    }
    if (number_of_layers_ == 0) {
        std::cerr << "ERROR: Number of layers should be positive!\n" << std::endl;
        exit(1);
    }
    if (number_of_frames_to_be_averaged_ == 0) {
        std::cerr << "ERROR: Number of averaged frames should be positive!\n" << std::endl;
        exit(1);
    }
    if (frame_interval_time_ < 0) {
        std::cerr << "ERROR: Frame interval time should be positive!\n" << std::endl;
        exit(1);
    }
    if (layer_left_point_ < 0) {
        std::cerr << "ERROR: Left point of the layer should be zero or positive!\n" << std::endl;
        exit(1);
    }
    if (layer_right_point_ - layer_left_point_ < 0) {
        std::cerr << "ERROR: Right point of the layer should be larger than the left point!\n" << std::endl;
        exit(1);
    }

    layer_width_ = (layer_right_point_ - layer_left_point_) / number_of_layers_;

    std::ifstream trajectory_file_(trajectory_file_name_);

    if (!trajectory_file_) {
        std::cerr << "ERROR: Trajectory file not found!\n" << std::endl;
        exit(1);
    }

    std::string line_stream_;

    for (int i_line = 0; i_line < 3; ++i_line) {
        getline(trajectory_file_, line_stream_);
    }
    trajectory_file_ >> number_of_atoms_;
    for (int i_line = 0; i_line < dimension_ + 3; ++i_line) {
        getline(trajectory_file_, line_stream_);
    }
    
    number_of_types_of_atoms_ = 0;
    atom_type_.resize(number_of_atoms_);

    for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
        trajectory_file_ >> line_stream_ >> atom_type_[i_atom];
        if (atom_type_[i_atom] > number_of_types_of_atoms_) {
            number_of_types_of_atoms_ = atom_type_[i_atom];
        }
        getline(trajectory_file_, line_stream_);
    }
    trajectory_file_.seekg(0, std::ios::beg);

    number_of_frames_ = interval_ * number_of_timepoints_ + number_of_frames_to_be_averaged_;

    time_table_.resize(number_of_timepoints_ + 1);
#pragma omp parallel for
    for (int i_timepoint = 0; i_timepoint <= number_of_timepoints_; ++i_timepoint) {
        time_table_[i_timepoint] = i_timepoint * interval_ * frame_interval_time_;
    }

    atom_layer_.resize(number_of_frames_);
    for (int i_frame = 0; i_frame < number_of_frames_; ++i_frame) {
        atom_layer_[i_frame].resize(number_of_atoms_);
    }

    atom_velocities_.resize(number_of_frames_);
    for (int i_frame = 0; i_frame < number_of_frames_; ++i_frame) {
        atom_velocities_[i_frame].resize(number_of_atoms_);
        for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
            atom_velocities_[i_frame][i_atom].resize(dimension_);
        }
    }

    number_of_atoms_for_average_.resize(number_of_layers_);
    for (int i_layer = 0; i_layer < number_of_layers_; ++i_layer){
        number_of_atoms_for_average_[i_layer].resize(number_of_types_of_atoms_);
        for (int i_type = 0; i_type < number_of_types_of_atoms_; ++i_type) {
            number_of_atoms_for_average_[i_layer][i_type].resize(number_of_timepoints_ + 1);
        }
    }

    velocity_auto_correlation_function_.resize(number_of_layers_);
    for (int i_layer = 0; i_layer < number_of_layers_; ++i_layer) {
        velocity_auto_correlation_function_[i_layer].resize(number_of_types_of_atoms_);
        for (int i_type = 0; i_type < number_of_types_of_atoms_; ++i_type) {
            velocity_auto_correlation_function_[i_layer][i_type].resize(number_of_timepoints_ + 1);
        }
    }

#pragma omp parallel for
    for (int i_layer = 0; i_layer < number_of_layers_; ++i_layer) {
        for (int i_type = 0; i_type < number_of_types_of_atoms_; ++i_type) {
            for (int i_timepoint = 0; i_timepoint <= number_of_timepoints_ + 1; ++i_timepoint) {
                number_of_atoms_for_average_[i_layer][i_type][i_timepoint] = 0;
                velocity_auto_correlation_function_[i_layer][i_type][i_timepoint] = 0.0;
            }
        }
    }

    for (int i_frame = 1; i_frame < start_frame_; ++i_frame) {
        for (int i_line = 0; i_line < number_of_atoms_ + dimension_ + 6; ++i_line) {
            getline(trajectory_file_, line_stream_);
        }
    }

    for (int i_frame = 0; i_frame < number_of_frames_; ++i_frame) {
        for (int i_line = 0; i_line < dimension_ + 6; ++i_line) {
            getline(trajectory_file_, line_stream_);
        }
        for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
            double temp_z_coordinate_;
            trajectory_file_ >> line_stream_ >> line_stream_;
            for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                trajectory_file_ >> atom_velocities_[i_frame][i_atom][i_dimension];
            }
            trajectory_file_ >> temp_z_coordinate_;
            atom_layer_[i_frame][i_atom] = int(floor(fabs(temp_z_coordinate_) - layer_left_point_) / layer_width_);
            getline(trajectory_file_, line_stream_);
        }
    }

    trajectory_file_.close();
}

void VelocityAutoCorrelationFunction::compute_vacf()
{
    std::cout << "Computing...\n" << std::endl;

    // Determine the dimension needed to be calculated if only 2 dimensions neeeded
    int dimension_0_, dimension_1_, dimension_2_;
    if (calculation_dimension_ < 10) {
        dimension_0_ = calculation_dimension_ - 1;
    }
    else if (calculation_dimension_ < 100) {
        dimension_1_ = calculation_dimension_ / 10 - 1;
        dimension_2_ = calculation_dimension_ % 10 - 1;
    }

#pragma omp parallel for
    for (int i_average = 0; i_average < number_of_frames_to_be_averaged_; ++i_average) {
        for (int i_timepoint = 0; i_timepoint <= number_of_timepoints_; ++i_timepoint) {
            for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
                if (atom_layer_[i_average][i_atom] == atom_layer_[i_average + i_timepoint * interval_][i_atom] &&
                    atom_layer_[i_average][i_atom] >= 0 && atom_layer_[i_average][i_atom] < number_of_layers_) {
                    double temp_addition_ = 0.0;
                    if (calculation_dimension_ < 10) {
                        temp_addition_ += atom_velocities_[i_average][i_atom][dimension_0_] *
                                          atom_velocities_[i_average + i_timepoint * interval_][i_atom][dimension_0_];
                    }
                    else if (calculation_dimension_ > 100) {
                        for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                            temp_addition_ += atom_velocities_[i_average][i_atom][i_dimension] *
                                              atom_velocities_[i_average + i_timepoint * interval_][i_atom][i_dimension];
                        }
                    }
                    else {
                        temp_addition_ += atom_velocities_[i_average][i_atom][dimension_1_] *
                                          atom_velocities_[i_average + i_timepoint * interval_][i_atom][dimension_1_] +
                                          atom_velocities_[i_average][i_atom][dimension_2_] *
                                          atom_velocities_[i_average + i_timepoint * interval_][i_atom][dimension_2_];
                    }
#pragma omp atomic
                    velocity_auto_correlation_function_[atom_layer_[i_average][i_atom]][atom_type_[i_atom] - 1][i_timepoint] += temp_addition_;
#pragma omp atomic
                    ++number_of_atoms_for_average_[atom_layer_[i_average][i_atom]][atom_type_[i_atom] - 1][i_timepoint];
                }
            }
        }
    }

#pragma omp parallel for
    for (int i_layer = 0; i_layer < number_of_layers_; ++i_layer) {
        for (int i_type = 0; i_type < number_of_types_of_atoms_; ++i_type) {
            for (int i_timepoint = 0; i_timepoint <= number_of_timepoints_; ++i_timepoint) {
                velocity_auto_correlation_function_[i_layer][i_type][i_timepoint] /= number_of_atoms_for_average_[i_layer][i_type][i_timepoint];
            }
        }
    }
}

void VelocityAutoCorrelationFunction::write_output_file()
{
    std::cout << "Writing output file...\n" << std::endl;

    std::ofstream output_file_(output_file_name_);
    output_file_.precision(output_precision_);

    for (int i_layer = 0; i_layer < number_of_layers_; ++i_layer) {
        if (print_headline_) {
            output_file_ << "Layer range: " << layer_left_point_ + i_layer * layer_width_ << " - "
                             << layer_left_point_ + (i_layer + 1) * layer_width_  << "\n";
            output_file_ << "time";
            for (int i_type = 0; i_type < number_of_types_of_atoms_; ++i_type) {
                output_file_ << " atom_" << i_type + 1;
            }
            output_file_ << "\n";
        }

        for (int i_timepoint = 0; i_timepoint <= number_of_timepoints_; ++i_timepoint) {
            output_file_ << time_table_[i_timepoint];
            for (int i_type = 0; i_type < number_of_types_of_atoms_; ++i_type) {
                output_file_ << " " << velocity_auto_correlation_function_[i_layer][i_type][i_timepoint];
            }
            output_file_ << "\n";
        }
    }

    output_file_.close();
}
