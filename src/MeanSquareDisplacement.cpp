#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifndef MEANSQUAREDISPLACEMENT_HPP
#define MEANSQUAREDISPLACEMENT_HPP
#include "MeanSquareDisplacement.hpp"
#endif

#ifndef OMP
#define OMP
#include "omp.h"
#endif

MeanSquareDisplacement::MeanSquareDisplacement():
        input_file_name_("r2_t.in"),
        trajectory_file_name_(""),
        output_file_name_("r2_t.txt"),
        time_scale_type_("linear"),
        print_headline_(false),
        is_trajectory_wrapped_(true),
        dimension_(3),
        calculation_dimension_(123),
        start_frame_(1),
        number_of_timepoints_(1),
        number_of_frames_to_be_averaged_(1),
        number_of_layers_(1),
        output_precision_(15),
        number_of_atoms_(0),
        number_of_types_of_atoms_(0),
        number_of_frames_(0),
        interval_(1.0),
        frame_interval_time_(0.001),
        layer_left_point_(0.0),
        layer_right_point_(0.0),
        layer_width_(0.0),
        time_table_({}),
        atom_type_({}),
        atom_layer_({}),
        atom_coordinates_({}),
        mean_square_displacement_({})
{
}

MeanSquareDisplacement::~MeanSquareDisplacement()
{
}

void MeanSquareDisplacement::read_commands(int argc, char **argv)
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

void MeanSquareDisplacement::read_input_file()
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
        if (input_word_ == "is_trajectory_wrapped") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            if (input_word_ == "true" || input_word_ == "yes") {
                is_trajectory_wrapped_ = true;
            }
            else if(input_word_ == "false" || input_word_ == "no") {
                is_trajectory_wrapped_ = false;
            }
            else {
                is_trajectory_wrapped_ = stoi(input_word_);
            }
            continue;
        }
        if (input_word_ == "time_scale_type") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            if (input_word_ == "linear" || input_word_ == "log") {
                time_scale_type_ = input_word_;
            }
            else {
                std::cerr << "ERROR: Cannot analyze time_scale_type: " << input_word_ << "!\n" << std::endl;
                exit(1);
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
            calculation_dimension_ = stoi(input_word_);
            continue;
        }
        if (input_word_ == "start_frame") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            start_frame_ = stoi(input_word_);
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
            number_of_timepoints_ = stoi(input_word_);
            continue;
        }
        if (input_word_ == "number_of_layers") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            number_of_layers_ = stoi(input_word_);
            continue;
        }
        if (input_word_ == "number_of_frames_to_be_averaged") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            number_of_frames_to_be_averaged_ = stoi(input_word_);
            continue;
        }
        if (input_word_ == "output_precision") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            output_precision_ = stoi(input_word_);
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
        
        std::cerr << "WARNING: no matching input type for: " << input_word_ << " disregarding this variable and continueing to next line\n" << std::endl;
        getline(input_file_, input_word_);
    }

    input_file_.close();
}

void MeanSquareDisplacement::initialize_rest_members()
{
    std::cout << "Initializing...\n" << std::endl;

    // check for legality
    if (calculation_dimension_ != 1 && calculation_dimension_ != 2 && calculation_dimension_ != 3 &&
        calculation_dimension_ != 12 && calculation_dimension_ != 13 && calculation_dimension_ != 23 && calculation_dimension_ != 123) {
        std::cerr << "ERROR: Calculation dimension can only be 1, 2, 3, 12, 13, 23 or 123!" << std::endl;
    }
    if (interval_ < 1.0) {
        std::cerr << "ERROR: Interval should be equal or larger than 1!\n" << std::endl;
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
    
    // construct time table
    compute_time_table();
    if (time_table_.size() < number_of_timepoints_) {
        number_of_timepoints_ = time_table_.size();
        std::cout << "NOTE: Reduce number_of_timepoints to " << number_of_timepoints_ <<" to avoid repetition!" << std::endl;
    }

    number_of_frames_ = time_table_[number_of_timepoints_ - 1] + number_of_frames_to_be_averaged_;
    layer_width_ = (layer_right_point_ - layer_left_point_) / number_of_layers_;
    
    std::cout << "NOTICE: First " << start_frame_ << " frames will be skipped and following " << number_of_frames_ << " frames will be read in.\n" << std::endl;

    std::ifstream trajectory_file_(trajectory_file_name_);

    if (!trajectory_file_) {
        std::cerr << "ERROR: Trajectory file not found!\n" << std::endl;
        exit(1);
    }

    std::string line_stream_;

    for (size_t i_line = 0; i_line < 3; ++i_line) {
        getline(trajectory_file_, line_stream_);
    }
    trajectory_file_ >> number_of_atoms_;
    for (size_t i_line = 0; i_line < dimension_ + 3; ++i_line) {
        getline(trajectory_file_, line_stream_);
    }

    number_of_types_of_atoms_ = 0;
    atom_type_.resize(number_of_atoms_);

    for (size_t i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
        trajectory_file_ >> line_stream_ >> atom_type_[i_atom];
        if (atom_type_[i_atom] > number_of_types_of_atoms_) {
            number_of_types_of_atoms_ = atom_type_[i_atom];
        }
        getline(trajectory_file_, line_stream_);
    }
    trajectory_file_.seekg(0, std::ios::beg);

    atom_layer_.resize(number_of_frames_, std::vector<int>(number_of_atoms_));
    atom_coordinates_.resize(number_of_frames_, std::vector<std::vector<double>>(number_of_atoms_, std::vector<double>(dimension_)));
    mean_square_displacement_.resize(number_of_layers_, std::vector<std::vector<double>>(number_of_timepoints_, std::vector<double>(number_of_types_of_atoms_, 0.0)));
    std::vector<std::vector<std::vector<double>>> box_boundaries_(number_of_frames_, std::vector<std::vector<double>>(dimension_, std::vector<double>(2)));

    // Jump useless frames
    for (size_t i_frame = 0; i_frame < start_frame_; ++i_frame) {
        for (size_t i_line = 0; i_line < number_of_atoms_ + dimension_ + 6; ++i_line) {
            getline(trajectory_file_, line_stream_);
        }
    }

    // Read in trajectories
    for (size_t i_frame = 0; i_frame < number_of_frames_; ++i_frame) {
        std::cout << "\r" << "Reading " << i_frame + 1 << "th frame..." << std::flush;
        
        for (size_t i_line = 0; i_line < 5; ++i_line) {
            getline(trajectory_file_, line_stream_);
        }
        for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
            trajectory_file_ >> box_boundaries_[i_frame][i_dimension][0] >> box_boundaries_[i_frame][i_dimension][1];
            getline(trajectory_file_, line_stream_);
        }
        getline(trajectory_file_, line_stream_);
        for (size_t i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
            trajectory_file_ >> line_stream_ >> line_stream_;
            for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                trajectory_file_ >> atom_coordinates_[i_frame][i_atom][i_dimension];
            }
            getline(trajectory_file_, line_stream_);
        }
    }

    trajectory_file_.close();
    std::cout << "\r" << "Complete reading in input file.\n" << std::endl;

    // Unwrap the trajectory
    if (is_trajectory_wrapped_) {
        for (size_t i_frame = 1; i_frame < number_of_frames_; ++i_frame) {
#pragma omp parallel for
            for (size_t i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
                for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                    atom_coordinates_[i_frame][i_atom][i_dimension] += round(
                            atom_coordinates_[i_frame - 1][i_atom][i_dimension] -
                            atom_coordinates_[i_frame][i_atom][i_dimension]);
                }
            }
        }
    }

    // Unscale the trajectory
#pragma omp parallel for
    for (size_t i_frame = 0; i_frame < number_of_frames_; ++i_frame) {
        for (size_t i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
            for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                double width_box_ = box_boundaries_[i_frame][i_dimension][1] - box_boundaries_[i_frame][i_dimension][0];
                atom_coordinates_[i_frame][i_atom][i_dimension] = atom_coordinates_[i_frame][i_atom][i_dimension] * width_box_ + box_boundaries_[i_frame][i_dimension][0];
            }
        }
    }

    // Compute atom layer
#pragma omp parallel for
    for (size_t i_frame = 0; i_frame < number_of_frames_; ++i_frame) {
        for (size_t i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
            atom_layer_[i_frame][i_atom] = int(floor((fabs(atom_coordinates_[i_frame][i_atom][2]) - layer_left_point_) / layer_width_));
        }
    }
}

void MeanSquareDisplacement::compute_msd()
{
    std::cout << "Computing...\n" << std::endl;

    // Determine the dimension needed to be calculated if only 2 dimensions neeeded
    int dimension_0_  = 0, dimension_1_ = 0, dimension_2_ = 0;
    if (calculation_dimension_ < 10) {
        dimension_0_ = calculation_dimension_ - 1;
    }
    else if (calculation_dimension_ < 100) {
        dimension_1_ = calculation_dimension_ / 10 - 1;
        dimension_2_ = calculation_dimension_ % 10 - 1;
    }

    std::vector< std::vector< std::vector< unsigned int > > > number_of_atoms_for_average_(number_of_layers_, std::vector< std::vector< unsigned int > > (number_of_timepoints_, std::vector< unsigned int > (number_of_types_of_atoms_, 0)));

    size_t status = 1;
#pragma omp parallel for
    for (size_t i_timepoint = 1; i_timepoint < number_of_timepoints_; ++i_timepoint) {
        unsigned int current_time = time_table_[i_timepoint];
        for (size_t i_average = 0; i_average < number_of_frames_to_be_averaged_; ++i_average) {
            for (size_t i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
                if (atom_layer_[i_average][i_atom] == atom_layer_[i_average + current_time][i_atom] && atom_layer_[i_average][i_atom] >= 0 && atom_layer_[i_average][i_atom] < static_cast<int>(number_of_layers_)) {
                    double temp_value_, temp_square_distance_ = 0.0;
                    if (calculation_dimension_ < 10) {
                        temp_value_ = atom_coordinates_[i_average + current_time][i_atom][dimension_0_] -
                                      atom_coordinates_[i_average][i_atom][dimension_0_];
                        temp_square_distance_ += temp_value_ * temp_value_;
                    }
                    else if (calculation_dimension_ < 100) {
                        temp_value_ = atom_coordinates_[i_average + current_time][i_atom][dimension_1_] -
                                      atom_coordinates_[i_average][i_atom][dimension_1_];
                        temp_square_distance_ += temp_value_ * temp_value_;
                        temp_value_ = atom_coordinates_[i_average + current_time][i_atom][dimension_2_] -
                                      atom_coordinates_[i_average][i_atom][dimension_2_];
                        temp_square_distance_ += temp_value_ * temp_value_;
                    }
                    else {
                        for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                            temp_value_ = atom_coordinates_[i_average + current_time][i_atom][i_dimension] -
                                          atom_coordinates_[i_average][i_atom][i_dimension];
                            temp_square_distance_ += temp_value_ * temp_value_;
                        }
                    }
                    mean_square_displacement_[atom_layer_[i_average][i_atom]][i_timepoint][atom_type_[i_atom] - 1] += temp_square_distance_;
                    ++number_of_atoms_for_average_[atom_layer_[i_average][i_atom]][i_timepoint][atom_type_[i_atom] - 1];
                }
            }
        }
#pragma omp critical
{
        print_status(status);
}
    }
    std::cout << "\n" << std::endl;

#pragma omp parallel for
    for (size_t i_layer = 0; i_layer < number_of_layers_; ++i_layer) {
        for (size_t i_timepoint = 1; i_timepoint < number_of_timepoints_; ++i_timepoint) {
            for (size_t i_type = 0; i_type < number_of_types_of_atoms_; ++i_type) {
                mean_square_displacement_[i_layer][i_timepoint][i_type] /= number_of_atoms_for_average_[i_layer][i_timepoint][i_type];
            }
        }
    }
}

void MeanSquareDisplacement::write_output_file()
{
    std::cout << "Writing output file...\n" << std::endl;

    std::ofstream output_file_(output_file_name_);
    output_file_.precision(output_precision_);

    for (size_t i_layer = 0; i_layer < number_of_layers_; ++i_layer) {
        if (print_headline_) {
            output_file_ << "layer = " << layer_left_point_ + (i_layer + 0.5) * layer_width_ << "\n" << "t";
            for (size_t i_type = 1; i_type <= number_of_types_of_atoms_; ++i_type) {
                output_file_ << " atom_" << i_type;
            }
            output_file_ << "\n";
        }
        for (size_t i_timepoint = 0; i_timepoint < number_of_timepoints_; ++i_timepoint){
            output_file_ << time_table_[i_timepoint] * frame_interval_time_;
            for (size_t i_type = 0; i_type < number_of_types_of_atoms_; ++i_type) {
                output_file_ << " " << mean_square_displacement_[i_layer][i_timepoint][i_type];
            }
            output_file_ << "\n";
        }
    }

    output_file_.close();
}

void MeanSquareDisplacement::compute_time_table()
{
    time_table_.push_back(0);
    double total_time = interval_;
    unsigned int previous_frame = 0;
    
    if (time_scale_type_ == "linear") {
        for (size_t i_time = 1; i_time < number_of_timepoints_; ++i_time) {
            time_table_.push_back(static_cast<unsigned int>(round(total_time)));
            total_time += interval_;
        }
    }
    else if (time_scale_type_ == "log") {
        for (size_t i_time = 1; i_time < number_of_timepoints_; ++i_time) {
            unsigned int current_frame = static_cast<unsigned int>(round(total_time));
            if (current_frame != previous_frame) {
                time_table_.push_back(current_frame);
                previous_frame = current_frame;
            }
            total_time *= interval_;
        }
    }
}

void MeanSquareDisplacement::print_status(size_t& status)
{
    std::cout << "\rCurrent progress of calculating the mean squared displacement is: ";
    std::cout << status * 100.0/(number_of_timepoints_ - 1);
    std::cout << " \%";
    std::cout << std::flush;
    ++status;
}
