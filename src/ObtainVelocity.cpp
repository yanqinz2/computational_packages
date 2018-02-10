#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifndef OBTAINVELOCITY_HPP
#define OBTAINVELOCITY_HPP
#include "ObtainVelocity.hpp"
#endif

#ifndef OMP
#define OMP
#include "omp.h"
#endif

ObtainVelocity::ObtainVelocity():
        trajectory_file_name_(""),
        output_file_name_(""),
        start_frame_(1),
        end_frame_(1),
        dimension_(3),
        output_precision_(15),
        frame_interval_time_(0.005)
{
}

ObtainVelocity::~ObtainVelocity()
{
}

void ObtainVelocity::read_commands(int argc, char* argv[])
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

void ObtainVelocity::initialize_rest_members()
{
    std::cout << "Initializing...\n" << std::endl;

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
    }
    input_file_.close();

    // check for legality
    if (start_frame_ == 0) {
        std::cerr << "ERROR: Start frame should be positive!\n" << std::endl;
        exit(1);
    }
    if (end_frame_ < start_frame_) {
        std::cerr << "ERROR: End frame should be greater than start frame!\n" << std::endl;
        exit(1);
    }
    if (frame_interval_time_ < 0.0) {
        std::cerr << "ERROR: Frame interval time should be positive!\n" << std::endl;
        exit(1);
    }

    std::ifstream trajectory_file_(trajectory_file_name_);

    if (!trajectory_file_) {
        std::cerr << "ERROR: Trajectory file not found!\n" << std::endl;
        exit(1);
    }

    std::string line_stream_;
    unsigned int atom_type_temp_;

    for (int i_line = 0; i_line < 3; ++i_line) {
        getline(trajectory_file_, line_stream_);
    }
    trajectory_file_ >> number_of_atoms_;

    trajectory_file_.close();

    atom_id_.resize(number_of_atoms_);

    atom_type_.resize(number_of_atoms_);

    left_point_of_the_box_.resize(dimension_);

    right_point_of_the_box_.resize(dimension_);

    box_length_.resize(dimension_);

    trajectory_front_.resize(number_of_atoms_);
    for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
        trajectory_front_[i_atom].resize(dimension_);
    }

    trajectory_back_.resize(number_of_atoms_);
    for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
        trajectory_back_[i_atom].resize(dimension_);
    }

    velocity_.resize(number_of_atoms_);
    for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
        velocity_[i_atom].resize(dimension_);
    }
}

void ObtainVelocity::compute_and_write_output_file()
{
    std::cout << "Computing...\n" << std::endl;

    unsigned int number_of_frames_ = end_frame_ - start_frame_ + 1;
    std::ifstream trajectory_file_(trajectory_file_name_);
    std::ofstream output_file_(output_file_name_);
    std::string line_stream_;

    output_file_.precision(output_precision_);

    // jump over needless frames
    for (int i_frame = 1; i_frame < start_frame_; ++i_frame) {
        for (int i_line = 0; i_line < number_of_atoms_ + dimension_ + 6; ++i_line) {
            getline(trajectory_file_, line_stream_);
        }
    }

    // read in first frame
    for (int i_line = 0; i_line < 5; ++i_line) {
        getline(trajectory_file_, line_stream_);
    }
    for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
        trajectory_file_ >> left_point_of_the_box_[i_dimension] >> right_point_of_the_box_[i_dimension];
        box_length_[i_dimension] = right_point_of_the_box_[i_dimension] - left_point_of_the_box_[i_dimension];
        getline(trajectory_file_, line_stream_);
    }
    getline(trajectory_file_, line_stream_);
    for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
        trajectory_file_ >> atom_id_[i_atom] >> atom_type_[i_atom];
        for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
            trajectory_file_ >> trajectory_back_[i_atom][i_dimension];
        }
        getline(trajectory_file_, line_stream_);
    }
#pragma omp parallel for
    for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
        for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
            trajectory_back_[i_atom][i_dimension] = trajectory_back_[i_atom][i_dimension]
                                                    * box_length_[i_dimension] + left_point_of_the_box_[i_dimension];
        }
    }

    // start compute velocity and output
    for (int i_frame = 1; i_frame < number_of_frames_; ++i_frame) {
        // read and output fisrt 9 lines
        for (int i_line = 0; i_line < 5; ++i_line) {
            getline(trajectory_file_, line_stream_);
            output_file_ << line_stream_ << "\n";
        }
        for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
            trajectory_file_ >> left_point_of_the_box_[i_dimension] >> right_point_of_the_box_[i_dimension];
            box_length_[i_dimension] = right_point_of_the_box_[i_dimension] - left_point_of_the_box_[i_dimension];
            getline(trajectory_file_, line_stream_);
            output_file_ << left_point_of_the_box_[i_dimension] << " " << right_point_of_the_box_[i_dimension] << "\n";
        }
        getline(trajectory_file_, line_stream_);
        output_file_ << line_stream_ << "\n";

        // read trajectory into front vessel
        for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
            trajectory_file_ >> line_stream_ >> line_stream_;
            for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                trajectory_file_ >> trajectory_front_[i_atom][i_dimension];
            }
            getline(trajectory_file_, line_stream_);
        }
#pragma omp parallel for
        for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
            for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                trajectory_front_[i_atom][i_dimension] = trajectory_front_[i_atom][i_dimension]
                                                         * box_length_[i_dimension] + left_point_of_the_box_[i_dimension];
            }
        }

        // compute the velocity
#pragma omp parallel for
        for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
            for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                velocity_[i_atom][i_dimension] = trajectory_front_[i_atom][i_dimension] - trajectory_back_[i_atom][i_dimension];
                if (velocity_[i_atom][i_dimension] > 0.5 * box_length_[i_dimension]) {
                    velocity_[i_atom][i_dimension]  -= box_length_[i_dimension];
                }
                else if (velocity_[i_atom][i_dimension] < - 0.5 * box_length_[i_dimension]) {
                    velocity_[i_atom][i_dimension] += box_length_[i_dimension];
                }
                velocity_[i_atom][i_dimension] /= frame_interval_time_;
            }
        }

        // output the velocity
        for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
            output_file_ << atom_id_[i_atom] << " " << atom_type_[i_atom];
            for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                output_file_ << " " << velocity_[i_atom][i_dimension];
            }
            output_file_ << " " << trajectory_front_[i_atom][2] << "\n";
        }

        trajectory_back_.swap(trajectory_front_);

        std::cout << "\r" << float(i_frame + 1) / number_of_frames_ * 100 << "% has been completed..." << std::flush;
    }

    std::cout << "\n" << std::endl;

    trajectory_file_.close();
    output_file_.close();
}