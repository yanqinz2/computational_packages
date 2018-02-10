#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifndef UNWRAPTRAJECTORY_HPP
#define UNWRAPTRAJECTORY_HPP
#include "UnwrapTrajectory.hpp"
#endif

#ifndef OMP
#define OMP
#include "omp.h"
#endif

UnwrapTrajectory::UnwrapTrajectory():
    trajectory_file_name_(""),
    output_file_name_(""),
    start_frame_(1),
    end_frame_(1),
    dimension_(3),
    output_precision_(15)
{
}

UnwrapTrajectory::~UnwrapTrajectory()
{
}

void UnwrapTrajectory::read_commands(int argc, char* argv[])
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

void UnwrapTrajectory::initialize_rest_members()
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
    }
    input_file_.close();

    // check for legality
    if (start_frame_ == 0) {
        std::cerr << "ERROR: Start frame should be equal or greater than 1!\n" << std::endl;
        exit(1);
    }
    if (end_frame_ < start_frame_) {
        std::cerr << "ERROR: End frame should be equal or greater than start frame!\n" << std::endl;
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

    trajectory_front_.resize(number_of_atoms_);
    for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
        trajectory_front_[i_atom].resize(dimension_);
    }

    trajectory_back_.resize(number_of_atoms_);
    for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
        trajectory_back_[i_atom].resize(dimension_);
    }
}

void UnwrapTrajectory::unwrap_trajectory()
{
    std::cout << "Unwrapping...\n" << std::endl;

    unsigned int number_of_frames_ = end_frame_ - start_frame_ + 1;

    std::ifstream trajectory_file_(trajectory_file_name_);

    std::ofstream output_file_(output_file_name_);
    output_file_.precision(output_precision_);

    std::string line_stream_;

    // jump over needless frames
    for (int i_frame = 1; i_frame < start_frame_; ++i_frame) {
        for (int i_line = 0; i_line < number_of_atoms_ + dimension_ + 6; ++i_line) {
            getline(trajectory_file_, line_stream_);
        }
    }

    // read and output the first frame
    for (int i_line = 0; i_line < dimension_ + 6; ++i_line) {
        getline(trajectory_file_, line_stream_);
        output_file_ << line_stream_ << "\n";
    }
    for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
        trajectory_file_ >> atom_id_[i_atom] >> atom_type_[i_atom];
        for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
            trajectory_file_ >> trajectory_back_[i_atom][i_dimension];
        }
        getline(trajectory_file_, line_stream_);
        output_file_ << atom_id_[i_atom] << " " << atom_type_[i_atom];
        for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
            output_file_ << " " << trajectory_back_[i_atom][i_dimension];
        }
        output_file_ << "\n";
    }

    // start unwrap
    for (int i_frame = 1; i_frame < number_of_frames_; ++i_frame) {
        // read and output fisrt 9 lines
        for (int i_line = 0; i_line < dimension_ + 6; ++i_line) {
            getline(trajectory_file_, line_stream_);
            output_file_ << line_stream_ << "\n";
        }

        // read trajectory into front vessel
        for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
            trajectory_file_ >> line_stream_ >> line_stream_;
            for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                trajectory_file_ >> trajectory_front_[i_atom][i_dimension];
            }
            getline(trajectory_file_, line_stream_);
        }

        // unwrap the front trajectory
#pragma omp parallel for
        for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
            for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                trajectory_front_[i_atom][i_dimension] +=
                        round(trajectory_back_[i_atom][i_dimension] - trajectory_front_[i_atom][i_dimension]);
            }
        }

        // output the unwrapped trajectory
        for (int i_atom = 0; i_atom < number_of_atoms_; ++i_atom) {
            output_file_ << atom_id_[i_atom] << " " << atom_type_[i_atom];
            for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                output_file_ << " " << trajectory_front_[i_atom][i_dimension];
            }
            output_file_ << "\n";
        }

        trajectory_back_.swap(trajectory_front_);

        std::cout << "\r" << float(i_frame + 1) / number_of_frames_ * 100 << "% has been completed..." << std::flush;
    }

    std::cout << "\n" << std::endl;

    trajectory_file_.close();
    output_file_.close();
}