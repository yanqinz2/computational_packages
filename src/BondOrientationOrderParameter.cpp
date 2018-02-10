#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifndef BONDORIENTATIONORDERPRAMETER_HPP
#define BONDORIENTATIONORDERPRAMETER_HPP
#include "BondOrientationOrderParameter.hpp"
#endif

#ifndef SPHERICAL_HARMONIC_HPP
#define SPHERICAL_HARMONIC_HPP
#include <boost/math/special_functions/spherical_harmonic.hpp>
#endif

#ifndef OMP
#define OMP
#include "omp.h"
#endif

BondOrientationOrderParameter::BondOrientationOrderParameter():
        print_headline_(false),
        trajectory_file_name_(""),
        output_file_name_(""),
        start_frame_(1),
        interval_(1),
        number_of_frames_(1),
        dimension_(3),
        number_of_bins_(50),
        output_precision_(15),
        object_atom_type_(1),
        measurement_atom_type_(1),
        bond_order_parameter_(4),
        layer_left_point_(0.0),
        layer_right_point_(0.0),
        cutoff_distance_(0.0)
{
}

BondOrientationOrderParameter::~BondOrientationOrderParameter()
{
}

void BondOrientationOrderParameter::read_commands(int argc, char* argv[])
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

void BondOrientationOrderParameter::read_input_file()
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
        if (input_word_ == "interval") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            interval_ = stod(input_word_);
            continue;
        }
        if (input_word_ == "number_of_frames") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            number_of_frames_ = stod(input_word_);
            continue;
        }
        if (input_word_ == "number_of_bins") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            number_of_bins_ = stod(input_word_);
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
        if (input_word_ == "object_atom_type") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            object_atom_type_ = stod(input_word_);
            continue;
        }
        if (input_word_ == "measurement_atom_type") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            measurement_atom_type_ = stod(input_word_);
            continue;
        }
        if (input_word_ == "bond_order_parameter") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            bond_order_parameter_ = stod(input_word_);
            continue;
        }
        if (input_word_ == "left_point_of_layer") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            layer_left_point_ = stod(input_word_);
            continue;
        }
        if (input_word_ == "right_point_of_layer") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            layer_right_point_ = stod(input_word_);
            continue;
        }
        if (input_word_ == "cutoff_distance") {
            input_file_ >> input_word_;
            if (input_word_[0] == '=') {
                input_file_ >> input_word_;
            }
            cutoff_distance_ = stod(input_word_);
            continue;
        }
    }

    input_file_.close();
}

void BondOrientationOrderParameter::initialize_rest_members()
{
    std::cout << "Initializing...\n" << std::endl;

    if (start_frame_ == 0) {
        std::cerr << "ERROR: Start frame should be positive!\n" << std::endl;
        exit(1);
    }
    if (interval_ == 0) {
        std::cerr << "ERROR: Interval should be positive!\n" << std::endl;
        exit(1);
    }
    if (number_of_frames_ == 0) {
        std::cerr << "ERROR: Number of frames should be positive!\n" << std::endl;
        exit(1);
    }
    if (number_of_bins_ == 0) {
        std::cerr << "ERROR: Number of bins should be positive!\n" << std::endl;
        exit(1);
    }
    if (layer_left_point_ < 0) {
        std::cerr << "ERROR: Left point of the layer should be zero or positive!\n" << std::endl;
        exit(1);
    }
    if (layer_right_point_ < layer_left_point_) {
        std::cerr << "ERROR: Right point of the layer should be larger than the left!\n" << std::endl;
        exit(1);
    }

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
        getline(trajectory_file_, line_stream_);
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

    if (object_atom_type_) {
        number_object_atom_ = number_of_atoms_of_each_type_[object_atom_type_ - 1];
    }
    else {
        number_object_atom_ = number_of_atoms_;
    }

    if (measurement_atom_type_) {
        number_measurement_atom_ = number_of_atoms_of_each_type_[measurement_atom_type_ - 1];
    }
    else {
        number_measurement_atom_ = number_of_atoms_;
    }

    atom_type_.resize(number_of_atoms_);

    bond_orientation_parameter_.resize(number_object_atom_);

    probability_of_bop_.resize(number_of_bins_);

    object_atom_coordinate_.resize(number_object_atom_);
    for (int i_atom = 0; i_atom < number_object_atom_; ++i_atom) {
        object_atom_coordinate_[i_atom].resize(dimension_);
    }

    measurement_atom_coordinate_.resize(number_measurement_atom_);
    for (int i_atom =0; i_atom < number_measurement_atom_; ++i_atom) {
        measurement_atom_coordinate_[i_atom].resize(dimension_);
    }

    real_term_.resize(number_object_atom_);
    for (int i_object = 0; i_object < number_object_atom_; ++i_object) {
        real_term_[i_object].resize(2 * bond_order_parameter_ + 1);
    }

    imaginary_term_.resize(number_object_atom_);
    for (int i_object = 0; i_object < number_object_atom_; ++i_object) {
        imaginary_term_[i_object].resize(2 * bond_order_parameter_ + 1);
    }

    orientated_distance_.resize(number_object_atom_);
    for (int i_object = 0; i_object < number_object_atom_; ++i_object) {
        orientated_distance_[i_object].resize(dimension_);
    }
}

void BondOrientationOrderParameter::compute_bop()
{
    std::cout << "Computing...\n" << std::endl;

    std::ifstream trajectory_file_(trajectory_file_name_);
    std::string line_stream_;

    for (int i_bin = 0; i_bin < number_of_bins_; ++i_bin) {
        probability_of_bop_[i_bin] = 0.0;
    }

    // jump over needless frames
    for (int i_frame = 1; i_frame < start_frame_; ++i_frame) {
        for (int i_line = 0; i_line < number_of_atoms_ + dimension_ + 6; ++i_line) {
            getline(trajectory_file_, line_stream_);
        }
    }

    unsigned int statistic_atoms_ = 0;

    for (int i_frame = 0; i_frame < number_of_frames_; ++i_frame) {
        unsigned int temp_atom_type_;

        // read in necessary trajectory and jump over interval frames
        for (int i_line = 0; i_line < dimension_ + 6; ++i_line) {
            getline(trajectory_file_, line_stream_);
        }

        for (int i_atom = 0, i_object = 0, i_measurement = 0; i_atom < number_of_atoms_; ++i_atom) {
            trajectory_file_ >> line_stream_ >> temp_atom_type_;

            if (object_atom_type_) {
                if (measurement_atom_type_) {
                    if (object_atom_type_ == measurement_atom_type_) {
                        if (temp_atom_type_ == object_atom_type_) {
                            for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                                trajectory_file_ >> object_atom_coordinate_[i_object][i_dimension];
                                object_atom_coordinate_[i_object][i_dimension] =
                                        object_atom_coordinate_[i_object][i_dimension] * box_length_[i_dimension]
                                        + left_point_of_box_[i_dimension];
                                measurement_atom_coordinate_[i_measurement][i_dimension] = object_atom_coordinate_[i_object][i_dimension];
                            }
                            ++i_object;
                            ++i_measurement;
                        }
                        else {
                            if (temp_atom_type_ == object_atom_type_) {
                                for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                                    trajectory_file_ >> object_atom_coordinate_[i_object][i_dimension];
                                    object_atom_coordinate_[i_object][i_dimension] =
                                            object_atom_coordinate_[i_object][i_dimension] * box_length_[i_dimension]
                                            + left_point_of_box_[i_dimension];
                                }
                                ++i_object;
                            }
                            if (temp_atom_type_ == measurement_atom_type_) {
                                for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                                    trajectory_file_ >> measurement_atom_coordinate_[i_measurement][i_dimension];
                                    measurement_atom_coordinate_[i_measurement][i_dimension] =
                                            measurement_atom_coordinate_[i_measurement][i_dimension] * box_length_[i_dimension]
                                            + left_point_of_box_[i_dimension];
                                }
                                ++i_measurement;
                            }
                        }
                    }
                }
                else {
                    for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                        trajectory_file_ >> measurement_atom_coordinate_[i_measurement][i_dimension];
                        measurement_atom_coordinate_[i_measurement][i_dimension] =
                                measurement_atom_coordinate_[i_measurement][i_dimension] * box_length_[i_dimension]
                                + left_point_of_box_[i_dimension];
                    }
                    if (temp_atom_type_ == object_atom_type_) {
                        for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                            object_atom_coordinate_[i_object][i_dimension] =
                                    measurement_atom_coordinate_[i_measurement][i_dimension];
                        }
                        ++i_object;
                    }
                    ++i_measurement;
                }
            }
            else {
                if (measurement_atom_type_) {
                    for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                        trajectory_file_ >> object_atom_coordinate_[i_object][i_dimension];
                        object_atom_coordinate_[i_object][i_dimension] =
                                object_atom_coordinate_[i_object][i_dimension] * box_length_[i_dimension]
                                + left_point_of_box_[i_dimension];
                    }
                    if (temp_atom_type_ == measurement_atom_type_) {
                        for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                            measurement_atom_coordinate_[i_measurement][i_dimension] =
                                    object_atom_coordinate_[i_object][i_dimension];
                        }
                        ++i_measurement;
                    }
                    ++i_object;
                }
                else {
                    for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                        trajectory_file_ >> object_atom_coordinate_[i_object][i_dimension];
                        object_atom_coordinate_[i_object][i_dimension] =
                                object_atom_coordinate_[i_object][i_dimension] * box_length_[i_dimension]
                                + left_point_of_box_[i_dimension];
                        measurement_atom_coordinate_[i_measurement][i_dimension] = object_atom_coordinate_[i_object][i_dimension];
                    }
                    ++i_object;
                    ++i_measurement;
                }
            }

            getline(trajectory_file_, line_stream_);
        }

        for (int i_line = 0; i_line < (interval_ - 1) * (number_of_atoms_ + dimension_ + 6); ++i_line) {
            getline(trajectory_file_, line_stream_);
        }

#pragma omp parallel
        {
            // compute bond orientation order parameter
#pragma omp for
            for (int i_atom = 0; i_atom < number_object_atom_; ++i_atom) {
                bond_orientation_parameter_[i_atom] = 0.0;
            }
#pragma omp for
            for (int i_object = 0; i_object < number_object_atom_; ++i_object) {
                unsigned int neighbor_atoms_ = 0;
                double atom_distance_, theta_, phi_;
                if (fabs(object_atom_coordinate_[i_object][2]) < layer_left_point_
                    || fabs(object_atom_coordinate_[i_object][2]) > layer_right_point_) {
                    // mark for those object atoms outside of the statistical region
                    bond_orientation_parameter_[i_object] = - 1.0;
                    continue;
                }
                for (int m = 0; m < 2 * bond_order_parameter_ + 1; ++m) {
                    real_term_[i_object][m] = imaginary_term_[i_object][m] = 0.0;
                }
                for (int i_measurement = 0; i_measurement < number_measurement_atom_; ++i_measurement) {
                    atom_distance_ = 0.0;

                    // compute orientated distance and distance
                    for (int i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                        if (measurement_atom_coordinate_[i_measurement][i_dimension]
                            - object_atom_coordinate_[i_object][i_dimension] > 0.5 * box_length_[i_dimension]) {
                            orientated_distance_[i_object][i_dimension] = measurement_atom_coordinate_[i_measurement][i_dimension]
                                                                - object_atom_coordinate_[i_object][i_dimension] -
                                                                box_length_[i_dimension];
                        } else if (measurement_atom_coordinate_[i_measurement][i_dimension]
                                   - object_atom_coordinate_[i_object][i_dimension] < - 0.5 * box_length_[i_dimension]) {
                            orientated_distance_[i_object][i_dimension] = measurement_atom_coordinate_[i_measurement][i_dimension]
                                                                - object_atom_coordinate_[i_object][i_dimension] +
                                                                box_length_[i_dimension];
                        } else {
                            orientated_distance_[i_object][i_dimension] = measurement_atom_coordinate_[i_measurement][i_dimension]
                                                                - object_atom_coordinate_[i_object][i_dimension];
                        }
                        if (fabs(orientated_distance_[i_object][i_dimension]) > cutoff_distance_) {
                            atom_distance_ = 0.0;
                            break;
                        }
                        atom_distance_ += orientated_distance_[i_object][i_dimension] * orientated_distance_[i_object][i_dimension];
                    }

                    if (atom_distance_ < cutoff_distance_ * cutoff_distance_ && atom_distance_ > 0.1) {
                        atom_distance_ = sqrt(atom_distance_);
                        theta_ = acos(orientated_distance_[i_object][2] / atom_distance_);
                        phi_ = atan(orientated_distance_[i_object][1] / orientated_distance_[i_object][0]);

                        for (int m = 0; m <= 2 * bond_order_parameter_; ++m) {
                            real_term_[i_object][m] += boost::math::spherical_harmonic_r(bond_order_parameter_,
                                                                               m - bond_order_parameter_, theta_, phi_);
                            imaginary_term_[i_object][m] += boost::math::spherical_harmonic_i(bond_order_parameter_,
                                                                                    m - bond_order_parameter_, theta_,
                                                                                    phi_);
                        }
                        ++neighbor_atoms_;
                    }
                }
                if (neighbor_atoms_ != 0 && neighbor_atoms_ != 1) {
                    for (int m = 0; m < 2 * bond_order_parameter_ + 1; ++m) {
                        bond_orientation_parameter_[i_object] +=
                                (real_term_[i_object][m] * real_term_[i_object][m] + imaginary_term_[i_object][m] * imaginary_term_[i_object][m]);
                    }
                    bond_orientation_parameter_[i_object] /= neighbor_atoms_ * neighbor_atoms_;
                    bond_orientation_parameter_[i_object] =
                            sqrt(4 * M_PI * bond_orientation_parameter_[i_object] / (2 * bond_order_parameter_ + 1));
                } else {
                    // mark for those with no neighbor atoms or have only one neighbor
                    bond_orientation_parameter_[i_object] = -1.0;
                }
            }
        }

        for (int i_object = 0; i_object < number_object_atom_; ++i_object) {
            if (bond_orientation_parameter_[i_object] > 0.0) {
                if (bond_orientation_parameter_[i_object] >= 1.0) {
                    ++probability_of_bop_[number_of_bins_ - 1];
                    ++statistic_atoms_;
                } else {
                    ++probability_of_bop_[int(bond_orientation_parameter_[i_object] * number_of_bins_)];
                    ++statistic_atoms_;
                }
            }
        }

        std::cout << "\r" << float(i_frame + 1) / number_of_frames_ * 100 << "% has been completed..." << std::flush;
    }
    
    trajectory_file_.close();
    
    std::cout << "\n" << std::endl;

#pragma omp parallel for
    for (int i_bin = 0; i_bin < number_of_bins_; ++i_bin) {
        probability_of_bop_[i_bin] /= statistic_atoms_;
    }
}

void BondOrientationOrderParameter::write_output_file()
{
    std::cout << "Writing output file...\n" << std::endl;

    std::ofstream output_file_(output_file_name_);
    output_file_.precision(output_precision_);

    if (print_headline_) {
        output_file_ << "BOP     Probability\n";
    }

    for (int i_bin = 0; i_bin < number_of_bins_; ++i_bin) {
        output_file_ << (i_bin + 0.5) / number_of_bins_ << " " << probability_of_bop_[i_bin] << "\n";
    }

    output_file_.close();
}