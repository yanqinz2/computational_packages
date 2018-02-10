#include <fstream>
#include <string>
#include <vector>

class NonGaussianParameter
{
public:
    NonGaussianParameter();
    virtual ~NonGaussianParameter();

    void read_commands(int argc, char* argv[]);
    void read_input_file();
    void initialize_rest_members();
    void compute_ngp();
    void write_output_file();

protected:
    bool print_headline_;

    std::string input_file_name_;
    std::string trajectory_file_name_;
    std::string output_file_name_;

    unsigned int dimension_;
    unsigned int calculation_dimension_;
    unsigned int start_frame_;
    unsigned int interval_;
    unsigned int number_of_timepoints_;
    unsigned int number_of_frames_to_be_averaged_;
    unsigned int number_of_layers_;
    unsigned int output_precision_;
    unsigned int number_of_atoms_;
    unsigned int number_of_types_of_atoms_;

    std::vector< unsigned int > atom_type_;

    // number_of_atoms_for_average_[layer][type][timepoint]
    std::vector< std::vector< std::vector< unsigned int > > > number_of_atoms_for_average_;

    double frame_interval_time_;
    double layer_left_point_;
    double layer_right_point_;
    double layer_width_;

    std::vector< double > left_point_of_box_;
    std::vector< double > right_point_of_box_;
    std::vector< double > box_length_;
    std::vector< double > time_table_;

    // coordinate_[atom][dimension]
    std::vector< std::vector< double > > base_atom_coordinate_;
    std::vector< std::vector< double > > measurement_atom_coordinate_;

    // mean_displacement_[layer][type][timepoint]
    std::vector< std::vector< std::vector< double > > > mean_square_displacement_;
    std::vector< std::vector< std::vector< double > > > mean_quartic_displacement_;
    std::vector< std::vector< std::vector< double > > > non_gaussian_parameter_;

    std::ifstream* file_pointer_;

private:
};