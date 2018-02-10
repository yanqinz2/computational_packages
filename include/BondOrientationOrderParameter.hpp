#include <string>
#include <vector>

class BondOrientationOrderParameter
{
public:
    BondOrientationOrderParameter();
    virtual ~BondOrientationOrderParameter();

    void read_commands(int argc, char* argv[]);
    void read_input_file();
    void initialize_rest_members();
    void compute_bop();
    void write_output_file();

protected:
    bool print_headline_;

    std::string input_file_name_;
    std::string trajectory_file_name_;
    std::string output_file_name_;

    unsigned int start_frame_;
    unsigned int interval_;
    unsigned int number_of_frames_;
    unsigned int number_of_atoms_;
    unsigned int dimension_;
    unsigned int number_of_types_of_atoms_;
    unsigned int number_of_bins_;
    unsigned int output_precision_;
    unsigned int object_atom_type_;
    unsigned int measurement_atom_type_;
    unsigned int number_object_atom_;
    unsigned int number_measurement_atom_;
    unsigned int bond_order_parameter_;

    std::vector< unsigned int > number_of_atoms_of_each_type_;
    std::vector< unsigned int > atom_type_;

    double layer_left_point_;
    double layer_right_point_;
    double cutoff_distance_;

    std::vector< double > left_point_of_box_;
    std::vector< double > right_point_of_box_;
    std::vector< double > box_length_;
    std::vector< double > bond_orientation_parameter_;
    std::vector< double > probability_of_bop_;

    // coordinate_[atom][dimension] theta_/phi_[object_atom][measurement_atom]
    std::vector< std::vector< double > > object_atom_coordinate_;
    std::vector< std::vector< double > > measurement_atom_coordinate_;
    std::vector< std::vector< double > > real_term_;
    std::vector< std::vector< double > > imaginary_term_;
    std::vector< std::vector< double > > orientated_distance_;
};