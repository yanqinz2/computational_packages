#include <string>
#include <vector>

class AtomNumberDensity
{
public:
    AtomNumberDensity();
    virtual ~AtomNumberDensity();

    void read_commands(int argc, char* argv[]);
    void read_input_file();
    void initialize_rest_members();
    void compute_atom_number_density();
    void write_output_file();

protected:
    // protected members
    bool print_headline_;

    std::string input_file_name_;
    std::string trajectory_file_name_;
    std::string output_file_name_;

    unsigned int start_frame_;
    unsigned int end_frame_;
    unsigned int number_of_frames_;
    unsigned int number_of_atoms_;
    unsigned int dimension_;
    unsigned int number_of_types_of_atoms_;
    unsigned int number_of_layers_;
    unsigned int output_precision_;

    std::vector< unsigned int > number_of_atoms_of_each_type_;
    std::vector< unsigned int > atom_type_;

    double width_of_layers_;

    std::vector< double > left_point_of_box_;
    std::vector< double > right_point_of_box_;
    std::vector< double > box_length_;

    // trajectory_[i_atom][i_dimension] number_density_of_each_layer_[i_layer][i_type]
    std::vector< std::vector< double > > trajectory_;
    std::vector< std::vector< double > > number_density_of_each_layer_;

private:
};
