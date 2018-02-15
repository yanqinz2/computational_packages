#include <string>
#include <vector>

class VelocityAutoCorrelationFunction
{
public:
    VelocityAutoCorrelationFunction();
    virtual ~VelocityAutoCorrelationFunction();

    void read_commands(int argc, char* argv[]);
    void read_input_file();
    void initialize_rest_members();
    void compute_vacf();
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
    std::vector< int > base_atom_layer_;
    std::vector< int > measurement_atom_layer_;

    // number_of_atoms_for_average_[layer][type][timepoint]
    std::vector< std::vector< std::vector< unsigned int > > > number_of_atoms_for_average_;

    double frame_interval_time_;
    double layer_left_point_;
    double layer_right_point_;
    double layer_width_;

    std::vector< double > time_table_;
    std::vector< double > base_atom_z_coordinate_;
    std::vector< double > measurement_atom_z_coordinate_;

    // coordinate_[atom][dimension] mean_square_displacement_[layer][type][timepoint]
    std::vector< std::vector< double > > base_atom_velocity_;
    std::vector< std::vector< double > > measurement_atom_velocity_;
    std::vector< std::vector< std::vector< double > > > velocity_auto_correlation_function_;

    std::FILE** file_pointer_;

private:
};