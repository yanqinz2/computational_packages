#include <fstream>
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
    unsigned int number_of_frames_;

    std::vector< unsigned int > atom_type_;

    // number_of_atoms_for_average_[layer][type][timepoint]
    std::vector< std::vector< std::vector< unsigned int > > > number_of_atoms_for_average_;

    std::vector< std::vector< int > >atom_layer_;

    double frame_interval_time_;
    double layer_left_point_;
    double layer_right_point_;
    double layer_width_;

    std::vector< double > time_table_;

    // atom_velocities_[frame][atom][dimension] mean_square_displacement_[layer][type][timepoint]
    std::vector< std::vector< std::vector< double > > > atom_velocities_;
    std::vector< std::vector< std::vector< double > > > velocity_auto_correlation_function_;

private:
};