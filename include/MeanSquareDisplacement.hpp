#include <string>
#include <vector>

class MeanSquareDisplacement
{
public:
    MeanSquareDisplacement();
    virtual ~MeanSquareDisplacement();

    void read_commands(int argc, char* argv[]);
    void read_input_file();
    void initialize_rest_members();
    void compute_msd();
    void write_output_file();

protected:
    std::string input_file_name_;
    std::string trajectory_file_name_;
    std::string output_file_name_;
    std::string time_scale_type_;
    
    bool print_headline_;
    bool is_trajectory_wrapped_;

    unsigned int dimension_;
    unsigned int calculation_dimension_;
    unsigned int start_frame_;
    unsigned int number_of_timepoints_;
    unsigned int number_of_frames_to_be_averaged_;
    unsigned int number_of_layers_;
    unsigned int output_precision_;
    unsigned int number_of_atoms_;
    unsigned int number_of_types_of_atoms_;
    unsigned int number_of_frames_;

    double interval_;
    double frame_interval_time_;
    double layer_left_point_;
    double layer_right_point_;
    double layer_width_;

    std::vector<unsigned int> time_table_;
    std::vector<unsigned int> atom_type_;                                                  // [atom]
    std::vector<std::vector<int>> atom_layer_;                                           // [frame][atom]
    std::vector<std::vector<std::vector<double>>> atom_coordinates_;                  // [frame][atom][dimension]
    std::vector<std::vector<std::vector<double>>>  mean_square_displacement_;         // [layer][timepoint][type]

private:
    void compute_time_table();
    void print_status(size_t&);
};
