#include <string>
#include <vector>

class ObtainVelocity
{
public:
    ObtainVelocity();
    virtual ~ObtainVelocity();

    void read_commands(int argc, char* argv[]);
    void initialize_rest_members();
    void compute_and_write_output_file();

protected:
    std::string input_file_name_;
    std::string trajectory_file_name_;
    std::string output_file_name_;

    unsigned int start_frame_;
    unsigned int end_frame_;
    unsigned int dimension_;
    unsigned int output_precision_;
    unsigned int number_of_atoms_;

    double frame_interval_time_;

    std::vector< unsigned int > atom_id_;
    std::vector< unsigned int > atom_type_;
    std::vector< double > left_point_of_the_box_;
    std::vector< double > right_point_of_the_box_;
    std::vector< double > box_length_;

    std::vector< std::vector < double > > trajectory_front_;
    std::vector< std::vector < double > > trajectory_back_;
    std::vector< std::vector < double > > velocity_;

private:
};