#include <string>
#include <vector>

class UnwrapTrajectory
{
public:
    UnwrapTrajectory();
    virtual ~UnwrapTrajectory();

    void read_commands(int argc, char* argv[]);
    void initialize_rest_members();
    void unwrap_trajectory();

protected:
    std::string input_file_name_;
    std::string trajectory_file_name_;
    std::string output_file_name_;

    unsigned int start_frame_;
    unsigned int end_frame_;
    unsigned int dimension_;
    unsigned int output_precision_;
    unsigned int number_of_atoms_;

    std::vector< unsigned int > atom_id_;
    std::vector< unsigned int > atom_type_;

    std::vector< std::vector < double > > trajectory_front_;
    std::vector< std::vector < double > > trajectory_back_;

private:
};