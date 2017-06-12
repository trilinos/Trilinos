#ifndef STK_STK_UTIL_STK_UTIL_COMMAND_LINE_COMMANDLINEPARSERUTILS_HPP_
#define STK_STK_UTIL_STK_UTIL_COMMAND_LINE_COMMANDLINEPARSERUTILS_HPP_
#include <mpi.h>
#include <string>

namespace stk {

class CommandLineParserParallel;

std::string angle_it(const std::string &s);
std::string bracket_it(const std::string &s);
std::string dash_it(const std::string &s);
std::string get_quick_error(const std::string &execName, const std::string &quickExample);
std::string get_version(const std::string &executableName);
void parse_command_line(int argc,
                        const char** argv,
                        const std::string& quickExample,
                        const std::string& longExample,
                        stk::CommandLineParserParallel& commandLine,
                        MPI_Comm comm);
namespace parallel {
void print_and_exit(const std::string &msg, MPI_Comm comm);
void require(bool requirement, const std::string &msg, MPI_Comm comm);
bool does_file_exist(const std::string& filename);
void require_file_exists(const std::string& inFile, const std::string& execName, const std::string& quickExample, MPI_Comm comm);
}


}

#endif /* STK_STK_UTIL_STK_UTIL_COMMAND_LINE_COMMANDLINEPARSERUTILS_HPP_ */
