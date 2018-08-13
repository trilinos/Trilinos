#ifndef STK_STK_BALANCE_STK_BALANCE_INTERNAL_BALANCECOMMANDLINE_HPP_
#define STK_STK_BALANCE_STK_BALANCE_INTERNAL_BALANCECOMMANDLINE_HPP_
#include <mpi.h>
#include <stk_util/command_line/CommandLineParser.hpp>
#include <string>
#include <stk_balance/internal/balanceDefaults.hpp>

namespace stk { class CommandLineParserParallel; }

namespace stk {
namespace balance {

struct CommandLineOptions
{
    stk::CommandLineOption infile{"infile", "i", "undecomposed serial input mesh file"};
    stk::CommandLineOption outputDirectory{"outputDirectory", "o", "output directory for decomposition"};
    stk::CommandLineOption smDefaults{"sm", "", "Use settings suitable for solving Solid Mechanics problems "
                                      "(Graph vertex weights, graph edge weight, and search tolerances "
                                      "appropriate for contact)"};
    stk::CommandLineOption sdDefaults{"sd", "", "Use settings suitable for solving Structural Dynamics"
                                      " problems (handle spider elements)"};
};

struct ParsedOptions
{
    std::string inFile;
    std::string outputDirectory;
    enum AppTypeDefaults appTypeDefaults;
};

std::string get_quick_example(const std::string &execName,
                              const std::string &infileName,
                              const std::string &outputDirectoryName,
                              MPI_Comm comm);

ParsedOptions parse_balance_command_line(int argc,
                                 const char**argv,
                                 const std::string &execName,
                                 MPI_Comm comm);

void print_running_msg(const std::string &execName, const ParsedOptions &balanceOptions, MPI_Comm comm);

}
}

#endif /* STK_STK_BALANCE_STK_BALANCE_INTERNAL_BALANCECOMMANDLINE_HPP_ */

