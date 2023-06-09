#include "stk_util/command_line/CommandLineParserUtils.hpp"
#include "stk_util/command_line/CommandLineParser.hpp"          // for CommandLineParser, Comman...
#include "stk_util/command_line/CommandLineParserParallel.hpp"  // for CommandLineParserParallel
#include "stk_util/environment/Env.hpp"                         // for outputP0
#include "stk_util/parallel/Parallel.hpp"                       // for MPI_Finalize, parallel_ma...
#include "stk_util/registry/ProductRegistry.hpp"                // for get_version
#include "stk_util/util/ReportHandler.hpp"                      // for ThrowRequireMsg
#include "stk_util/util/string_utils.hpp"                       // for tailname
#include <cstdlib>                                              // for exit
#include <fstream>                                              // for endl, basic_ostream, ostream
#include <string>                                               // for allocator, string, operator+

namespace stk {

std::string get_quick_error(const std::string &execName, const std::string &quickExample)
{
    std::string s = quickExample;
    s += "Use '" + execName + " --help' for more information.";
    return s;
}

void parse_command_line(int argc,
                        const char** argv,
                        const std::string& quickExample,
                        const std::string& longExample,
                        stk::CommandLineParserParallel& commandLine,
                        MPI_Comm comm)
{
    std::string execName = stk::tailname(argv[0]);
    stk::CommandLineParser::ParseState state = commandLine.parse(argc, argv);
    if(state == stk::CommandLineParser::ParseVersionOnly)
        stk::parallel::print_and_exit(stk::get_version(execName), comm);

    if(state == stk::CommandLineParser::ParseHelpOnly)
    {
        std::string usage = quickExample + commandLine.get_usage() + longExample;
        stk::parallel::print_and_exit(usage, comm);
    }
    STK_ThrowRequireMsg(state == stk::CommandLineParser::ParseComplete, stk::get_quick_error(execName, quickExample));
}

namespace parallel {

void print_and_exit(const std::string &msg, MPI_Comm comm)
{
    if(stk::parallel_machine_rank(comm) == 0)
        sierra::Env::outputP0() << msg << std::endl;
    stk::parallel_machine_finalize();
    std::exit(0);
}

}


}
