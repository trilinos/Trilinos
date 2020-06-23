#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/command_line/CommandLineParser.hpp>
#include <stk_util/command_line/CommandLineParserParallel.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_util/registry/ProductRegistry.hpp>
#include <stk_util/environment/Env.hpp>
#include <stk_util/environment/FileUtils.hpp>
#include <stk_util/util/string_utils.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <fstream>
#include <string>

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
    ThrowRequireMsg(state == stk::CommandLineParser::ParseComplete, stk::get_quick_error(execName, quickExample));
}

namespace parallel {

void print_and_exit(const std::string &msg, MPI_Comm comm)
{
    if(stk::parallel_machine_rank(comm) == 0)
        sierra::Env::outputP0() << msg << std::endl;
    MPI_Finalize();
    std::exit(0);
}

}


}
