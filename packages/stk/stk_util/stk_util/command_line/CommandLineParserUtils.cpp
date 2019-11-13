#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/command_line/CommandLineParser.hpp>
#include <stk_util/command_line/CommandLineParserParallel.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_util/registry/ProductRegistry.hpp>
#include <stk_util/environment/FileUtils.hpp>
#include <stk_util/util/string_utils.hpp>
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
    stk::parallel::require(state == stk::CommandLineParser::ParseComplete, stk::get_quick_error(execName, quickExample), comm);
}

namespace parallel {

void print_and_exit(const std::string &msg, MPI_Comm comm)
{
    if(stk::parallel_machine_rank(comm) == 0)
        std::cerr << msg << std::endl;
    MPI_Finalize();
    std::exit(0);
}

void require(bool requirement, const std::string &msg, MPI_Comm comm)
{
    if(!requirement)
        print_and_exit(msg, comm);
}

bool does_file_exist(const std::string& filename)
{
    bool exists = true;
    if(!std::ifstream(filename))
        exists = false;
    return exists;
}

void require_file_exists(const std::string& inFile, const std::string& execName, const std::string& quickExample, MPI_Comm comm)
{
    require(does_file_exist(inFile),
                           "Error: input file does not exist.\n" + stk::get_quick_error(execName, quickExample),
                           comm);
}

}


}
