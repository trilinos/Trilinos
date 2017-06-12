#ifndef STK_UTIL_ENVIRONMENT_COMMANDLINEPARSERUTILS_HPP
#define STK_UTIL_ENVIRONMENT_COMMANDLINEPARSERUTILS_HPP

#include <stk_util/environment/FileUtils.hpp>
#include "CommandLineParser.hpp"
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/registry/ProductRegistry.hpp>

namespace stk {

class CommandLineParserParallel : public CommandLineParser
{
public:
    CommandLineParserParallel(MPI_Comm c) : CommandLineParser(), comm(c) {}
    explicit CommandLineParserParallel(const std::string &usagePreamble, MPI_Comm c) : CommandLineParser(usagePreamble), comm(c) {}
    virtual void print_message(const std::string &msg)
    {
        if(stk::parallel_machine_rank(comm) == 0)
            CommandLineParser::print_message(msg);
    }
protected:
    MPI_Comm comm;
};

}

#endif //STK_UTIL_ENVIRONMENT_COMMANDLINEPARSERUTILS_HPP
