// Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
 // 
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
 //     * Neither the name of Sandia Corporation nor the names of its
 //       contributors may be used to endorse or promote products derived
 //       from this software without specific prior written permission.
 // 
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
