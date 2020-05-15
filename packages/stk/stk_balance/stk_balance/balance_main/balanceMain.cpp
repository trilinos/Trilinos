 // Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#include "mpi.h"
#include <stk_balance/balance.hpp>

#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/Inputs.hpp>
#include <stk_balance/internal/balanceCommandLine.hpp>
#include <stk_balance/internal/balanceDefaults.hpp>

#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_util/environment/FileUtils.hpp>
#include <stk_util/util/string_utils.hpp>
#include <Kokkos_Core.hpp>

#include <string>
#include <iostream>
#include <fstream>

#include "stk_util/environment/EnvData.hpp"

int main(int argc, const char**argv)
{
    MPI_Init(&argc, const_cast<char***>(&argv));
    MPI_Comm comm = MPI_COMM_WORLD;

    if(stk::parallel_machine_rank(comm) != 0)
    {
        stk::EnvData::instance().m_outputP0 = &stk::EnvData::instance().m_outputNull;
    }

    Kokkos::initialize(argc, const_cast<char**>(argv));

    std::string execName = stk::tailname(argv[0]);
    stk::balance::CommandLineOptions cmdLineOptions;
    std::string quickExample = stk::balance::get_quick_example(execName, cmdLineOptions.infile.name, cmdLineOptions.outputDirectory.name, comm);
    stk::balance::ParsedOptions balanceOptions;
    try
    {
        balanceOptions = stk::balance::parse_balance_command_line(argc,
                                                                  argv,
                                                                  execName,
                                                                  comm);
    }
    catch(const std::exception &e)
    {
        std::string errorMessage = e.what() + stk::get_quick_error(execName, quickExample);
        stk::parallel::print_and_exit(errorMessage, comm);
    }

    stk::parallel::require_file_exists(balanceOptions.m_inFile, execName, quickExample, comm);

    stk::balance::print_running_msg(execName, balanceOptions, comm);
    try {
        stk::balance::run_stk_rebalance(balanceOptions, comm);
    }
    catch(std::exception& e)
    {
        std::cerr<<e.what()<<std::endl;
    }

    Kokkos::finalize_all();

    MPI_Finalize();
    return 0;
}


