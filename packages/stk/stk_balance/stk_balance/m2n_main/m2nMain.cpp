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
//

#include <stk_util/parallel/Parallel.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/m2n/m2nRebalance.hpp>
#include <stk_balance/setup/M2NParser.hpp>
#include "stk_balance/internal/LogUtils.hpp"

#include <stk_util/environment/Env.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/util/human_bytes.hpp>

int main(int argc, const char**argv)
{
    MPI_Init(&argc, const_cast<char***>(&argv));
    MPI_Comm comm = MPI_COMM_WORLD;
    stk::balance::initialize_environment(comm, argv);

    stk::balance::M2NBalanceSettings balanceSettings;

    stk::balance::M2NParser parser(comm);
    parser.parse_command_line_options(argc, argv, balanceSettings);

    stk::balance::m2n::set_output_streams(comm, balanceSettings);
    stk::balance::m2n::rebalance_m2n(balanceSettings, comm);

    size_t hwmMax = 0, hwmMin = 0, hwmAvg = 0;
    stk::get_memory_high_water_mark_across_processors(comm, hwmMax, hwmMin, hwmAvg);
    sierra::Env::outputP0() << "Memory HWM across procs, max/min/avg: "
                            << stk::human_bytes(hwmMax) << " / "
                            << stk::human_bytes(hwmMin) << " / "
                            << stk::human_bytes(hwmAvg) << std::endl;
    MPI_Finalize();
    return 0;
}
