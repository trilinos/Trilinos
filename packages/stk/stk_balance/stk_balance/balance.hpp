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
 //
#ifndef STKBALANCE_HPP
#define STKBALANCE_HPP

namespace stk { namespace balance { class BalanceSettings; } }
namespace stk { namespace mesh { class BulkData; } }

#include <stk_mesh/base/Selector.hpp>
#include <vector>
#include <string>
#include "mpi.h"
#include <stk_balance/internal/balanceDefaults.hpp>

namespace stk
{
namespace balance
{

bool balanceStkMesh(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData);
bool balanceStkMesh(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData, const std::vector<stk::mesh::Selector>& selectors);
void run_stk_rebalance(const std::string& outputDirectory, const std::string& exodusFilename, stk::balance::AppTypeDefaults appType, MPI_Comm comm);
void run_stk_balance_with_settings(const std::string& outputDirectory, const std::string& exodusFilename, MPI_Comm comm, stk::balance::BalanceSettings& graphOptions);
void initial_decomp_and_balance(stk::mesh::BulkData &bulk,
                                stk::balance::BalanceSettings& graphOptions,
                                const std::string& exodusFilename,
                                const std::string& outputFilename);

}
}
#endif
