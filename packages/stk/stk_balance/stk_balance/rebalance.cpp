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
#include "rebalance.hpp"
#include "stk_balance/internal/InputMesh.hpp"
#include "stk_balance/internal/OutputMesh.hpp"
#include "stk_balance/internal/privateDeclarations.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include <vector>

namespace stk {
namespace balance {

void rebalance(stk::io::StkMeshIoBroker & ioBroker, const stk::balance::BalanceSettings & balanceSettings)
{
  stk::balance::InputMesh inputMesh(ioBroker, balanceSettings);

  const std::vector<std::vector<unsigned>> targetSubdomainsForEachBatch = inputMesh.get_output_subdomains_for_each_batch();

  for (const std::vector<unsigned>& targetSubdomains : targetSubdomainsForEachBatch) {
    stk::balance::OutputMesh outputMesh(inputMesh, targetSubdomains);
    outputMesh.transfer_and_write();
  }

  internal::logMessage(ioBroker.bulk_data().parallel(), "Finished writing output mesh");
}

}
}
