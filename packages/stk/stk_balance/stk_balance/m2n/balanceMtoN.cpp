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
#include "balanceMtoN.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>
#include <stk_balance/internal/entityDataToField.hpp>
#include <stk_balance/m2n/M2NDecomposer.hpp>

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/Comm.hpp>
#include "MxNutils.hpp"
#include "MtoNRebalancer.hpp"

namespace stk {
namespace balance {
namespace m2n {

using DecomposerPtr = std::shared_ptr<M2NDecomposer>;

DecomposerPtr make_decomposer(stk::mesh::BulkData& bulkData,
                              const stk::balance::M2NBalanceSettings& balanceSettings)
{
  DecomposerPtr decomposer;
  if (balanceSettings.get_use_nested_decomp()) {
    decomposer = std::make_shared<M2NDecomposerNested>(bulkData, balanceSettings);
  }
  else {
    decomposer = std::make_shared<M2NDecomposer>(bulkData, balanceSettings);
  }
  return decomposer;
}

bool rebalanceMtoN(stk::io::StkMeshIoBroker& ioBroker,
                   const stk::balance::M2NBalanceSettings & balanceSettings,
                   int numSteps,
                   double timeStep)
{
    DecomposerPtr decomposer = make_decomposer(ioBroker.bulk_data(), balanceSettings);
    MtoNRebalancer m2nRebalancer(ioBroker, *decomposer, balanceSettings);
    m2nRebalancer.rebalance(numSteps, timeStep);

    return true;
}
}}}
