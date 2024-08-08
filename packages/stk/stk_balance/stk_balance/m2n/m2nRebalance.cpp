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
#include "m2nRebalance.hpp"
#include <stk_balance/m2n/M2NInputMesh.hpp>
#include <stk_balance/m2n/M2NOutputMesh.hpp>
#include <stk_balance/balanceUtils.hpp>
#include "stk_balance/internal/LogUtils.hpp"
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/FillMesh.hpp>
#include "stk_util/environment/Env.hpp"
#include "stk_util/environment/EnvData.hpp"
#include "stk_util/environment/OutputLog.hpp"

namespace stk {
namespace balance {
namespace m2n {

void set_output_streams(MPI_Comm comm, const stk::balance::M2NBalanceSettings & balanceSettings)
{
  if (stk::parallel_machine_rank(comm) == 0) {
    const std::string & logName = balanceSettings.get_log_filename();
    if (logName == "cout"  || logName == "cerr") {
      stk::EnvData::instance().m_outputP0 = stk::get_log_ostream(logName);
    }
    else {
      stk::bind_output_streams("log=\"" + logName + "\"");
      stk::EnvData::instance().m_outputP0 = stk::get_log_ostream("log");
    }
  }
  else {
    stk::EnvData::instance().m_outputP0 = &stk::EnvData::instance().m_outputNull;
  }

  Ioss::Utils::set_output_stream(sierra::Env::outputP0());
}

void print_running_message(const stk::balance::M2NBalanceSettings & balanceSettings, MPI_Comm comm)
{
  if (stk::parallel_machine_rank(comm) == 0) {
    std::ostream diag_stream(std::cout.rdbuf());
    stk::register_ostream(diag_stream, "diag_stream");

    const std::string & logName = balanceSettings.get_log_filename();
    const bool usingLogFile = not (logName == "cout" || logName == "cerr");
    if (usingLogFile) {
      stk::bind_output_streams("diag_stream>log");
      stk::bind_output_streams("diag_stream>+cout");
    }
    else {
      stk::bind_output_streams("diag_stream>" + logName);
    }

    const int inputRanks = stk::parallel_machine_size(comm);
    const int outputRanks = balanceSettings.get_num_output_processors();
    diag_stream << "stk_balance_m2n converting from " << inputRanks << " to " << outputRanks << " MPI ranks" << std::endl;
    if (usingLogFile) {
      diag_stream << "        Log file:  " << logName << std::endl;
    }

    diag_stream << "     Input files:  " << balanceSettings.get_input_filename();
    if (inputRanks > 1) {
      diag_stream << "." << inputRanks << ".*";
    }
    diag_stream << std::endl;

    diag_stream << "    Output files:  " << balanceSettings.get_input_filename();
    if (outputRanks > 1) {
      diag_stream << "." << outputRanks << ".*";
    }
    diag_stream << std::endl;

    stk::unregister_ostream(diag_stream);
  }
}

void m2nRebalance(stk::io::StkMeshIoBroker& ioBroker, const stk::balance::M2NBalanceSettings& balanceSettings)
{
  stk::balance::m2n::InputMesh inputMesh(ioBroker, balanceSettings);

  const std::vector<std::vector<unsigned>> targetSubdomainsForEachBatch = inputMesh.get_output_subdomains_for_each_batch();

  for (const std::vector<unsigned>& targetSubdomains : targetSubdomainsForEachBatch) {
    stk::balance::m2n::OutputMesh outputMesh(inputMesh, targetSubdomains);
    outputMesh.transfer_and_write();
  }
}

void rebalance_m2n(stk::balance::M2NBalanceSettings &balanceSettings, MPI_Comm comm)
{
  print_running_message(balanceSettings, comm);
  print_banner(sierra::Env::outputP0());

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(comm).create();
  stk::io::StkMeshIoBroker ioBroker;
  stk::io::fill_mesh_preexisting(ioBroker, balanceSettings.get_input_filename(), *bulk);

  stk::balance::m2n::m2nRebalance(ioBroker, balanceSettings);
}

}}}
