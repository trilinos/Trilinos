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

#include "MeshFixtureRebalance.hpp"
#include <stk_util/parallel/OutputStreams.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_balance/rebalance.hpp>
#include <vector>
#include <stk_unit_test_utils/BuildMesh.hpp>

namespace {
using stk::unit_test_util::build_mesh;

class RebalanceFileOutput : public MeshFixtureRebalance
{
public:
  void rebalance_mesh(int numFinalProcs, const std::string & decompMethod = "rcb")
  {
    m_balanceSettings.set_is_rebalancing(true);
    m_balanceSettings.set_output_filename(get_output_file_name());
    m_balanceSettings.set_num_input_processors(stk::parallel_machine_size(get_comm()));
    m_balanceSettings.set_num_output_processors(numFinalProcs);
    m_balanceSettings.setDecompMethod(decompMethod);

    stk::set_outputP0(&stk::outputNull());
    stk::balance::rebalance(m_ioBroker, m_balanceSettings);
    stk::reset_default_output_streams();
  }
};

std::vector<std::pair<stk::mesh::EntityId, int>> getSharingInfo(stk::mesh::BulkData& bulkData)
{
  stk::mesh::EntityVector sharedNodes;
  const bool sortById = true;
  stk::mesh::get_entities(bulkData, stk::topology::NODE_RANK, bulkData.mesh_meta_data().globally_shared_part(), sharedNodes, sortById);
  std::vector<std::pair<stk::mesh::EntityId, int>> nodeSharingInfo;
  nodeSharingInfo.reserve(8*sharedNodes.size());

  std::vector<int> sharingProcs;
  for(stk::mesh::Entity sharedNode : sharedNodes)
  {
    bulkData.comm_shared_procs(bulkData.entity_key(sharedNode), sharingProcs);
    for(unsigned j=0;j<sharingProcs.size();++j)
      nodeSharingInfo.push_back(std::make_pair(bulkData.identifier(sharedNode), sharingProcs[j]));
  }

  return nodeSharingInfo;
}

void verify_node_sharing_info(const std::vector<std::pair<stk::mesh::EntityId, int>> &nodeSharingInfo, const std::string& filename)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3, MPI_COMM_WORLD);
  stk::io::fill_mesh(filename, *bulk);

  std::vector<std::pair<stk::mesh::EntityId, int>> nodeSharingInfoAfter = getSharingInfo(*bulk);

  EXPECT_TRUE(nodeSharingInfo == nodeSharingInfoAfter);
}

TEST_F(RebalanceFileOutput, CheckSharingInformation)
{
  if (stk::parallel_machine_size(get_comm()) != 4) return;

  setup_initial_mesh("1x4x4");
  std::vector<std::pair<stk::mesh::EntityId, int>> nodeSharingInfo = getSharingInfo(get_bulk());

  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(get_bulk(), counts);
  int global_num_nodes = counts[stk::topology::NODE_RANK];
  int global_num_elems = counts[stk::topology::ELEM_RANK];

  const std::string outputFilename = "TemporaryOutputFile.g";
  stk::io::OutputParams params(get_bulk());
  stk::io::write_file_for_subdomain(outputFilename,
                                    get_bulk().parallel_rank(),
                                    get_bulk().parallel_size(),
                                    global_num_nodes,
                                    global_num_elems,
                                    params,
                                    nodeSharingInfo);

  verify_node_sharing_info(nodeSharingInfo, outputFilename);

  clean_up_temporary_files();
    
  if (get_parallel_rank() == 0) {
    for (int i = 0; i < get_parallel_size(); ++i) {
      std::string suffix = "." + std::to_string(get_parallel_size()) + "." + std::to_string(i);
      std::string filename = outputFilename + suffix;
      std::string inputFilename = get_input_file_name() + suffix;
      unlink(filename.c_str());
      unlink(inputFilename.c_str());
    }
  }
}

}
