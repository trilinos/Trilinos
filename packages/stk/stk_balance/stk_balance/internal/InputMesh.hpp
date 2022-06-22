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
#ifndef INPUTMESH_HPP
#define INPUTMESH_HPP

#include <stk_mesh/base/Types.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace io { class StkMeshIoBroker; } }
namespace stk { namespace balance { class BalanceSettings; } }
namespace stk { namespace balance { class Decomposer; } }

namespace stk {
namespace balance {

class InputMesh
{
public:
  InputMesh(stk::io::StkMeshIoBroker& ioBroker,
            const stk::balance::BalanceSettings& balanceSettings);
  ~InputMesh();

  std::vector<std::vector<unsigned>> get_output_subdomains_for_each_batch() const;
  std::string get_subdomain_part_name(unsigned subdomainId) const;
  const stk::mesh::BulkData& get_bulk() const { return m_bulk; }
  stk::io::StkMeshIoBroker& get_io_broker() const { return m_ioBroker; }
  const std::vector<unsigned>& get_owner_for_each_final_subdomain() const { return m_ownerForEachFinalSubdomain; }
  unsigned get_num_output_processors() const;
  std::string get_output_file_name() const;
  int get_global_num_nodes() const { return m_globalNumNodes; }
  int get_global_num_elements() const { return m_globalNumElems; }
  stk::io::EntitySharingInfo get_node_sharing_info(unsigned mySubdomain,
                                                   const std::vector<unsigned>& subdomainsInBatch) const;
  const BalanceSettings & get_balance_settings() const { return m_balanceSettings; }

private:
  void initialize_mesh_counts();
  void compute_output_partition();
  void compute_owner_for_each_output_subdomain();
  void store_final_decomp_on_elements();
  void move_elements_into_output_subdomain_parts();
  void declare_all_output_subdomain_parts();
  stk::mesh::EntityVector get_nodes_shared_between_subdomains(int this_subdomain_index,
                                                              int other_subdomain_index) const;
  bool is_connected_to_element(stk::mesh::Entity entity);
  stk::mesh::EntityVector get_all_orphans();

  stk::io::StkMeshIoBroker& m_ioBroker;
  stk::mesh::BulkData& m_bulk;
  stk::mesh::MetaData& m_meta;
  const BalanceSettings& m_balanceSettings;
  Decomposer* m_decomposer;
  int m_globalNumNodes;
  int m_globalNumElems;
  stk::balance::DecompositionChangeList m_decomp;
  std::vector<unsigned> m_ownerForEachFinalSubdomain;
  stk::mesh::PartVector m_subdomainParts;
};

}
}

#endif // INPUTMESH_HPP
