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


#ifndef stk_mesh_impl_GetOrdPermImpl_hpp
#define stk_mesh_impl_GetOrdPermImpl_hpp

#include <stk_util/stk_config.h>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace stk {
namespace mesh {

typedef std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> OrdinalAndPermutation;

namespace impl {

typedef std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ConnectivityAndOrdinal;


inline
void mark_if_negative_permutation(stk::EquivalentPermutation &result, stk::topology sub_topology)
{
    if (result.is_equivalent && result.permutation_number >= sub_topology.num_positive_permutations())
    {
        result.is_equivalent = false;
    }
}

inline
void set_ordinal_and_permutation_if_equivalent(const stk::EquivalentPermutation &result,
                                               unsigned ordinal,
                                               ConnectivityAndOrdinal &ordinalAndPermutation)
{
    if (result.is_equivalent == true)
    {
        ordinalAndPermutation.first = static_cast<stk::mesh::ConnectivityOrdinal>(ordinal);
        ordinalAndPermutation.second = static_cast<stk::mesh::Permutation>(result.permutation_number);
    }
}

class ShellPermutationFilter
{
public:
    ShellPermutationFilter(const stk::mesh::BulkData& mesh,
                           stk::mesh::Entity parent_entity,
                           stk::mesh::EntityRank to_rank)
    : m_mesh(mesh), m_entity(parent_entity), m_toRank(to_rank), m_filterForShell(false)    {
        stk::topology elemTopology = m_mesh.bucket(m_entity).topology();
        m_filterForShell = elemTopology.is_shell() && to_rank == mesh.mesh_meta_data().side_rank();
    }

    ~ShellPermutationFilter() {}

    template<typename NodeArrayType>
    bool set_ordinal_and_permutation(const NodeArrayType& nodes_of_sub_rank,
                                     stk::mesh::Entity nodes_of_sub_topology[],
                                     stk::topology sub_topology,
                                     unsigned ordinal,
                                     ConnectivityAndOrdinal &ordinalAndPermutation) const
    {     
        stk::EquivalentPermutation result = sub_topology.is_equivalent(nodes_of_sub_rank.data(), nodes_of_sub_topology);

        constexpr unsigned numLegacyShellSidesTo_TEMPORARILY_PreservePreviousBehavior = 2;
        if(m_filterForShell && 
           ordinal < numLegacyShellSidesTo_TEMPORARILY_PreservePreviousBehavior) {
            mark_if_negative_permutation(result, sub_topology);
        }

        set_ordinal_and_permutation_if_equivalent(result, ordinal, ordinalAndPermutation);
        return result.is_equivalent;
    }

private:
    ShellPermutationFilter();

    const stk::mesh::BulkData& m_mesh;
    stk::mesh::Entity m_entity;
    stk::mesh::EntityRank m_toRank;
    bool m_filterForShell;
};

template<typename PermutationFilter, typename NodeArrayType>
OrdinalAndPermutation
get_ordinal_and_permutation_with_filter(const stk::mesh::BulkData& mesh,
                                        stk::mesh::Entity parent_entity,
                                        stk::mesh::EntityRank to_rank, 
                                        const NodeArrayType& nodes_of_sub_rank,
                                        PermutationFilter &pFilter)
{
  ConnectivityAndOrdinal ordinalAndPermutation = std::make_pair(stk::mesh::INVALID_CONNECTIVITY_ORDINAL, stk::mesh::INVALID_PERMUTATION);

  unsigned nodes_of_sub_rank_size = nodes_of_sub_rank.size();

  const Entity* elemNodes = mesh.begin_nodes(parent_entity);
  stk::topology elemTopology = mesh.bucket(parent_entity).topology();
  unsigned num_entities_of_sub_topology = elemTopology.num_sub_topology(to_rank);
  const unsigned max_nodes_possible = 100;
  stk::mesh::Entity nodes_of_sub_topology[max_nodes_possible];

  for (unsigned i=0;i<num_entities_of_sub_topology;++i) {
    stk::topology sub_topology = elemTopology.sub_topology(to_rank, i);
    unsigned num_nodes = sub_topology.num_nodes();

    if (num_nodes !=  nodes_of_sub_rank_size) {
      continue;
    }

    STK_ThrowRequireMsg(num_nodes<=max_nodes_possible, "Program error. Exceeded expected array dimensions. Contact sierra-help for support.");
    elemTopology.sub_topology_nodes(elemNodes, to_rank, i, nodes_of_sub_topology);

    pFilter.set_ordinal_and_permutation(nodes_of_sub_rank, nodes_of_sub_topology, sub_topology, i, ordinalAndPermutation);
    if (ordinalAndPermutation.first != stk::mesh::INVALID_CONNECTIVITY_ORDINAL) {
      break;
    }
  }

  return ordinalAndPermutation;
}

} // namespace impl
} // namespace mesh
} // namespace stk

#endif
