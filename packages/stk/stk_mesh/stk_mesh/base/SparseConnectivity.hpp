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

#ifndef STK_MESH_SPARSECONNECTIVITY_HPP
#define STK_MESH_SPARSECONNECTIVITY_HPP

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>      // for EntityRank
#include <stk_mesh/base/Connectivity.hpp>
#include <stk_topology/topology.hpp>
#include <vector>

namespace stk {
namespace mesh {

class SparseConnectivity {
public:
  SparseConnectivity(unsigned numRanks = stk::topology::NUM_RANKS);
  virtual ~SparseConnectivity(){}

  void update_num_ranks(unsigned numRanks);
  void update_size_of_entity_index_space(unsigned numEntities);

  static bool has_permutation(EntityRank fromRank, EntityRank toRank);

  unsigned num_connectivity(Entity fromEntity, EntityRank connectedRank) const;

  const Entity* begin_connectivity(Entity fromEntity, EntityRank connectedRank) const;
  const Entity* end_connectivity(Entity fromEntity, EntityRank connectedRank) const;

  const ConnectivityOrdinal* begin_ordinals(Entity fromEntity, EntityRank connectedRank) const;
  const ConnectivityOrdinal* end_ordinals(Entity fromEntity, EntityRank connectedRank) const;

  const Permutation* begin_permutations(Entity fromEntity, EntityRank connectedRank) const;
  const Permutation* end_permutations(Entity fromEntity, EntityRank connectedRank) const;

  bool add_connectivity(EntityRank fromRank,
                        Entity fromEntity,
                        EntityRank connectedRank,
                        Entity connectedEntity,
                        ConnectivityOrdinal ordinal,
                        Permutation permutation);

  bool replace_or_add_connectivity(EntityRank fromRank,
                                   Entity fromEntity,
                                   EntityRank connectedRank,
                                   Entity connectedEntity,
                                   ConnectivityOrdinal ordinal);

  bool remove_connectivity(EntityRank fromRank,
                           Entity fromEntity,
                           EntityRank connectedRank,
                           Entity connectedEntity,
                           ConnectivityOrdinal ordinal);

  bool replace_permutation(EntityRank fromRank,
                           Entity fromEntity,
                           EntityRank connectedRank,
                           ConnectivityOrdinal ordinal,
                           Permutation newPermutation);

  PairIterEntity get_connectivity(Entity fromEntity, EntityRank connectedRank) const;

  const Connectivity& get_connectivity(Entity entity) const { return m_meshConnectivity[entity.local_offset()]; }
  Connectivity& get_connectivity(Entity entity) { return m_meshConnectivity[entity.local_offset()]; }

private:

  std::vector<Connectivity> m_meshConnectivity;
  std::vector<std::vector<unsigned>> m_entityPermutationOffset;
  std::vector<std::vector<std::vector<Permutation>>> m_permutations;
  unsigned m_numRanks;
};

inline
unsigned SparseConnectivity::num_connectivity(Entity fromEntity,
                                              EntityRank connectedRank) const
{ 
  return m_meshConnectivity[fromEntity.local_offset()].num_connectivity(connectedRank);
}

inline
const Entity* SparseConnectivity::begin_connectivity(Entity fromEntity,
                                                     EntityRank connectedRank) const
{
  return m_meshConnectivity[fromEntity.local_offset()].begin_connectivity(connectedRank);
}

inline
const Entity* SparseConnectivity::end_connectivity(Entity fromEntity,
                                                   EntityRank connectedRank) const
{
  return m_meshConnectivity[fromEntity.local_offset()].end_connectivity(connectedRank);
}

inline
PairIterEntity SparseConnectivity::get_connectivity(Entity fromEntity, EntityRank connectedRank) const
{
  return m_meshConnectivity[fromEntity.local_offset()].get_connectivity(connectedRank);
}

inline
const ConnectivityOrdinal* SparseConnectivity::begin_ordinals(Entity fromEntity,
                                                              EntityRank connectedRank) const
{
  return m_meshConnectivity[fromEntity.local_offset()].begin_ordinals(connectedRank);
}

inline
const ConnectivityOrdinal* SparseConnectivity::end_ordinals(Entity fromEntity,
                                                            EntityRank connectedRank) const
{
  return m_meshConnectivity[fromEntity.local_offset()].end_ordinals(connectedRank);
}

} // namespace mesh
} // namespace stk

#endif //STK_MESH_SPARSECONNECTIVITY_HPP

