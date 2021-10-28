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

#include <stk_mesh/base/SparseConnectivity.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <algorithm>

namespace stk {
namespace mesh {

SparseConnectivity::SparseConnectivity(unsigned numRanks)
 : m_numRanks(numRanks)
{
  update_num_ranks(numRanks);
}

void SparseConnectivity::update_num_ranks(unsigned numRanks)
{
  const unsigned minNumRanks = stk::topology::NUM_RANKS;
  m_numRanks = std::max(numRanks, minNumRanks);

  m_entityPermutationOffset.resize(m_numRanks);
  m_permutations.resize(m_numRanks);
  for(unsigned rank=0; rank<m_numRanks; ++rank) {
    m_entityPermutationOffset[rank].clear();
    m_permutations[rank].clear();
    m_permutations[rank].push_back(std::vector<Permutation>());
  }
}

void SparseConnectivity::update_size_of_entity_index_space(unsigned numEntities)
{
  m_meshConnectivity.resize(numEntities);
}

bool SparseConnectivity::has_permutation(EntityRank fromRank, EntityRank toRank)
{
  return (fromRank > stk::topology::NODE_RANK && fromRank < stk::topology::CONSTRAINT_RANK)
      && (toRank   > stk::topology::NODE_RANK && toRank   < stk::topology::CONSTRAINT_RANK);
}

const Permutation* SparseConnectivity::begin_permutations(Entity fromEntity,
                                                          EntityRank connectedRank) const
{
  if (static_cast<unsigned>(connectedRank) < m_entityPermutationOffset.size() &&
      fromEntity.local_offset() < m_entityPermutationOffset[connectedRank].size()) {
    const int entityOffset = m_entityPermutationOffset[connectedRank][fromEntity.local_offset()];
    return m_permutations[connectedRank][entityOffset].data();
  }
  return nullptr;
}

const Permutation* SparseConnectivity::end_permutations(Entity fromEntity,
                                                        EntityRank connectedRank) const
{
  if (static_cast<unsigned>(connectedRank) < m_entityPermutationOffset.size() &&
      fromEntity.local_offset() < m_entityPermutationOffset[connectedRank].size()) {
    const int entityOffset = m_entityPermutationOffset[connectedRank][fromEntity.local_offset()];
    const std::vector<Permutation>& perms = m_permutations[connectedRank][entityOffset];
    return (perms.data()+perms.size());
  }
  return nullptr;
}

bool SparseConnectivity::add_connectivity(EntityRank fromRank,
                                          Entity fromEntity,
                                          EntityRank connectedRank,
                                          Entity connectedEntity,
                                          ConnectivityOrdinal ordinal,
                                          Permutation permutation)
{
  Connectivity& conn = m_meshConnectivity[fromEntity.local_offset()];
  int insertLocation = conn.add_connectivity(connectedRank,
           connectedEntity, ordinal, (connectedRank > fromRank));

  if (insertLocation >= 0 && has_permutation(fromRank, connectedRank)) {
    if (fromEntity.local_offset() >= m_entityPermutationOffset[connectedRank].size()) {
      m_entityPermutationOffset[connectedRank].resize(fromEntity.local_offset()+1);
      m_entityPermutationOffset[connectedRank][fromEntity.local_offset()] = 0;
    }
    unsigned& permOffset = m_entityPermutationOffset[connectedRank][fromEntity.local_offset()];
    if (permOffset == 0) {
      permOffset = m_permutations[connectedRank].size();
      m_permutations[connectedRank].push_back(std::vector<Permutation>(1, permutation));
      return true;
    }
    std::vector<Permutation>& permutations = m_permutations[connectedRank][permOffset];
    permutations.insert(permutations.begin()+insertLocation, permutation);
  }

  return insertLocation >= 0;
}

bool SparseConnectivity::replace_or_add_connectivity(EntityRank fromRank,
                                                     Entity fromEntity,
                                                     EntityRank connectedRank,
                                                     Entity connectedEntity,
                                                     ConnectivityOrdinal ordinal)
{
  Connectivity& conn = m_meshConnectivity[fromEntity.local_offset()];
  int insertLocation = conn.replace_or_add_connectivity(connectedRank,
           connectedEntity, ordinal);

  return insertLocation >= 0;
}

bool SparseConnectivity::remove_connectivity(EntityRank fromRank,
                                             Entity fromEntity,
                                             EntityRank connectedRank,
                                             Entity connectedEntity,
                                             ConnectivityOrdinal ordinal)
{
  Connectivity& conn = m_meshConnectivity[fromEntity.local_offset()];
  int location = conn.remove_connectivity(connectedRank, connectedEntity, ordinal);

  if (location >= 0) {
    if (has_permutation(fromRank, connectedRank)) {
      const unsigned permOffset = m_entityPermutationOffset[connectedRank][fromEntity.local_offset()];
      std::vector<Permutation>& permutations = m_permutations[connectedRank][permOffset];
      permutations.erase(permutations.begin()+location);
    }
    return true;
  }
  return false;
}

bool SparseConnectivity::replace_permutation(EntityRank fromRank,
                                             Entity fromEntity,
                                             EntityRank connectedRank,
                                             ConnectivityOrdinal ordinal,
                                             Permutation newPermutation)
{
  if (fromEntity.local_offset() < m_entityPermutationOffset[connectedRank].size()) {
    const ConnectivityOrdinal* ordsBegin = begin_ordinals(fromEntity, connectedRank);
    const ConnectivityOrdinal* ordsEnd = end_ordinals(fromEntity, connectedRank);

    bool foundIt = false;
    const bool downward = fromRank > connectedRank;
    unsigned idx = ordsEnd - ordsBegin;
    if (downward) {
      const ConnectivityOrdinal* iter = std::lower_bound(ordsBegin, ordsEnd, ordinal);
      idx = iter - ordsBegin;
      foundIt = (iter != ordsEnd && *iter == ordinal);
    }
    else {
      const ConnectivityOrdinal* iter = std::find(ordsBegin, ordsEnd, ordinal);
      idx = iter - ordsBegin;
      foundIt = (iter != ordsEnd);
    }
    if (foundIt) {
      const int entityPermOffset = m_entityPermutationOffset[connectedRank][fromEntity.local_offset()];
      std::vector<Permutation>& perms = m_permutations[connectedRank][entityPermOffset];
      perms[idx] = newPermutation;
      return true;
    }
  }
  return false;
}

} // namespace mesh
} // namespace stk

