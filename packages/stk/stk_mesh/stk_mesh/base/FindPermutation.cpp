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

#include <stk_mesh/base/FindPermutation.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace stk {
namespace mesh {

Permutation find_permutation(const BulkData& bulk,
                             const stk::topology& entityTopology,
                             const Entity* entityNodes,
                             const stk::topology& sideTopology,
                             const Entity* sideNodes,
                             unsigned sideOrdinal)
{
  Entity expectedNodes[100];
  switch (sideTopology.rank())
  {    
  case stk::topology::EDGE_RANK:
      entityTopology.edge_nodes(entityNodes, sideOrdinal, expectedNodes);
      break;
  case stk::topology::FACE_RANK:
      entityTopology.face_nodes(entityNodes, sideOrdinal, expectedNodes);
      break;
  default:
      return INVALID_PERMUTATION;
  }

  stk::EquivalentPermutation equivPerm = sideTopology.is_equivalent(expectedNodes, sideNodes);
  return equivPerm.is_equivalent ? static_cast<Permutation>(equivPerm.permutation_number) : INVALID_PERMUTATION;
}

bool check_permutation(const BulkData& bulk,
                       Entity entity,
                       Entity subEntity,
                       unsigned subOrdinal,
                       Permutation expectedPerm)
{
  const stk::topology entityTopo = bulk.mesh_index(entity).bucket->topology();
  const stk::topology subTopo    = bulk.mesh_index(subEntity).bucket->topology();
  const Entity* entityNodes      = bulk.begin_nodes(entity);
  const Entity* subEntityNodes   = bulk.begin_nodes(subEntity);

  Permutation computedPerm = find_permutation(bulk, entityTopo, entityNodes, subTopo, subEntityNodes, subOrdinal);

  return computedPerm == expectedPerm;
}

}} // namespace stk::mesh

