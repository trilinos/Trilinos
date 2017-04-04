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

#ifndef stk_mesh_FEMHelpers_hpp
#define stk_mesh_FEMHelpers_hpp

#include <stddef.h>                     // for NULL
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, etc
#include <stk_mesh/base/Types.hpp>      // for PartVector, EntityId, etc
#include <vector>                       // for vector, etc
#include <stk_topology/topology.hpp>
#include "stk_mesh/base/Entity.hpp"     // for Entity
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Part; } }

namespace stk {
namespace mesh {

Entity declare_element( BulkData & mesh ,
                        PartVector & parts , // parts[0] expected to have topology
                        const EntityId elem_id ,
                        const EntityIdVector & node_ids );

inline
Entity declare_element( BulkData & mesh ,
                        Part & partWithTopology ,
                        const EntityId elem_id ,
                        const EntityIdVector & node_ids )
{
  PartVector vec(1, &partWithTopology);
  return declare_element(mesh, vec, elem_id, node_ids);
}

/** \brief  Create (or find) an element edge.
 *
 *  The element must be a member of a Part with a topology.
 */
Entity declare_element_edge( BulkData & mesh ,
			     const stk::mesh::EntityId global_side_id ,
			     Entity elem ,
			     const unsigned local_side_id ,
			     const stk::mesh::PartVector& parts = stk::mesh::PartVector());



/** \brief  Create (or find) an element side.
 *
 *  The element must be a member of a Part with a topology.
 */
Entity connect_side_to_element_with_ordinal( BulkData & mesh ,
                               Entity elem ,
                               Entity side ,
                               const unsigned local_side_id ,
                               stk::mesh::Part* part = NULL);

/** \brief finds oridinal and permutation of an entity relative to a parent entity
 *
 * This assumes parent is no higher rank than element and no less than edge and
 * that child is of less rank than parent.
 *
 *
 */
typedef std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> OrdinalAndPermutation;

OrdinalAndPermutation get_ordinal_and_permutation(const stk::mesh::BulkData& mesh,
                                                  stk::mesh::Entity parent_entity,
                                                  stk::mesh::EntityRank to_rank,
                                                  const stk::mesh::EntityVector &nodes_of_sub_rank);

std::pair<bool, unsigned> sub_rank_equivalent(const stk::mesh::BulkData& mesh, stk::mesh::Entity element, unsigned ordinal, stk::mesh::EntityRank subRank,
                                                            const stk::mesh::Entity* subRankNodes);

std::pair<bool, unsigned> side_equivalent(const stk::mesh::BulkData& mesh, stk::mesh::Entity element, unsigned sideOrdinal, const stk::mesh::Entity* candidateSideNodes);

bool is_side_equivalent(const stk::mesh::BulkData& mesh, stk::mesh::Entity element, unsigned sideOrdinal, const stk::mesh::Entity* candidateSideNodes);

bool is_edge_equivalent(const stk::mesh::BulkData& mesh, stk::mesh::Entity element, unsigned edgeOrdinal, const stk::mesh::Entity* candidateEdgeNodes);

NAMED_PAIR(EquivAndPositive, bool, is_equiv, bool, is_positive)

EquivAndPositive is_side_equivalent_and_positive(const stk::mesh::BulkData& mesh, stk::mesh::Entity element, unsigned sideOrdinal, const stk::mesh::Entity* candidateSideNodes);

EquivAndPositive is_equivalent_and_positive(const stk::mesh::BulkData& mesh, stk::mesh::Entity element, unsigned ordinal, stk::mesh::EntityRank subRank, const stk::mesh::Entity* candidateNodes);
/**
 * Given an entity, subcell_rank, and subcell_id, return the nodes
 * that make up the subcell in a correct order for the given polarity.
 *
 * \param entity
 * \param subcell_rank
 * \param subcell_identifier
 * \param subcell_nodes EntityVector output of the subcell nodes
 * \return topology of the requested subcell
 */
stk::topology get_subcell_nodes(const BulkData& mesh,
    const Entity entity ,
    EntityRank         subcell_rank ,
    unsigned           subcell_identifier ,
    EntityVector     & subcell_nodes
    );

/** \brief  Given an entity and collection of nodes, return the
 *          local id of the subcell that contains those nodes in the
 *          correct orientation.
 */
int get_entity_subcell_id( const BulkData& mesh, const Entity entity ,
                           const EntityRank          subcell_rank,
                           stk::topology side_topology,
                           const EntityVector      & side_nodes );

inline
void get_parts_with_topology(stk::topology topology,
                             stk::mesh::BulkData& mesh,
                             stk::mesh::PartVector& parts,
                             bool skip_topology_root_parts=false)
{
  parts.clear();

  const stk::mesh::MetaData & fem_meta = mesh.mesh_meta_data();

  const stk::mesh::PartVector& all_parts = fem_meta.get_parts();

  stk::mesh::PartVector::const_iterator
    iter = all_parts.begin(),
    iter_end = all_parts.end();

  for(; iter!=iter_end; ++iter) {
    stk::mesh::Part* part =  *iter;
    if (fem_meta.get_topology(*part) == topology) {
      if (skip_topology_root_parts && stk::mesh::is_topology_root_part(*part)) {
        continue;
      }
      parts.push_back(part);
    }
  }
}

stk::mesh::Entity get_side_entity_for_elem_side_pair(const stk::mesh::BulkData &bulk, Entity elem, int sideOrdinal);
stk::mesh::Entity get_side_entity_for_elem_id_side_pair_of_rank(const stk::mesh::BulkData &bulk, int64_t elemId, int sideOrdinal, stk::mesh::EntityRank sideRank);

/** \---} */

} //namespace mesh
} //namespace stk
#endif
