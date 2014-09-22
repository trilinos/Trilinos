/*------------------------------------------------------------------------*/
/*                 Copyright (c) 2013, Sandia Corporation.
/*                 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*                 the U.S. Governement retains certain rights in this software.
/*                 
/*                 Redistribution and use in source and binary forms, with or without
/*                 modification, are permitted provided that the following conditions are
/*                 met:
/*                 
/*                     * Redistributions of source code must retain the above copyright
/*                       notice, this list of conditions and the following disclaimer.
/*                 
/*                     * Redistributions in binary form must reproduce the above
/*                       copyright notice, this list of conditions and the following
/*                       disclaimer in the documentation and/or other materials provided
/*                       with the distribution.
/*                 
/*                     * Neither the name of Sandia Corporation nor the names of its
/*                       contributors may be used to endorse or promote products derived
/*                       from this software without specific prior written permission.
/*                 
/*                 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*                 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*                 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*                 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*                 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*                 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*                 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*                 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*                 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*                 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*                 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*                 
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_mesh/base/BoundaryAnalysis.hpp>
#include <stddef.h>                     // for NULL, size_t
#include <algorithm>                    // for sort, binary_search, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, EntityLess
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/FEMHelpers.hpp>  // for get_entity_subcell_id, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, get_cell_topology
#include <stk_mesh/base/Part.hpp>       // for Part
#include <vector>                       // for vector, etc
#include "Shards_CellTopologyData.h"    // for CellTopologyData
#include "stk_mesh/base/CellTopology.hpp"  // for CellTopology
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Relation.hpp"
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator|
#include "stk_mesh/base/Types.hpp"      // for EntityVector, EntityRank
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowAssert
#include "stk_util/util/NamedPair.hpp"  // for EntitySide::first_type, etc


namespace stk {
namespace mesh {

namespace {

void filter_superimposed_entities(const BulkData& mesh, const Entity entity, EntityVector & entities)
{
  // Get the node entities for the nodes that make up the entity, we'll
  // use this to check for superimposed entities
  EntityVector entity_nodes;
  entity_nodes.reserve(mesh.num_nodes(entity));
  Entity const *rels_itr = mesh.begin_nodes(entity);
  Entity const *rels_end = mesh.end_nodes(entity);
  for ( ; rels_itr != rels_end; ++rels_itr ) {
    entity_nodes.push_back(*rels_itr);
  }
  std::sort(entity_nodes.begin(), entity_nodes.end());

  // Make sure to remove the all superimposed entities from the list
  unsigned num_nodes_in_orig_entity = entity_nodes.size();
  EntityVector current_nodes;
  current_nodes.resize(num_nodes_in_orig_entity);
  EntityVector::iterator itr = entities.begin();
  while ( itr != entities.end() )
  {
    Entity current_entity = *itr;
    size_t num_node_rels = mesh.num_nodes(current_entity);

    if (current_entity == entity) {
      // Superimposed with self by definition
      itr = entities.erase(itr);
    }
    else if (
        num_node_rels != num_nodes_in_orig_entity
        ) {
      // current_entity has a different number of nodes than entity, they
      // cannot be superimposed
      ++itr;
    }
    else {
      Entity const *node_rels = mesh.begin_nodes(current_entity);
      for (unsigned i = 0; i < num_node_rels; ++i)
      {
        current_nodes[i] = node_rels[i];
      }
      std::sort(current_nodes.begin(), current_nodes.end());

      bool entities_are_superimposed = entity_nodes == current_nodes;
      if (entities_are_superimposed) {
        itr = entities.erase(itr);
      }
      else {
        ++itr;
      }
    }
  }
}

/** \brief  Get the entities adjacent to the input entity.
 *
 *  The adjacent entities are of the same rank as the input entity.
 *  Adjacency is defined by the input entity sharing a common
 *  sub-cell with the adjacent entities.
 *
 *  subcell_rank defines the rank of the (potentially) common subcell
 *  subcell_identifier defined the local id of the common subcell
 *  adjacent_entities is an output parameter that contains pairs that
 *     have the adjacent entity and the local id of the common subcell
 *     with respect to the adjacent entity.
 */
void get_adjacent_entities( const BulkData& mesh,  const Entity entity ,
                            EntityRank subcell_rank ,
                            unsigned subcell_identifier ,
                            std::vector< EntitySideComponent> & adjacent_entities)
{
  adjacent_entities.clear();

  // Get nodes that make up the subcell we're looking at
  EntityVector subcell_nodes;
  const CellTopologyData * subcell_topology = get_subcell_nodes(mesh, entity,
                                                                subcell_rank,
                                                                subcell_identifier,
                                                                subcell_nodes);
  ThrowAssert(subcell_topology != NULL);


  // Given the nodes related to the subcell, find all entities
  // with the same rank that have a relation to all of these nodes
  EntityVector potentially_adjacent_entities;

  get_entities_through_relations(mesh, subcell_nodes,
                                 mesh.entity_rank(entity),
                                 potentially_adjacent_entities);

  // We don't want to include entities that are superimposed with
  // the input entity
  filter_superimposed_entities(mesh, entity, potentially_adjacent_entities);

  // Add the local ids, from the POV of the adj entitiy, to the return value.
  // Reverse the nodes so that the adjacent entity has them in the positive
  // orientation
  std::reverse(subcell_nodes.begin(),subcell_nodes.end());

  for (EntityVector::const_iterator eitr = potentially_adjacent_entities.begin();
       eitr != potentially_adjacent_entities.end(); ++eitr) {
    int local_subcell_num = get_entity_subcell_id(mesh, *eitr,
                                                  subcell_rank,
                                                  *subcell_topology,
                                                  subcell_nodes);
    if ( local_subcell_num != -1) {
      adjacent_entities.push_back(EntitySideComponent(*eitr, local_subcell_num));
    }
  }
}

} // unnamed namespace

void boundary_analysis(const BulkData& bulk_data,
                       const EntityVector & entities_closure,
                       EntityRank closure_rank,
                       EntitySideVector& boundary)
{
  const Selector locally_used = MetaData::get(bulk_data).locally_owned_part()
                              | MetaData::get(bulk_data).globally_shared_part();

  // find an iterator that points to the last item in the closure that is of a
  // lower-order than the closure_rank
  EntityVector::const_iterator itr =
    std::lower_bound(entities_closure.begin(),
                     entities_closure.end(),
                     EntityKey(closure_rank, 0),
                     EntityLess(bulk_data));

  // iterate over all the entities in the closure up to the iterator we computed above
  for ( ; itr != entities_closure.end() && bulk_data.entity_rank(*itr) == closure_rank; ++itr) {
    // some temporaries for clarity
    std::vector<EntitySideComponent > adjacent_entities;
    Entity curr_entity = *itr;
    const CellTopologyData* celltopology = get_cell_topology(bulk_data.bucket(curr_entity)).getCellTopologyData();
    if (celltopology == NULL) {
      continue;
    }

    EntityRank subcell_rank = closure_rank == stk::topology::ELEMENT_RANK ? bulk_data.mesh_meta_data().side_rank() : static_cast<EntityRank>(closure_rank - 1);

    // iterate over the subcells of the current entity
    for (unsigned nitr = 0; nitr < celltopology->subcell_count[subcell_rank]; ++nitr) {
      // find the entities (same rank as subcell) adjacent to this subcell
      unsigned subcell_identifier = nitr;
      get_adjacent_entities(bulk_data, curr_entity,
                            subcell_rank,
                            subcell_identifier,
                            adjacent_entities);

      // try to figure out if we want to keep ((curr_entity, subcell_identifier),
      //                                       (adjacent_entity, [k]))
      // if so, push into boundary

      // it is a keeper if adjacent entities[k] is not in the entities closure
      // AND if either curr_entity OR adjacent entities[k] is a member of the
      //               bulk_data.locally_used

      if (adjacent_entities.empty()) {
        EntitySide keeper;
        keeper.inside.entity = curr_entity;
        keeper.inside.side_ordinal = subcell_identifier;
        keeper.outside.entity = Entity();
        keeper.outside.side_ordinal = 0;
        boundary.push_back(keeper);
        continue;
      }

      // iterate over adjacent entities (our neighbors)
      for (std::vector<EntitySideComponent >::const_iterator
           adj_itr = adjacent_entities.begin();
           adj_itr != adjacent_entities.end(); ++adj_itr) {
        // grab a reference to this neighbor for clarity
        const Entity neighbor = adj_itr->entity;

        // see if this neighbor is in the closure, if so, not a keeper
        bool neighbor_is_in_closure =
          std::binary_search(entities_closure.begin(),
                             entities_closure.end(),
                             neighbor,
                             EntityLess(bulk_data));

        if (neighbor_is_in_closure) {
          continue;
        }

        // if neighbor or curr_entity is locally-used, add it to keeper
        if ( locally_used( bulk_data.bucket(neighbor)) || locally_used( bulk_data.bucket(curr_entity) ) ) {
          EntitySide keeper;
          keeper.inside.entity = curr_entity;
          keeper.inside.side_ordinal = subcell_identifier;
          keeper.outside = *adj_itr;

          boundary.push_back(keeper);
        }
      }
    }
  }
}


}
}
