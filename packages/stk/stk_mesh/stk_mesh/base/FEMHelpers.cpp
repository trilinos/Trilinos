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


#include <stk_mesh/base/FEMHelpers.hpp>
#include <sstream>                      // for operator<<, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Types.hpp>      // for PartVector, EntityId, etc
#include <string>                       // for char_traits, operator<<
#include "stk_mesh/base/CellTopology.hpp"  // for CellTopology
#include "stk_mesh/base/MetaData.hpp"   // for get_cell_topology, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_topology/topology.hpp"    // for topology::rank
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowErrorMsgIf, etc

namespace stk {
namespace mesh {

namespace {

void verify_declare_element_side(
    const BulkData & mesh,
    const Entity elem,
    const unsigned local_side_id
    )
{
  const CellTopologyData * const elem_top = get_cell_topology(mesh.bucket(elem)).getCellTopologyData();

  const CellTopologyData * const side_top =
    ( elem_top && local_side_id < elem_top->side_count )
    ? elem_top->side[ local_side_id ].topology : NULL ;

  ThrowErrorMsgIf( elem_top && local_side_id >= elem_top->side_count,
    "For elem " << mesh.identifier(elem) << ", local_side_id " << local_side_id << ", " <<
    "local_side_id exceeds " << elem_top->name << ".side_count = " << elem_top->side_count );

  ThrowErrorMsgIf( side_top == NULL,
    "For elem " << mesh.identifier(elem) << ", local_side_id " << local_side_id << ", " <<
    "No element topology found");
}

void verify_declare_element_edge(
    const BulkData & mesh,
    const Entity elem,
    const unsigned local_edge_id
    )
{
  const CellTopologyData * const elem_top = get_cell_topology( mesh.bucket(elem) ).getCellTopologyData();

  const CellTopologyData * const edge_top =
    ( elem_top && local_edge_id < elem_top->edge_count )
    ? elem_top->edge[ local_edge_id ].topology : NULL ;

  ThrowErrorMsgIf( elem_top && local_edge_id >= elem_top->edge_count,
    "For elem " << mesh.identifier(elem) << ", local_edge_id " << local_edge_id << ", " <<
    "local_edge_id exceeds " << elem_top->name << ".edge_count = " << elem_top->edge_count );

  ThrowErrorMsgIf( edge_top == NULL,
    "For elem " << mesh.identifier(elem) << ", local_edge_id " << local_edge_id << ", " <<
    "No element topology found");
}

} // unnamed namespace

Entity declare_element( BulkData & mesh ,
                        PartVector & parts ,
                        const EntityId elem_id ,
                        const EntityId node_id[] )
{
  MetaData & fem_meta = MetaData::get(mesh);
  const CellTopologyData * const top = fem_meta.get_cell_topology( *parts[0] ).getCellTopologyData();

  ThrowErrorMsgIf(top == NULL,
                  "Part " << parts[0]->name() << " does not have a local topology");

  PartVector empty ;

  const EntityRank entity_rank = stk::topology::ELEMENT_RANK;

  Entity elem = mesh.declare_entity( entity_rank, elem_id, parts );

  const EntityRank node_rank = stk::topology::NODE_RANK;

  Permutation perm = static_cast<Permutation>(0);
  OrdinalVector ordinal_scratch;
  ordinal_scratch.reserve(64);
  PartVector part_scratch;
  part_scratch.reserve(64);

  for ( unsigned i = 0 ; i < top->node_count ; ++i ) {
    //declare node if it doesn't already exist
    Entity node = mesh.get_entity( node_rank , node_id[i]);
    if ( !mesh.is_valid(node) ) {
      node = mesh.declare_entity( node_rank , node_id[i], empty );
    }

    mesh.declare_relation( elem , node , i, perm, ordinal_scratch, part_scratch );
  }
  return elem ;
}

Entity declare_element_side(
  BulkData & mesh,
  Entity elem ,
  Entity side,
  const unsigned local_side_id ,
  Part * part )
{
  verify_declare_element_side(mesh, elem, local_side_id);

  const CellTopologyData * const elem_top = get_cell_topology( mesh.bucket(elem) ).getCellTopologyData();

  ThrowErrorMsgIf( elem_top == NULL,
      "Element[" << mesh.identifier(elem) << "] has no defined topology" );

  const CellTopologyData * const side_top = elem_top->side[ local_side_id ].topology;

  ThrowErrorMsgIf( side_top == NULL,
      "Element[" << mesh.identifier(elem) << "], local_side_id = " <<
      local_side_id << ", side has no defined topology" );

  const unsigned * const side_node_map = elem_top->side[ local_side_id ].node ;

  PartVector add_parts ;
  Permutation perm = static_cast<Permutation>(0);
  OrdinalVector ordinal_scratch;
  ordinal_scratch.reserve(64);
  PartVector part_scratch;
  part_scratch.reserve(64);

  if ( part ) { add_parts.push_back( part ); }

  mesh.change_entity_parts(side, add_parts);

  mesh.declare_relation( elem , side , local_side_id, perm, ordinal_scratch, part_scratch );

  Entity const *elem_nodes = mesh.begin_nodes(elem);

  for ( unsigned i = 0 ; i < side_top->node_count ; ++i )
  {
    Entity node = elem_nodes[ side_node_map[i] ];
    mesh.declare_relation( side , node , i, perm, ordinal_scratch, part_scratch );
  }

  return side ;
}

Entity declare_element_edge(
  BulkData & mesh ,
  Entity elem ,
  Entity edge,
  const unsigned local_edge_id ,
  Part * part )
{
  const CellTopologyData * const elem_top = get_cell_topology(mesh.bucket(elem)).getCellTopologyData();

  ThrowErrorMsgIf( elem_top == NULL,
      "Element[" << mesh.identifier(elem) << "] has no defined topology" );

  const CellTopologyData * const edge_top = elem_top->edge[ local_edge_id ].topology;

  ThrowErrorMsgIf( edge_top == NULL,
      "Element[" << mesh.identifier(elem) << "], local_edge_id = " <<
      local_edge_id << ", edge has no defined topology" );

  const unsigned * const edge_node_map = elem_top->edge[ local_edge_id ].node ;

  PartVector add_parts ;
  Permutation perm = static_cast<Permutation>(0);
  OrdinalVector ordinal_scratch;
  ordinal_scratch.reserve(64);
  PartVector part_scratch;
  part_scratch.reserve(64);

  if ( part ) { add_parts.push_back( part ); }

  mesh.change_entity_parts(edge, add_parts);

  mesh.declare_relation( elem , edge , local_edge_id, perm, ordinal_scratch, part_scratch );

  Entity const *elem_nodes = mesh.begin_nodes(elem);

  for ( unsigned i = 0 ; i < edge_top->node_count ; ++i )
  {
    Entity node = elem_nodes[ edge_node_map[i] ];
    mesh.declare_relation( edge , node , i, perm, ordinal_scratch, part_scratch );
  }

  return edge ;
}

Entity declare_element_side(
  BulkData & mesh ,
  const stk::mesh::EntityId global_side_id ,
  Entity elem ,
  const unsigned local_side_id ,
  Part * part ,
  bool check_pre_existing )
{
  verify_declare_element_side(mesh, elem, local_side_id);

  const CellTopologyData * const elem_top = get_cell_topology(mesh.bucket(elem)).getCellTopologyData();

  ThrowErrorMsgIf( elem_top == NULL,
      "Element[" << mesh.identifier(elem) << "] has no defined topology");

  const CellTopologyData * const shards_side_top = elem_top->side[ local_side_id ].topology;
  const stk::topology side_top = get_topology(shards_side_top);

  ThrowErrorMsgIf( shards_side_top == NULL,
		   "Element[" << mesh.identifier(elem) << "], local_side_id = " <<
		   local_side_id << ", side has no defined topology" );

  PartVector empty_parts ;
  Entity side;
  if (check_pre_existing) {
    side = mesh.get_entity( side_top.rank(), global_side_id);
    if (!mesh.is_valid(side)) {
      side = mesh.declare_entity( side_top.rank() , global_side_id, empty_parts );
      declare_element_side(mesh, elem, side, local_side_id, part);
    }
  }
  else {
    side = mesh.declare_entity( side_top.rank() , global_side_id, empty_parts );
    declare_element_side(mesh, elem, side, local_side_id, part);
  }
  return side;
}

Entity declare_element_edge(
  BulkData & mesh ,
  const stk::mesh::EntityId global_edge_id ,
  Entity elem ,
  const unsigned local_edge_id ,
  Part * part,
  bool check_pre_existing )
{
  verify_declare_element_edge(mesh, elem, local_edge_id);

  const CellTopologyData * const elem_top = get_cell_topology(mesh.bucket(elem) ).getCellTopologyData();

  ThrowErrorMsgIf( elem_top == NULL,
      "Element[" << mesh.identifier(elem) << "] has no defined topology");


  const CellTopologyData * const shards_edge_top = elem_top->edge[ local_edge_id ].topology;
  const stk::topology edge_top = get_topology(shards_edge_top);

  ThrowErrorMsgIf( shards_edge_top == NULL,
      "Element[" << mesh.identifier(elem) << "], local_edge_id = " <<
      local_edge_id << ", edge has no defined topology" );

  PartVector empty_parts ;
  Entity edge;
  if (check_pre_existing) {
    edge = mesh.get_entity(edge_top.rank(), global_edge_id);
    if (!mesh.is_valid(edge)) {
      edge = mesh.declare_entity( edge_top.rank() , global_edge_id, empty_parts );
      declare_element_edge(mesh, elem, edge, local_edge_id, part);
    }
  }
  else {
    edge = mesh.declare_entity( edge_top.rank() , global_edge_id, empty_parts );
    declare_element_edge(mesh, elem, edge, local_edge_id, part);
  }
  return edge;
}

const CellTopologyData * get_subcell_nodes(const BulkData& mesh, const Entity entity ,
                                           EntityRank subcell_rank ,
                                           unsigned subcell_identifier ,
                                           EntityVector & subcell_nodes)
{
  ThrowAssert(subcell_rank <= stk::topology::ELEMENT_RANK);

  subcell_nodes.clear();

  // get cell topology
  const CellTopologyData* celltopology = get_cell_topology(mesh.bucket(entity)).getCellTopologyData();

  //error checking
  {
    //no celltopology defined
    if (celltopology == NULL) {
      return NULL;
    }

    // valid ranks fall within the dimension of the cell topology
    const bool bad_rank = subcell_rank >= celltopology->dimension;
    ThrowInvalidArgMsgIf( bad_rank, "subcell_rank is >= celltopology dimension\n");

    // subcell_identifier must be less than the subcell count
    const bool bad_id = subcell_identifier >= celltopology->subcell_count[subcell_rank];
    ThrowInvalidArgMsgIf( bad_id,   "subcell_id is >= subcell_count\n");
  }

  // Get the cell topology of the subcell
  const CellTopologyData * subcell_topology =
    celltopology->subcell[subcell_rank][subcell_identifier].topology;

  const int num_nodes_in_subcell = subcell_topology->node_count;

  // For the subcell, get it's local nodes ids
  const unsigned* subcell_node_local_ids =
    celltopology->subcell[subcell_rank][subcell_identifier].node;

  Entity const *node_relations = mesh.begin_nodes(entity);
  subcell_nodes.reserve(num_nodes_in_subcell);

  for (int i = 0; i < num_nodes_in_subcell; ++i )
  {
    subcell_nodes.push_back( node_relations[subcell_node_local_ids[i]] );
  }

  return subcell_topology;
}


int get_entity_subcell_id( const BulkData& mesh,
                           const Entity entity ,
                           const EntityRank subcell_rank,
                           const CellTopologyData & subcell_topology,
                           const std::vector<Entity>& subcell_nodes )
{
  ThrowAssert(subcell_rank <= stk::topology::ELEMENT_RANK);

  const int INVALID_SIDE = -1;

  unsigned num_nodes = subcell_topology.node_count;

  if (num_nodes != subcell_nodes.size()) {
    return INVALID_SIDE;
  }

  // get topology of elem
  const CellTopologyData* entity_topology = get_cell_topology(mesh.bucket(entity)).getCellTopologyData();
  if (entity_topology == NULL) {
    return INVALID_SIDE;
  }

  // get nodal relations for entity
  Entity const *node_rels = mesh.begin_nodes(entity);
  const int num_permutations = subcell_topology.permutation_count;

  // Iterate over the subcells of entity...
  for (unsigned local_subcell_ordinal = 0;
      local_subcell_ordinal < entity_topology->subcell_count[subcell_rank];
      ++local_subcell_ordinal) {

    // get topological data for this subcell
    const CellTopologyData& curr_subcell_topology =
      *entity_topology->subcell[subcell_rank][local_subcell_ordinal].topology;

    // If topologies are not the same, there is no way the subcells are the same
    if (&subcell_topology == &curr_subcell_topology) {

      const unsigned* const subcell_node_map = entity_topology->subcell[subcell_rank][local_subcell_ordinal].node;

      // Taking all positive permutations into account, check if this subcell
      // has the same nodes as the subcell_nodes argument. Note that this
      // implementation preserves the node-order so that we can take
      // entity-orientation into account.
      for (int p = 0; p < num_permutations; ++p) {

        if (curr_subcell_topology.permutation[p].polarity ==
            CELL_PERMUTATION_POLARITY_POSITIVE) {

          const unsigned * const perm_node =
            curr_subcell_topology.permutation[p].node ;

          bool all_match = true;
          for (unsigned j = 0 ; j < num_nodes; ++j ) {
            if (subcell_nodes[j] != node_rels[subcell_node_map[perm_node[j]]]) {
              all_match = false;
              break;
            }
          }

          // all nodes were the same, we have a match
          if ( all_match ) {
            return local_subcell_ordinal ;
          }
        }
      }
    }
  }

  return INVALID_SIDE;
}


}
}
