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
#include "stk_mesh/base/MetaData.hpp"   // for get_topology, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_topology/topology.hpp"    // for topology, etc
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
  stk::topology elem_top = mesh.bucket(elem).topology();
  stk::topology invalid = stk::topology::INVALID_TOPOLOGY;

  stk::topology side_top =
    ( (elem_top!=stk::topology::INVALID_TOPOLOGY) && (local_side_id < elem_top.num_sides()) )
    ? elem_top.side_topology(local_side_id) : invalid ;

  ThrowErrorMsgIf( elem_top!=stk::topology::INVALID_TOPOLOGY && local_side_id >= elem_top.num_sides(),
    "For elem " << mesh.identifier(elem) << ", local_side_id " << local_side_id << ", " <<
    "local_side_id exceeds " << elem_top.name() << ".num_sies() = " << elem_top.num_sides() );

  ThrowErrorMsgIf( side_top == stk::topology::INVALID_TOPOLOGY,
    "For elem " << mesh.identifier(elem) << ", local_side_id " << local_side_id << ", " <<
    "No element topology found");
}

void verify_declare_element_edge(
    const BulkData & mesh,
    const Entity elem,
    const unsigned local_edge_id
    )
{
  stk::topology elem_top = mesh.bucket(elem).topology();
  stk::topology invalid = stk::topology::INVALID_TOPOLOGY;

  stk::topology edge_top =
    ( elem_top!=stk::topology::INVALID_TOPOLOGY && local_edge_id < elem_top.num_edges() )
    ? elem_top.edge_topology() : invalid;

  ThrowErrorMsgIf( elem_top!=stk::topology::INVALID_TOPOLOGY && local_edge_id >= elem_top.num_edges(),
    "For elem " << mesh.identifier(elem) << ", local_edge_id " << local_edge_id << ", " <<
    "local_edge_id exceeds " << elem_top.name() << ".edge_count = " << elem_top.num_edges() );

  ThrowErrorMsgIf( edge_top == stk::topology::INVALID_TOPOLOGY,
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
  stk::topology top = fem_meta.get_topology( *parts[0] );

  ThrowErrorMsgIf(top == stk::topology::INVALID_TOPOLOGY,
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

  for ( unsigned i = 0 ; i < top.num_nodes(); ++i ) {
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

  stk::topology elem_top = mesh.bucket(elem).topology();

  ThrowErrorMsgIf( elem_top == stk::topology::INVALID_TOPOLOGY,
      "Element[" << mesh.identifier(elem) << "] has no defined topology" );

  stk::topology side_top = elem_top.side_topology( local_side_id );

  ThrowErrorMsgIf( side_top == stk::topology::INVALID_TOPOLOGY,
      "Element[" << mesh.identifier(elem) << "], local_side_id = " <<
      local_side_id << ", side has no defined topology" );

  std::vector<unsigned> side_node_map(side_top.num_nodes());
  elem_top.side_node_ordinals(local_side_id, side_node_map.begin());

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

  for ( unsigned i = 0 ; i < side_top.num_nodes() ; ++i )
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
  stk::topology elem_top = mesh.bucket(elem).topology();

  ThrowErrorMsgIf( elem_top == stk::topology::INVALID_TOPOLOGY,
      "Element[" << mesh.identifier(elem) << "] has no defined topology" );

  stk::topology edge_top = elem_top.edge_topology();

  ThrowErrorMsgIf( edge_top == stk::topology::INVALID_TOPOLOGY,
      "Element[" << mesh.identifier(elem) << "], local_edge_id = " <<
      local_edge_id << ", edge has no defined topology" );

  std::vector<unsigned> edge_node_map(edge_top.num_nodes());
  elem_top.edge_node_ordinals( local_edge_id, edge_node_map.begin() );

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

  for ( unsigned i = 0 ; i < edge_top.num_nodes(); ++i )
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

  stk::topology elem_top = mesh.bucket(elem).topology();

  ThrowErrorMsgIf( elem_top == stk::topology::INVALID_TOPOLOGY,
      "Element[" << mesh.identifier(elem) << "] has no defined topology");

  stk::topology side_top = elem_top.side_topology( local_side_id );

  ThrowErrorMsgIf( side_top == stk::topology::INVALID_TOPOLOGY,
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

  stk::topology elem_top = mesh.bucket(elem).topology();

  ThrowErrorMsgIf( elem_top == stk::topology::INVALID_TOPOLOGY,
      "Element[" << mesh.identifier(elem) << "] has no defined topology");


  stk::topology edge_top = elem_top.edge_topology();

  ThrowErrorMsgIf( edge_top == stk::topology::INVALID_TOPOLOGY,
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

stk::topology get_subcell_nodes(const BulkData& mesh, const Entity entity ,
                                           EntityRank subcell_rank ,
                                           unsigned subcell_identifier ,
                                           EntityVector & subcell_nodes)
{
  ThrowAssert(subcell_rank <= stk::topology::ELEMENT_RANK);

  subcell_nodes.clear();

  // get cell topology
  stk::topology celltopology = mesh.bucket(entity).topology();

  //error checking
  {
    //no celltopology defined
    if (celltopology == stk::topology::INVALID_TOPOLOGY) {
      return celltopology;
    }

    // valid ranks fall within the dimension of the cell topology
    const bool bad_rank = subcell_rank >= celltopology.dimension();
    ThrowInvalidArgMsgIf( bad_rank, "subcell_rank is >= celltopology dimension\n");

    // subcell_identifier must be less than the subcell count
    const bool bad_id = subcell_identifier >= celltopology.num_sub_topology(subcell_rank);
    ThrowInvalidArgMsgIf( bad_id,   "subcell_id is >= subcell_count\n");
  }

  // Get the cell topology of the subcell
  stk::topology subcell_topology =
    celltopology.sub_topology(subcell_rank, subcell_identifier);

  const int num_nodes_in_subcell = subcell_topology.num_nodes();

  // For the subcell, get it's local nodes ids
  std::vector<unsigned> subcell_node_local_ids(num_nodes_in_subcell);
  celltopology.sub_topology_node_ordinals(subcell_rank, subcell_identifier, subcell_node_local_ids.begin());

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
                           stk::topology subcell_topology,
                           const std::vector<Entity>& subcell_nodes )
{
  ThrowAssert(subcell_rank <= stk::topology::ELEMENT_RANK);
  const int INVALID_SIDE = -1;

  stk::topology entity_topology = mesh.bucket(entity).topology();
  ThrowAssert(entity_topology.num_nodes() == mesh.num_nodes(entity));

  const Entity* entity_nodes = mesh.begin_nodes(entity);
  std::vector<Entity> topology_sub_nodes(subcell_nodes.size());

  for(size_t i=0; i<entity_topology.num_sub_topology(subcell_rank); ++i) {
      entity_topology.sub_topology_nodes(entity_nodes, subcell_rank, i, topology_sub_nodes.begin());
      if (subcell_topology.equivalent(topology_sub_nodes, subcell_nodes).first) {
          return i;
      }
  }

  return INVALID_SIDE;
}

}
}
