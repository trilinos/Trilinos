/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards
 */

#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <cassert>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>

#include <stk_util/util/StaticAssert.hpp>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

namespace {

void verify_declare_element_side(
    const BulkData & mesh,
    const Entity & elem,
    const unsigned local_side_id
    )
{
#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const CellTopologyData * const elem_top = get_cell_topology( elem );
#else // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const CellTopologyData * const elem_top = fem::get_cell_topology( elem ).getCellTopologyData();
#endif // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS

  const CellTopologyData * const side_top =
    ( elem_top && local_side_id < elem_top->side_count )
    ? elem_top->side[ local_side_id ].topology : NULL ;

  ThrowErrorMsgIf( &mesh != & (elem.bucket().mesh()),
    "For elem " << print_entity_key(elem) <<
    ", Bulkdata for 'elem' and mesh are different");

  ThrowErrorMsgIf( elem_top && local_side_id >= elem_top->side_count,
    "For elem " << print_entity_key(elem) << ", local_side_id " << local_side_id << ", " <<
    "local_side_id exceeds " << elem_top->name << ".side_count = " << elem_top->side_count );

  ThrowErrorMsgIf( side_top == NULL,
    "For elem " << print_entity_key(elem) << ", local_side_id " << local_side_id << ", " <<
    "No element topology found");
}

}

//----------------------------------------------------------------------

Entity & declare_element_side(
  Entity & elem ,
  Entity & side,
  const unsigned local_side_id ,
  Part * part )
{
  BulkData & mesh = side.bucket().mesh();

  verify_declare_element_side(mesh, elem, local_side_id);

#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const CellTopologyData * const elem_top = get_cell_topology( elem );
#else // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const CellTopologyData * const elem_top = fem::get_cell_topology( elem ).getCellTopologyData();
#endif // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS

  ThrowErrorMsgIf( elem_top == NULL,
      "Element[" << elem.identifier() << "] has no defined topology" );

  const CellTopologyData * const side_top = elem_top->side[ local_side_id ].topology;

  ThrowErrorMsgIf( side_top == NULL,
      "Element[" << elem.identifier() << "], local_side_id = " <<
      local_side_id << ", side has no defined topology" );

  const unsigned * const side_node_map = elem_top->side[ local_side_id ].node ;

  PartVector add_parts ;

  if ( part ) { add_parts.push_back( part ); }

  mesh.change_entity_parts(side, add_parts);

  mesh.declare_relation( elem , side , local_side_id );

  PairIterRelation rel = elem.relations( fem::NODE_RANK );

  for ( unsigned i = 0 ; i < side_top->node_count ; ++i ) {
    Entity & node = * rel[ side_node_map[i] ].entity();
    mesh.declare_relation( side , node , i );
  }

  return side ;
}

//----------------------------------------------------------------------

Entity & declare_element_side(
  BulkData & mesh ,
  const stk::mesh::EntityId global_side_id ,
  Entity & elem ,
  const unsigned local_side_id ,
  Part * part )
{
  verify_declare_element_side(mesh, elem, local_side_id);

#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const CellTopologyData * const elem_top = get_cell_topology( elem );
#else // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const CellTopologyData * const elem_top = fem::get_cell_topology( elem ).getCellTopologyData();
#endif // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS

  ThrowErrorMsgIf( elem_top == NULL,
      "Element[" << elem.identifier() << "] has no defined topology");

  const CellTopologyData * const side_top = elem_top->side[ local_side_id ].topology;

  ThrowErrorMsgIf( side_top == NULL,
      "Element[" << elem.identifier() << "], local_side_id = " <<
      local_side_id << ", side has no defined topology" );

  PartVector empty_parts ;

  Entity & side = mesh.declare_entity( side_top->dimension , global_side_id, empty_parts );
  return declare_element_side( elem, side, local_side_id, part);
}

//----------------------------------------------------------------------

bool element_side_polarity( const Entity & elem ,
                            const Entity & side , int local_side_id )
{
  // 09/14/10:  TODO:  tscoffe:  Will this work in 1D?
  // 09/14/10:  TODO:  tscoffe:  We need an exception here if we don't get a TopologicalMetaData off of MetaData or we need to take one on input.
#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const bool is_side = side.entity_rank() != Edge;
  const CellTopologyData * const elem_top = get_cell_topology( elem );
#else // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  stk::mesh::fem::FEMInterface &fem = stk::mesh::fem::get_fem_interface(elem);
  const bool is_side = side.entity_rank() != fem::edge_rank(fem);
  const CellTopologyData * const elem_top = fem::get_cell_topology( elem ).getCellTopologyData();
#endif // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS

  const unsigned side_count = ! elem_top ? 0 : (
                                is_side ? elem_top->side_count
                                        : elem_top->edge_count );

  ThrowErrorMsgIf( elem_top == NULL,
                   "For Element[" << elem.identifier() << "], element has no defined topology");

  ThrowErrorMsgIf( local_side_id < 0 || static_cast<int>(side_count) <= local_side_id,
    "For Element[" << elem.identifier() << "], " <<
    "side: " << print_entity_key(side) << ", " <<
    "local_side_id = " << local_side_id <<
    " ; unsupported local_side_id");

  const CellTopologyData * const side_top =
    is_side ? elem_top->side[ local_side_id ].topology
            : elem_top->edge[ local_side_id ].topology ;

  const unsigned * const side_map =
    is_side ? elem_top->side[ local_side_id ].node
            : elem_top->edge[ local_side_id ].node ;

  const PairIterRelation elem_nodes = elem.relations( fem::NODE_RANK );
  const PairIterRelation side_nodes = side.relations( fem::NODE_RANK );

  bool good = true ;
  for ( unsigned j = 0 ; good && j < side_top->node_count ; ++j ) {
    good = side_nodes[j].entity() == elem_nodes[ side_map[j] ].entity();
  }
  return good ;
}

int get_entity_subcell_id( const Entity & entity ,
                           const EntityRank subcell_rank,
                           const CellTopologyData * subcell_topology,
                           const std::vector<Entity*>& subcell_nodes )
{
  const int INVALID_SIDE = -1;

  unsigned num_nodes = subcell_topology->node_count;

  if (num_nodes != subcell_nodes.size()) {
    return INVALID_SIDE;
  }

  // get topology of elem
#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const CellTopologyData* entity_topology = get_cell_topology(entity);
#else // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const CellTopologyData* entity_topology = fem::get_cell_topology(entity).getCellTopologyData();
#endif // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  if (entity_topology == NULL) {
    return INVALID_SIDE;
  }

  // get nodal relations for entity
  PairIterRelation relations = entity.relations(fem::NODE_RANK);

  const int num_permutations = subcell_topology->permutation_count;

  // Iterate over the subcells of entity...
  for (unsigned local_subcell_ordinal = 0;
      local_subcell_ordinal < entity_topology->subcell_count[subcell_rank];
      ++local_subcell_ordinal) {

    // get topological data for this subcell
    const CellTopologyData* curr_subcell_topology =
      entity_topology->subcell[subcell_rank][local_subcell_ordinal].topology;

    // If topologies are not the same, there is no way the subcells are the same
    if (subcell_topology == curr_subcell_topology) {

      const unsigned* const subcell_node_map = entity_topology->subcell[subcell_rank][local_subcell_ordinal].node;

      // Taking all positive permutations into account, check if this subcell
      // has the same nodes as the subcell_nodes argument. Note that this
      // implementation preserves the node-order so that we can take
      // entity-orientation into account.
      for (int p = 0; p < num_permutations; ++p) {

        if (curr_subcell_topology->permutation[p].polarity ==
            CELL_PERMUTATION_POLARITY_POSITIVE) {

          const unsigned * const perm_node =
            curr_subcell_topology->permutation[p].node ;

          bool all_match = true;
          for (unsigned j = 0 ; j < num_nodes; ++j ) {
            if (subcell_nodes[j] !=
                relations[subcell_node_map[perm_node[j]]].entity()) {
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

//----------------------------------------------------------------------

const CellTopologyData * get_elem_side_nodes( const Entity & elem,
                          RelationIdentifier side_ordinal,
                          EntityVector & side_nodes
                        )
{
  side_nodes.clear();

  const CellTopologyData * const elem_top = fem::get_cell_topology( elem ).getCellTopologyData();
  if (elem_top == NULL) {
    return NULL;
  }

  const CellTopologyData * const side_top = elem_top->side[ side_ordinal ].topology;
  if (side_top == NULL) {
    return NULL;
  }

  const unsigned * const side_node_map = elem_top->side[ side_ordinal ].node ;

  PairIterRelation relations = elem.relations( fem::NODE_RANK );

  // Find positive polarity permutation that starts with lowest entity id:
  const int num_permutations = side_top->permutation_count;
  int lowest_entity_id_permutation = 0;
  for (int p = 0; p < num_permutations; ++p) {

    if (side_top->permutation[p].polarity ==
        CELL_PERMUTATION_POLARITY_POSITIVE) {

      const unsigned * const pot_lowest_perm_node =
        side_top->permutation[p].node ;

      const unsigned * const curr_lowest_perm_node =
        side_top->permutation[lowest_entity_id_permutation].node;

      unsigned first_node_index = 0;

      unsigned current_lowest_entity_id =
        relations[side_node_map[curr_lowest_perm_node[first_node_index]]].entity()->identifier();

      unsigned potential_lowest_entity_id =
        relations[side_node_map[pot_lowest_perm_node[first_node_index]]].entity()->identifier();

      if ( potential_lowest_entity_id < current_lowest_entity_id ) {
        lowest_entity_id_permutation = p;
      }
    }
  }

  const unsigned * const perm_node =
    side_top->permutation[lowest_entity_id_permutation].node;

  side_nodes.reserve(side_top->node_count);
  for ( unsigned i = 0 ; i < side_top->node_count ; ++i ) {
    Entity * node = relations[side_node_map[perm_node[i]]].entity();
    side_nodes.push_back(node);
  }

  return side_top;

}
}// namespace mesh
}// namespace stk

