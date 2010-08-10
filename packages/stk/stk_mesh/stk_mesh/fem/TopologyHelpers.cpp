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
#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>

#include <stk_util/util/StaticAssert.hpp>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

const CellTopologyData * get_cell_topology( const Part & p )
{ return p.attribute<CellTopologyData>(); }

void set_cell_topology( Part & p , const CellTopologyData * singleton )
{
  static const char method[] = "stk::mesh::set_cell_topology" ;

  MetaData & m = p.mesh_meta_data();

  const CellTopologyData * t = NULL ;

  if ( p.mesh_meta_data().entity_rank_count() <= p.primary_entity_rank() ||
       singleton == NULL ||
       singleton != ( t = m.declare_attribute_no_delete(p,singleton) ) ) {
    std::ostringstream msg ;
    msg << method << "( " << p.name();
    msg << " entity_rank(" << p.primary_entity_rank() << ") , " ;
    if ( singleton ) { msg << singleton->name ; }
    else             { msg << "NULL" ; }
    msg << " ) ERROR" ;
    if ( t ) { msg << "Existing topology = " << t->name ; }
    throw std::runtime_error( msg.str() );
  }
}

//----------------------------------------------------------------------


const CellTopologyData * get_cell_topology( const Bucket & bucket )
{
  const CellTopologyData * top = NULL ;
  PartVector parts ;
  bucket.supersets( parts );

  PartVector::iterator i = parts.begin() ;

  for ( ; NULL == top && i != parts.end() ; ++i ) {
    if ( bucket.entity_rank() == (**i).primary_entity_rank() ) {
      top = get_cell_topology( **i );
    }
  }

  bool ok = true ;

  for ( ; ok && i != parts.end() ; ++i ) {
    if ( bucket.entity_rank() == (**i).primary_entity_rank() ) {
      const CellTopologyData * const tmp = get_cell_topology( **i );
      ok = ((tmp == NULL) || (tmp == top)) ;
    }
  }

  if ( ! ok ) {
    std::ostringstream msg ;
    msg << "stk::mesh::get_cell_topology( Bucket[" ;
    for ( i = parts.begin() ; i != parts.end() ; ++i ) {
      if ( bucket.entity_rank() == (**i).primary_entity_rank() ) {
        const CellTopologyData * const tmp = get_cell_topology( **i );
        msg << " " << (*i)->name();
        if ( tmp ) { msg << "->" << tmp->name ; }
        msg << " ] ) FAILED WITH MULTIPLE LOCAL TOPOLOGIES" ;
        throw std::runtime_error( msg.str() );
      }
    }
  }
  return top ;
}

const CellTopologyData * get_cell_topology( const Entity & entity )
{ return get_cell_topology( entity.bucket() ); }

//----------------------------------------------------------------------

namespace {

void verify_declare_element_side(
    const BulkData & mesh,
    const Entity & elem,
    const unsigned local_side_id
    )
{
  static const char method[] = "stk::mesh::declare_element_side" ;
  const CellTopologyData * const elem_top = get_cell_topology( elem );

  const CellTopologyData * const side_top =
    ( elem_top && local_side_id < elem_top->side_count )
    ? elem_top->side[ local_side_id ].topology : NULL ;

  bool different_bulk_data =  &mesh != & (elem.bucket().mesh());

  bool no_elem_top = NULL;
  if(side_top == NULL){
    no_elem_top = 1;
  }else {
    no_elem_top = 0;
  }

  bool bad_side_id = NULL;
  if( elem_top && local_side_id >= elem_top->side_count){
      bad_side_id = 1;
  }


  if ( different_bulk_data || no_elem_top || bad_side_id ) {
     std::ostringstream msg ;
     msg << method << "( ";
     print_entity_key( msg , mesh.mesh_meta_data() , elem.key() );
     msg << " , "
         << local_side_id
         << " ) FAILED" ;
     if ( different_bulk_data ) {
       msg << " Bulkdata for 'elem' and 'side' are different" ;
     } else if (no_elem_top && !bad_side_id ){
       msg << " No element topology found" ;

     } else if ( no_elem_top && bad_side_id ){
       msg << " No element topology found and " ;
       msg << " cell side id ";
       msg << local_side_id ;
       msg << " exceeds " ;
       msg << elem_top->name ;
       msg << ".side_count = " ;
       msg << elem_top->side_count ;
     }
     throw std::runtime_error( msg.str() );
   }
}

}

//----------------------------------------------------------------------

Entity & declare_element_side(
  Entity & elem ,
  Entity & side,
  const unsigned local_side_id ,
  Part * part )
{
  static const char method[] = "stk::mesh::declare_element_side" ;

  BulkData & mesh = side.bucket().mesh();

  verify_declare_element_side(mesh, elem, local_side_id);

  const CellTopologyData * const elem_top = get_cell_topology( elem );

  if (elem_top == NULL) {
    std::ostringstream msg ;
    msg << method ;
    msg << "( Element[" << elem.identifier() << "] has no defined topology" ;
    throw std::runtime_error( msg.str() );
  }

  const CellTopologyData * const side_top = elem_top->side[ local_side_id ].topology;

  if (side_top == NULL) {
    std::ostringstream msg ;
    msg << method ;
    msg << "( Element[" << elem.identifier() << "]" ;
    msg << " , local_side_id = " << local_side_id << " ) FAILED: " ;
    msg << " Side has no defined topology" ;
    throw std::runtime_error( msg.str() );
  }

  const unsigned * const side_node_map = elem_top->side[ local_side_id ].node ;

  // This is dangerous if the unsigned enums are changed. Try to catch at compile time...
  enum { DimensionMappingAssumption_OK =
           StaticAssert< stk::mesh::Edge == 1 && stk::mesh::Face == 2 >::OK };

  PartVector add_parts ;

  if ( part ) { add_parts.push_back( part ); }

  //\TODO refactor: is 'dimension' the right thing to use for EntityRank here???
  mesh.change_entity_parts(side, add_parts);

  mesh.declare_relation( elem , side , local_side_id );

  PairIterRelation rel = elem.relations( Node );

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
  static const char method[] = "stk::mesh::declare_element_side" ;
  verify_declare_element_side(mesh, elem, local_side_id);

  const CellTopologyData * const elem_top = get_cell_topology( elem );

  if (elem_top == NULL) {
    std::ostringstream msg ;
    msg << method ;
    msg << "( Element[" << elem.identifier() << "] has no defined topology" ;
    throw std::runtime_error( msg.str() );
  }

  const CellTopologyData * const side_top = elem_top->side[ local_side_id ].topology;

  if (side_top == NULL) {
    std::ostringstream msg ;
    msg << method ;
    msg << "( Element[" << elem.identifier() << "]" ;
    msg << " , local_side_id = " << local_side_id << " ) FAILED: " ;
    msg << " Side has no defined topology" ;
    throw std::runtime_error( msg.str() );
  }


  PartVector empty_parts ;

  Entity & side = mesh.declare_entity( side_top->dimension , global_side_id, empty_parts );
  return declare_element_side( elem, side, local_side_id, part);
}

//----------------------------------------------------------------------

bool element_side_polarity( const Entity & elem ,
                            const Entity & side , int local_side_id )
{
  static const char method[] = "stk::mesh::element_side_polarity" ;

  const bool is_side = side.entity_rank() != Edge ;

  const CellTopologyData * const elem_top = get_cell_topology( elem );

  const unsigned side_count = ! elem_top ? 0 : (
                                is_side ? elem_top->side_count
                                        : elem_top->edge_count );

  if ( NULL == elem_top ||
       local_side_id < 0 ||
       static_cast<int>(side_count) <= local_side_id ) {
    const MetaData & meta_data = elem.bucket().mesh().mesh_meta_data();
    std::ostringstream msg ;
    msg << method ;
    msg << "( Element[" << elem.identifier() << "]" ;
    msg << " , " << meta_data.entity_rank_names()[ side.entity_rank() ];
    msg << "[" << side.identifier() << "]" ;
    msg << " , local_side_id = " << local_side_id << " ) FAILED: " ;
    if ( NULL == elem_top ) {
      msg << " Element has no defined topology" ;
    }
    else {
      msg << " Unsupported local_side_id" ;
    }
    throw std::runtime_error( msg.str() );
  }

  const CellTopologyData * const side_top =
    is_side ? elem_top->side[ local_side_id ].topology
            : elem_top->edge[ local_side_id ].topology ;

  const unsigned * const side_map =
    is_side ? elem_top->side[ local_side_id ].node
            : elem_top->edge[ local_side_id ].node ;

  const PairIterRelation elem_nodes = elem.relations( Node );
  const PairIterRelation side_nodes = side.relations( Node );

  bool good = true ;
  for ( unsigned j = 0 ; good && j < side_top->node_count ; ++j ) {
    good = side_nodes[j].entity() == elem_nodes[ side_map[j] ].entity();
  }
  return good ;
}


int element_local_side_id( const Entity & elem ,
                           const CellTopologyData * side_topology,
                           const std::vector<Entity*>& side_nodes )
{

  // get topology of elem
  const CellTopologyData* elem_topology = get_cell_topology(elem);
  if (elem_topology == NULL) {
    return -1;
  }

  // get nodal relations for elem
  PairIterRelation relations = elem.relations(Node);

  const unsigned subcell_rank = elem_topology->dimension - 1;

  const int num_permutations = side_topology->permutation_count;

  // Iterate over the subcells of elem...
  for (unsigned itr = 0;
       itr < elem_topology->subcell_count[subcell_rank];
       ++itr) {

    // get topological data for this side
    const CellTopologyData* curr_side_topology =
      elem_topology->subcell[subcell_rank][itr].topology;
    unsigned num_nodes =
      elem_topology->subcell[subcell_rank][itr].topology->node_count;
    const unsigned* const side_node_map = elem_topology->side[itr].node;

    // If topologies are not the same, there is no way the sides are the same
    if (side_topology == curr_side_topology) {

      // Taking all positive permutations into account, check if this side
      // has the same nodes as the side_nodes argument
      for (int p = 0; p < num_permutations; ++p) {

        if (curr_side_topology->permutation[p].polarity ==
            CELL_PERMUTATION_POLARITY_POSITIVE) {

          const unsigned * const perm_node =
            curr_side_topology->permutation[p].node ;

          bool all_match = true;
          for (unsigned j = 0 ; j < num_nodes; ++j ) {
            if (side_nodes[j] !=
                relations[side_node_map[perm_node[j]]].entity()) {
              all_match = false;
              break;
            }
          }

          // all nodes were the same, we have a match
          if ( all_match ) {
            return itr ;
          }
        }
      }
    }
  }

  return -1;
}

//----------------------------------------------------------------------

}// namespace mesh
}// namespace stk

