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
#include <stk_mesh/fem/EntityTypes.hpp>
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

  if ( p.mesh_meta_data().entity_type_count() <= p.primary_entity_type() ||
       singleton == NULL ||
       singleton != ( t = m.declare_attribute_no_delete(p,singleton) ) ) {
    std::ostringstream msg ;
    msg << method << "( " << p.name();
    msg << " entity_type(" << p.primary_entity_type() << ") , " ;
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
    if ( bucket.entity_type() == (**i).primary_entity_type() ) {
      top = get_cell_topology( **i );
    }
  }

  bool ok = true ;
  
  for ( ; ok && i != parts.end() ; ++i ) {
    if ( bucket.entity_type() == (**i).primary_entity_type() ) {
      const CellTopologyData * const tmp = get_cell_topology( **i );
      ok = ((tmp == NULL) || (tmp == top)) ;
    }
  }
  
  if ( ! ok ) {
    std::ostringstream msg ;
    msg << "stk::mesh::get_cell_topology( Bucket[" ;
    for ( i = parts.begin() ; i != parts.end() ; ++i ) {
      if ( bucket.entity_type() == (**i).primary_entity_type() ) {
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

Entity & declare_element_side(
  BulkData & mesh ,
  const stk::mesh::EntityId global_side_id ,
  Entity & elem ,
  const unsigned local_side_id ,
  Part * part )
{
  static const char method[] = "stk::mesh::declare_element_side" ;

  const CellTopologyData * const elem_top = get_cell_topology( elem );

  const CellTopologyData * const side_top =
    ( elem_top && local_side_id < elem_top->side_count )
    ? elem_top->side[ local_side_id ].topology : NULL ;

  if ( NULL == side_top ) {
     std::ostringstream msg ;
     msg << method << "( mesh , "
         << global_side_id
         << " , " ;
     print_entity_key( msg , mesh.mesh_meta_data() , elem.key() );
     msg << " , "
         << local_side_id
         << " ) FAILED" ;
     if ( NULL == elem_top ) {
       msg << " Cannot discern element topology" ;
     }
     else {
       msg << " Cell side id exceeds " ;
       msg << elem_top->name ;
       msg << ".side_count = " ;
       msg << elem_top->side_count ;
     }
     throw std::runtime_error( msg.str() );
   }

  const unsigned * const side_node_map = elem_top->side[ local_side_id ].node ;

  // This is dangerous if the unsigned enums are changed. Try to catch at compile time...
  enum { DimensionMappingAssumption_OK =
           StaticAssert< stk::mesh::Edge == 1 && stk::mesh::Face == 2 >::OK };

  PartVector add_parts ;

  if ( part ) { add_parts.push_back( part ); }

  //\TODO refactor: is 'dimension' the right thing to use for EntityType here???
  Entity & side = mesh.declare_entity( side_top->dimension, global_side_id, add_parts );

  mesh.declare_relation( elem , side , local_side_id );

  PairIterRelation rel = elem.relations( Node );

  for ( unsigned i = 0 ; i < side_top->node_count ; ++i ) {
    Entity & node = * rel[ side_node_map[i] ].entity();
    mesh.declare_relation( side , node , i );
  }

  return side ;
}

//----------------------------------------------------------------------

bool element_side_polarity( const Entity & elem ,
                            const Entity & side , int local_side_id )
{
  static const char method[] = "stk::mesh::element_side_polarity" ;

  const bool is_side = side.entity_type() != Edge ;

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
    msg << " , " << meta_data.entity_type_names()[ side.entity_type() ];
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

void get_adjacent_entities( const Entity & entity ,
                            unsigned subcell_rank ,
                            unsigned subcell_identifier ,
                            std::vector<std::pair<Entity*, unsigned> > & adjacent_entities )
{
  adjacent_entities.clear();

  // get cell topology
  const CellTopologyData* celltopology = get_cell_topology(entity);
  if (celltopology == NULL) {
    return;
  }

  // valid ranks fall within the dimension of the cell topology
  bool bad_rank = subcell_rank >= celltopology->dimension;

  // local id should be < number of entities of the desired type
  // (if you have 4 edges, their ids should be 0-3)
  bool bad_id = false;
  if (!bad_rank) {
    bad_id = subcell_identifier >= celltopology->subcell_count[subcell_rank];
  }

  if (bad_rank || bad_id) {
    std::ostringstream msg;
    //parallel consisent throw
    if (bad_rank) {
      msg << "stk::mesh::get_adjacent_entities( const Entity& entity, unsigned subcell_rank, ... ) subcell_rank is >= celltopology dimension\n";
    }
    else if (bad_id) {
      msg << "stk::mesh::get_adjacent_entities( const Entity& entity, unsigned subcell_rank, unsigned subcell_identifier, ... ) subcell_identifier is >= subcell count\n";
    }

    throw std::runtime_error(msg.str());
  }

  // For the potentially common subcell, get it's nodes and num_nodes
  const unsigned* nodes = celltopology->subcell[subcell_rank][subcell_identifier].node;
  unsigned num_nodes = celltopology->subcell[subcell_rank][subcell_identifier].topology->node_count;

  // Get all the nodal relationships for this entity. We are guaranteed
  // that, if we make it this far, the entity is guaranteed to have
  // some relationship to nodes (we know it is a higher-order entity
  // than Node).
  PairIterRelation relations = entity.relations(Node);

  // Get the node entities that are related to entity
  std::vector<Entity*> node_entities;
  for (unsigned itr = 0; itr < num_nodes; ++itr) {
    node_entities.push_back(relations[nodes[itr]].entity());
  }

  // Given the nodes related to the original entity, find all entities
  // of similar rank that have some relation to one or more of these nodes
  std::vector<Entity*> elements;
  get_entities_through_relations(node_entities, entity.entity_type(), elements);

  // Make sure to remove the original entity from the list
  bool found = false;
  for (std::vector<Entity*>::iterator itr = elements.begin();
       itr != elements.end(); ++itr) {
    if (*itr == &entity) {
      elements.erase(itr);
      found = true;
      break;
    }
  }
  // The original entity should be related to the nodes of its subcells
  assert(found);

  // Add the local ids, from the POV of the adj entitiy, to the return value
  for (std::vector<Entity*>::const_iterator itr = elements.begin();
       itr != elements.end(); ++itr) {
    unsigned local_side_num = element_local_side_id(**itr, node_entities);
    adjacent_entities.push_back(std::pair<Entity*, unsigned>(*itr, local_side_num));
  }
}

int element_local_side_id( const Entity & elem ,
                           const Entity & side )
{
  return -1;
}

int element_local_side_id( const Entity & elem ,
                           const std::vector<Entity*>& entity_nodes )
{
  // sort the input nodes
  std::vector<Entity*> sorted_entity_nodes(entity_nodes);
  std::sort(sorted_entity_nodes.begin(), sorted_entity_nodes.end(), EntityLess());

  // get topology of elem
  const CellTopologyData* celltopology = get_cell_topology(elem);
  if (celltopology == NULL) {
    return -1;
  }

  // get nodal relations
  PairIterRelation relations = elem.relations(Node);

  const unsigned subcell_rank = celltopology->dimension - 1;

  // Iterate over the subcells of elem
  for (unsigned itr = 0; itr < celltopology->subcell_count[subcell_rank]; ++itr) {
    // get the nodes for this subcell
    const unsigned* nodes = celltopology->subcell[subcell_rank][itr].node;
    unsigned num_nodes = celltopology->subcell[subcell_rank][itr].topology->node_count;

    // Get the nodes in the subcell ???
    std::vector<Entity*> node_entities;
    for (unsigned nitr = 0; nitr < num_nodes; ++nitr) {
      node_entities.push_back(relations[nodes[nitr]].entity());
    }

    // check to see if this subcell exactly contains the nodes that were passed in
    std::sort(node_entities.begin(), node_entities.end(), EntityLess());
    if (node_entities == sorted_entity_nodes) {
      return itr;
    }
  }
  return -1;
}

//----------------------------------------------------------------------

}// namespace mesh
}// namespace stk

