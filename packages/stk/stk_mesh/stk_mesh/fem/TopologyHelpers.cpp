/**
 * @author H. Carter Edwards
 */

#include <stdexcept>
#include <sstream>

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

//----------------------------------------------------------------------

}// namespace mesh
}// namespace stk

