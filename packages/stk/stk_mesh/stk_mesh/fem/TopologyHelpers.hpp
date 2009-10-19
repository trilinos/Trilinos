#ifndef stk_mesh_TopologyHelpers_hpp
#define stk_mesh_TopologyHelpers_hpp

#include <stdexcept>
#include <Shards_CellTopologyTraits.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/fem/EntityTypes.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_bulk_data_element
 *  \{
 */

/// \todo REFACTOR: The functions in this file represent a "bridge"
///between the mesh and the Shards_CellTopologyData stuff. Does it
///belong here?

//----------------------------------------------------------------------
/** \brief Attach a CellTopology to a Part.
 *  There is at most one cell topology allowed.
 */
void set_cell_topology( Part & , const CellTopologyData * singleton );

/** \brief  Attach a CellTopology to a Part.
 *  There is at most one element topology allowed.
 */
template< class Traits >
void set_cell_topology( Part & p )
{ return set_cell_topology( p , shards::getCellTopologyData<Traits>() ); }

/** \brief  The the CellTopology attached to the Part, if any */
const CellTopologyData * get_cell_topology( const Part & );

/** \brief  The the CellTopology attached to at most one Part of the Bucket */
const CellTopologyData * get_cell_topology( const Bucket & );

/** \brief  The the CellTopology attached to at most one Part of the Entity */
const CellTopologyData * get_cell_topology( const Entity & );

//----------------------------------------------------------------------
template< class Traits >
void get_parts_with_topology(stk::mesh::BulkData& mesh,
                             stk::mesh::PartVector& parts)
{
  parts.clear();

  const stk::mesh::PartVector& all_parts = mesh.mesh_meta_data().get_parts();

  stk::mesh::PartVector::const_iterator
    iter = all_parts.begin(),
    iter_end = all_parts.end();
 
  const CellTopologyData* topology = shards::getCellTopologyData<Traits>();

  for(; iter!=iter_end; ++iter) {
    stk::mesh::Part* part =  *iter;
    if (stk::mesh::get_cell_topology(*part) == topology) {
      parts.push_back(part);
    }
  }
}

//----------------------------------------------------------------------
/** \brief  Declare an element member of a Part with a CellTopology
 *          and nodes conformal to that topology.
 */
template< typename IdType >
inline
Entity & declare_element( BulkData & mesh ,
                          Part & part ,
                          const IdType elem_id ,
                          const IdType node_id[] )
{
  const CellTopologyData * const top = get_cell_topology( part );
     
  if ( top == NULL ) {
    std::ostringstream msg ; 
    msg << "stk::mesh::declare_element( mesh , " ; 
    msg << part.name();
    msg << " , " ;
    msg << elem_id ;
    msg << " , node_id[] ) ERROR, Part does not have a local topology" ;
    throw std::runtime_error( msg.str() ); 
  }
 
  PartVector empty ;
  PartVector add( 1 ); add[0] = & part ;
 
  Entity & elem = mesh.declare_entity( Element, elem_id, add ); 
 
  for ( unsigned i = 0 ; i < top->node_count ; ++i ) {
    Entity & node = mesh.declare_entity( Node, node_id[i], empty );
    mesh.declare_relation( elem , node , i );
  }   
  return elem ; 
}         

/** \brief  Declare an element member of a Part with a CellTopology
 *          and nodes conformal to that topology.
 */
template< typename IdType >
Entity & declare_element( BulkData & mesh ,
                          Part & part ,
                          const IdType elem_id ,
                          Entity * node[] )
{
  const CellTopologyData * const top = get_cell_topology( part );
     
  if ( top == NULL ) {
    std::ostringstream msg ; 
    msg << "stk::mesh::declare_element( mesh , " ; 
    msg << part.name();
    msg << " , " ;
    msg << elem_id ;
    msg << " , node[] ) ERROR, Part does not have a local topology" ;
    throw std::runtime_error( msg.str() ); 
  }
 
  PartVector add( 1 ); add[0] = & part ;
 
  Entity & elem = mesh.declare_entity( Element, elem_id, add ); 
 
  for ( unsigned i = 0 ; i < top->node_count ; ++i ) {
    mesh.declare_relation( elem , *node[i] , i );
  }   
  return elem ; 
}         

//----------------------------------------------------------------------
/** \brief  Create (or find) an element side.
 *
 *  The element must be a member of a Part with a CellTopology.
 */
Entity & declare_element_side( BulkData & mesh ,
                               const stk::mesh::EntityId global_side_id ,
                               Entity & elem ,
                               const unsigned local_side_id ,
                               Part * part = NULL );

/** \brief  Determine the polarity of the local side,
 *          more efficient if the local_side_id is known.
 */
bool element_side_polarity( const Entity & elem ,
                            const Entity & side , int local_side_id = -1 );

/** \brief  Determine the local side identifier,
 *          return -1 if the side doesn't match the element.
 */
int element_local_side_id( const Entity & elem ,
                           const Entity & side );

//----------------------------------------------------------------------

/** \} */

}//namespace mesh
}//namespace stk

#endif

