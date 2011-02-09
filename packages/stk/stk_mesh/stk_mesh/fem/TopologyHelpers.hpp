/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_TopologyHelpers_hpp
#define stk_mesh_TopologyHelpers_hpp

#include <sstream>
#include <stdexcept>
#include <Shards_CellTopologyTraits.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/fem/FEMInterface.hpp>

#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
#include <stk_mesh/fem/TopologyHelpersDeprecated.hpp>
#endif // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_bulk_data_element
 *  \{
 */

/// \todo REFACTOR: The functions in this file represent a "bridge"
///between the mesh and the Shards_CellTopologyData stuff. Does it
///belong here?

//----------------------------------------------------------------------
template< class Traits >
void get_parts_with_topology(stk::mesh::BulkData& mesh,
                             stk::mesh::PartVector& parts)
{
  parts.clear();

  const stk::mesh::PartVector& all_parts = MetaData::get(mesh).get_parts();

  stk::mesh::PartVector::const_iterator
    iter = all_parts.begin(),
    iter_end = all_parts.end();

  const CellTopologyData* topology = shards::getCellTopologyData<Traits>();

  for(; iter!=iter_end; ++iter) {
    stk::mesh::Part* part =  *iter;
#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
    if (get_cell_topology(*part) == topology) {
      parts.push_back(part);
    }
#else // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
    if (fem::get_cell_topology(*part).getCellTopologyData() == topology) {
      parts.push_back(part);
    }
#endif // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
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
#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const CellTopologyData * const top = get_cell_topology( part );
#else // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const CellTopologyData * const top = fem::get_cell_topology( part ).getCellTopologyData();
#endif // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS

  ThrowErrorMsgIf(top == NULL,
                  "Part " << part.name() << " does not have a local topology");

  PartVector empty ;
  PartVector add( 1 ); add[0] = & part ;

#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const EntityRank entity_rank = element_rank_deprecated(MetaData::get(part));
#else // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const EntityRank entity_rank = fem::get_entity_rank(MetaData::get(part), top);
#endif // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS

  Entity & elem = mesh.declare_entity( entity_rank, elem_id, add );

  for ( unsigned i = 0 ; i < top->node_count ; ++i ) {
    //declare node if it doesn't already exist
    Entity * node = mesh.get_entity( fem::NODE_RANK , node_id[i]);
    if ( NULL == node) {
      node = & mesh.declare_entity( fem::NODE_RANK , node_id[i], empty );
    }

    mesh.declare_relation( elem , *node , i );
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
#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const CellTopologyData * const top = get_cell_topology( part );
#else // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const CellTopologyData * const top = fem::get_cell_topology( part ).getCellTopologyData();
#endif // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS

  ThrowErrorMsgIf(top == NULL,
                  "Part " << part.name() << " does not have a local topology");

  PartVector add( 1 ); add[0] = & part ;

#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const EntityRank entity_rank = element_rank_deprecated(MetaData::get(part));
#else // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
  const EntityRank entity_rank = top->dimension;
#endif // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS

  Entity & elem = mesh.declare_entity( entity_rank, elem_id, add );

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

Entity & declare_element_side( Entity & elem ,
                               Entity & side ,
                               const unsigned local_side_id ,
                               Part * part = NULL );

/** \brief  Determine the polarity of the local side,
 *          more efficient if the local_side_id is known.
 */
bool element_side_polarity( const Entity & elem ,
                            const Entity & side , int local_side_id = -1 );


/** \brief  Given an entity and collection of nodes, return the
 *          local id of the subcell that contains those nodes in the
 *          correct orientation.
 */
int get_entity_subcell_id( const Entity            & entity ,
                           const EntityRank          subcell_rank,
                           const CellTopologyData  * side_topology,
                           const EntityVector      & side_nodes );

/**
 * Given an entity, subcell_rank, and subcell_id, return the nodes
 * that make up the subcell in a correct order for the given polarity.
 *
 * \param entity
 * \param subcell_rank
 * \param subcell_indentifier
 * \param subcell_nodes EntityVector output of the subcell nodes
 * \param use_reverse_polarity
 * \return CellTopologyData * of the requested subcell
 */
const CellTopologyData * get_subcell_nodes(
    const Entity     & entity ,
    EntityRank         subcell_rank ,
    unsigned           subcell_identifier ,
    EntityVector     & subcell_nodes
    );



/** \} */

}//namespace mesh
}//namespace stk
#endif
