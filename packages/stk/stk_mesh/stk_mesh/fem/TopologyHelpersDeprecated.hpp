/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_TopologyHelpersDeprecated_hpp
#define stk_mesh_TopologyHelpersDeprecated_hpp

#include <sstream>
#include <stdexcept>
#include <Shards_CellTopologyTraits.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/fem/TopologicalMetaData.hpp>

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_bulk_data_element
 *  \{
 */

//----------------------------------------------------------------------
/** \brief Attach a CellTopology to a Part.
 *  There is at most one cell topology allowed.
 */
// DEPRECATED: 09/15/10 FEM TopologicalMetaData refactor
void set_cell_topology_deprecated( Part & , const CellTopologyData * singleton );
// DEPRECATED: 09/15/10 FEM TopologicalMetaData refactor
void set_cell_topology( Part & , const CellTopologyData * singleton );

/** \brief  Attach a CellTopology to a Part.
 *  There is at most one element topology allowed.
 */
// DEPRECATED: 09/15/10 FEM TopologicalMetaData refactor
template< class Traits >
void set_cell_topology_deprecated( Part & p )
{ return set_cell_topology_deprecated( p , shards::getCellTopologyData<Traits>() ); }
// DEPRECATED: 09/15/10 FEM TopologicalMetaData refactor
template< class Traits >
void set_cell_topology( Part & p )
{ return set_cell_topology( p , shards::getCellTopologyData<Traits>() ); }

/** \brief  The the CellTopology attached to the Part, if any */
// DEPRECATED: 09/15/10 FEM TopologicalMetaData refactor
const CellTopologyData * get_cell_topology_deprecated( const Part & );
// DEPRECATED: 09/15/10 FEM TopologicalMetaData refactor
const CellTopologyData * get_cell_topology( const Part & );

/** \brief  The the CellTopology attached to at most one Part of the Bucket */
// DEPRECATED: 09/15/10 FEM TopologicalMetaData refactor
const CellTopologyData * get_cell_topology_deprecated( const Bucket & );
// DEPRECATED: 09/15/10 FEM TopologicalMetaData refactor
const CellTopologyData * get_cell_topology( const Bucket & );

/** \brief  The the CellTopology attached to at most one Part of the Entity */
// DEPRECATED: 09/15/10 FEM TopologicalMetaData refactor
const CellTopologyData * get_cell_topology_deprecated( const Entity & );
// DEPRECATED: 09/15/10 FEM TopologicalMetaData refactor
const CellTopologyData * get_cell_topology( const Entity & );

// DEPRECATED: 09/15/10 FEM TopologicalMetaData refactor
inline
EntityRank element_rank_deprecated(const MetaData & meta)
{
  const TopologicalMetaData * top_data = meta.get_attribute< TopologicalMetaData >();
  return (top_data == NULL) ? (stk::mesh::Element) : (top_data->element_rank);
}

/** \} */

}//namespace mesh
}//namespace stk

#endif // stk_mesh_TopologyHelpersDeprecated_hpp

