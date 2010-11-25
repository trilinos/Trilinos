/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_TopologyHelpersDeprecated_hpp
#define stk_mesh_TopologyHelpersDeprecated_hpp

#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/FEMInterface.hpp>

// Migration to new Finite Element Mesh Cell Topology
//
// To turn off all the deprecated code:
//   *  Enable the following ifdef: SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
//   *  This can be enabled in Sierra with the following bake command:
//      bake stk_mesh -- cxxflags="-DSKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS"
//
// Wherever you call set_cell_topology:
//   *  #include<stk_mesh/fem/TopologicalMetaData.hpp>
//   *  Instantiate a TopologicalMetaData with MetaData and the spatial dimension
//   *  Call TopologicalMetaData.declare_part<shards::CELL_TOPOLOGY>("name")
//      which will declare the part and set the cell topology.
// Wherever you call get_cell_topology(...).getTopologyData():
//   *  Call the static function:  TopologicalMetaData::get_cell_topology(...).getTopologyData()
// Wherever you #include<stk_mesh/fem/FieldDeclarations.hpp>,
//   *  Copy these functions out into your own application and stop #including this file
// Wherever you #include<stk_mesh/fem/FieldTraits.hpp> (note:  fem, NOT base):
//   *  Change the #include to #include<stk_mesh/fem/CoordinateSystems.hpp>
//

// Note the following information about the changes:
//   *  EntityRank has been changed from a compile time enum to a runtime
//      value.
//   *  Part rank consistency with cell topology dimension is now enforced at
//      construction.
//   *  TopologyHelpers used to provide set_cell_topology and
//      get_cell_topology.  These are now provided by TopologicalMetaData.
//   *  The mechanism for storing cell topology on parts has changed.
//      Previously, the cell topology was stored on each part.  Now it is stored
//      as a map in the TopologicalMetaData which is stored as an attribute on
//      MetaData.
//   *  Previously, get_cell_topology was a free function, it is now a static
//      function on TopologicalMetaData.
//   *  Previously, set_cell_topology was a free function, it is now a regular
//      member function on TopologicalMetaData, which means you'll have to create
//      one of these objects in order to set cell topology on a part.
//   *  Previously, cell topologies were set on parts after the part was
//      declared.  Now, you must construct a part with a cell topology through
//      TopologicalMetaData.
//   *  stk_mesh/fem/FieldDeclarations.hpp is now deprecated due to its tight
//      coupling to specific application domains.  If you use any of these
//      functions, we advise you to copy them out into your own application
//      namespace.  They are very simple convenience wrappers around
//      stk::mesh::put_field.
//   *  stk_mesh/fem/FieldTraits.hpp is changing name to
//      stk_mesh/fem/CoordinateSystems.hpp and the FieldTraits file will remain
//      for one Trilinos release.
//   *  Note:  The Trilinos release after 10.6, presumably 10.8, mark when we
//      will delete the deprecated functionality.
//      stk_mesh/fem/FieldDeclarations.hpp and stk_mesh/fem/FieldTraits.hpp will
//      also be removed at this time.

#include <sstream>
#include <stdexcept>
#include <Shards_CellTopologyTraits.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace stk {
namespace mesh {


/** \addtogroup stk_mesh_bulk_data_element
 *  \{
 */

//----------------------------------------------------------------------
/** \brief Attach a CellTopology to a Part.
 *  There is at most one cell topology allowed.
 */
// DEPRECATED: 09/15/10 FEM refactor
void set_cell_topology_deprecated( Part & , const CellTopologyData * singleton );
// DEPRECATED: 09/15/10 FEM refactor
void set_cell_topology( Part & , const CellTopologyData * singleton );

/** \brief  Attach a CellTopology to a Part.
 *  There is at most one element topology allowed.
 */
// DEPRECATED: 09/15/10 FEM refactor
template< class Traits >
void set_cell_topology_deprecated( Part & p )
{ return set_cell_topology_deprecated( p , shards::getCellTopologyData<Traits>() ); }
// DEPRECATED: 09/15/10 FEM refactor
template< class Traits >
void set_cell_topology( Part & p )
{ return set_cell_topology( p , shards::getCellTopologyData<Traits>() ); }

/** \brief  The the CellTopology attached to the Part, if any */
// DEPRECATED: 09/15/10 FEM refactor
const CellTopologyData * get_cell_topology_deprecated( const Part & );
// DEPRECATED: 09/15/10 FEM refactor
const CellTopologyData * get_cell_topology( const Part & );

/** \brief  The the CellTopology attached to at most one Part of the Bucket */
// DEPRECATED: 09/15/10 FEM refactor
const CellTopologyData * get_cell_topology_deprecated( const Bucket & );
// DEPRECATED: 09/15/10 FEM refactor
const CellTopologyData * get_cell_topology( const Bucket & );

/** \brief  The the CellTopology attached to at most one Part of the Entity */
// DEPRECATED: 09/15/10 FEM refactor
const CellTopologyData * get_cell_topology_deprecated( const Entity & );
// DEPRECATED: 09/15/10 FEM refactor
const CellTopologyData * get_cell_topology( const Entity & );

// DEPRECATED: 09/15/10 FEM refactor
inline
EntityRank element_rank_deprecated(const MetaData & meta)
{
  const fem::FEMInterface* fem = meta.get_attribute<fem::FEMInterface>();
  if (fem) {
    return fem::element_rank(*fem);
  }
  else {
    return stk::mesh::Element;
  }
}

/** \} */


}//namespace mesh
}//namespace stk

#endif //  SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS

#endif // stk_mesh_TopologyHelpersDeprecated_hpp

