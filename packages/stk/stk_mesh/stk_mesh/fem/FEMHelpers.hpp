/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_FEMHelpers_hpp
#define stk_mesh_FEMHelpers_hpp

#include <stk_mesh/base/Types.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/CellTopology.hpp>

namespace stk {
namespace mesh {

class Bucket;
class Entity;

namespace fem {

/** \addtogroup stk_mesh_bulk_data_element
 *  \{
 */

//----------------------------------------------------------------------
/** \brief  Declare an element member of a Part with a CellTopology
 *          and nodes conformal to that topology.
 */
Entity & declare_element( BulkData & mesh ,
                          Part & part ,
                          const EntityId elem_id ,
                          const EntityId node_id[] );


/** \brief  Create (or find) an element side.
 *
 *  The element must be a member of a Part with a CellTopology.
 */
Entity & declare_element_side( BulkData & mesh ,
                               const stk::mesh::EntityId global_side_id ,
                               Entity & elem ,
                               const unsigned local_side_id ,
                               Part * part = NULL );

/** \brief  Create (or find) an element side.
 *
 *  The element must be a member of a Part with a CellTopology.
 */
Entity & declare_element_side( Entity & elem ,
                               Entity & side ,
                               const unsigned local_side_id ,
                               Part * part = NULL );

// TODO: Move to FEMBulkData once we create it. Remove the _new part of the func name
CellTopology get_cell_topology_new( const Bucket & bucket);

// TODO: Move to FEMBulkData once we create it. Remove the _new part of the func name
CellTopology get_cell_topology_new( const Entity & entity);

/** \brief  Declare a part with a given cell topology. This is just a convenient
            function that wraps FEMMetaData's declare_part.
 */
template< class Top >
Part &declare_part(FEMMetaData& meta_data, const std::string &name) {
  return meta_data.declare_part(name, shards::getCellTopologyData<Top>());
}

/** \} */

} //namespace fem
} //namespace mesh
} //namespace stk
#endif
