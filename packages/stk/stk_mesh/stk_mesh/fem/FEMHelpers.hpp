/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_FEMHelpers_hpp
#define stk_mesh_FEMHelpers_hpp

#include <sstream>
#include <stdexcept>
#include <Shards_CellTopologyTraits.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/fem/FEMInterface.hpp>

namespace stk {
namespace mesh {
namespace fem {

/** \addtogroup stk_mesh_bulk_data_element
 *  \{
 */

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
  stk::mesh::fem::FEMMetaData & fem_meta = stk::mesh::fem::FEMMetaData::get(mesh);
  const CellTopologyData * const top = fem_meta.get_cell_topology( part ).getCellTopologyData();

  ThrowErrorMsgIf(top == NULL,
                  "Part " << part.name() << " does not have a local topology");

  PartVector empty ;
  PartVector add( 1 ); add[0] = & part ;

  const EntityRank entity_rank = fem_meta.element_rank();

  Entity & elem = mesh.declare_entity( entity_rank, elem_id, add );

  const EntityRank node_rank = fem_meta.node_rank();

  for ( unsigned i = 0 ; i < top->node_count ; ++i ) {
    //declare node if it doesn't already exist
    Entity * node = mesh.get_entity( node_rank , node_id[i]);
    if ( NULL == node) {
      node = & mesh.declare_entity( node_rank , node_id[i], empty );
    }

    mesh.declare_relation( elem , *node , i );
  }
  return elem ;
}

/** \} */

}//namespace fem
}//namespace mesh
}//namespace stk
#endif
