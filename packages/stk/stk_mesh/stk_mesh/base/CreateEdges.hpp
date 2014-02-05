/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_CreateEdges_hpp
#define stk_mesh_CreateEdges_hpp

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace stk {
namespace mesh {

/** Create all the edges in the mesh and attach them to
 * existing elements and defined faces.
 *
 * This is a parallel collective function (it should be called on all
 * processors at the same time
 *
 */
void create_edges(  BulkData & mesh, const Selector & element_selector );

inline
void create_edges( BulkData & mesh )
{
  create_edges(mesh, mesh.mesh_meta_data().universal_part());
}

}
}
#endif

