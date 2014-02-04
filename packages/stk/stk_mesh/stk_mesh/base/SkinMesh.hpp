/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_SkinMesh_hpp
#define stk_mesh_SkinMesh_hpp

#include <stk_mesh/base/Types.hpp>      // for PartVector
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Selector; } }

namespace stk { namespace mesh {


/**
 * Skin the entire mesh.
 */
void skin_mesh( BulkData & mesh, Selector const& element_selector, PartVector const& skin_parts = PartVector());

void skin_mesh( BulkData & mesh, PartVector const& skin_parts = PartVector());

}} // namespace stk::mesh
#endif
