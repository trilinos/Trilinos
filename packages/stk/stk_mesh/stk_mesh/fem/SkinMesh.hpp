/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_SkinMesh_hpp
#define stk_mesh_SkinMesh_hpp


namespace stk {
namespace mesh {

class BulkData;
class Part;

void skin_mesh( BulkData & mesh, unsigned closure_rank, Part * part = NULL );

}
}
#endif
