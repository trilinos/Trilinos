/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
/*
 * CreateFaces.hpp created on Feb 20, 2014 by tcfishe
 */

#ifndef stk_mesh_CreateFaces_hpp
#define stk_mesh_CreateFaces_hpp

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace stk {
  namespace mesh {

    /** Create faces for all elements in "element_selector" and attach them to
     * existing elements.
     *
     * This is a parallel collective function (it should be called on all
     * processors at the same time
     *
     */
    void create_faces(  BulkData & mesh, const Selector & element_selector );

    /** Create faces for all elements in the mesh and attach them to
     * existing elements.
     *
     * This is a parallel collective function (it should be called on all
     * processors at the same time
     *
     */
    inline
    void create_faces( BulkData & mesh )
    {
      create_faces(mesh, mesh.mesh_meta_data().universal_part());
    }
  }
}

#endif
