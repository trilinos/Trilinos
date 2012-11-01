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

namespace stk {
namespace mesh {

/** Create all the edges in the mesh and their relations,
 * i.e. given a mesh with element and node relations
 *  * create all the edges and setup the relations
 *
 * Expensive function that should only be called once
 *
 * This is a parallel collective function (it should be called on all
 * processors at the same time
 *
 * \param mesh A consisent mesh with Element and node relations
 * \param add_parts Newly created edges will be added to the add_parts and
 *                  the root part for the given edge topology.
 */
void create_edges(
    BulkData & mesh,
    const PartVector & add_parts = PartVector()
    );

}
}
#endif

