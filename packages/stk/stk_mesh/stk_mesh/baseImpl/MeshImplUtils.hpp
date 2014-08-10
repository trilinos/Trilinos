/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_MeshUtils_hpp
#define stk_mesh_MeshUtils_hpp

//----------------------------------------------------------------------

#include <stk_mesh/base/Types.hpp>

#include <vector>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace impl {

//----------------------------------------------------------------------
//These functions are not part of the public API of stk-mesh.
//They are intended for use internally in the implementation of
//stk-mesh capabilities.
//----------------------------------------------------------------------

void find_elements_these_nodes_have_in_common(BulkData& mesh, unsigned numNodes, const Entity* nodes, std::vector<Entity>& elems);

void find_locally_owned_elements_these_nodes_have_in_common(BulkData& mesh, unsigned numNodes, const Entity* nodes, std::vector<Entity>& elems);

bool find_element_edge_ordinal_and_equivalent_nodes(BulkData& mesh, Entity element, unsigned numEdgeNodes, const Entity* edgeNodes, unsigned& elemEdgeOrdinal, Entity* elemEdgeNodes);

} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

