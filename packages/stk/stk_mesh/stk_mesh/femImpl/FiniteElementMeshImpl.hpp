/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_femImpl_FiniteElementMeshImpl_hpp
#define stk_femImpl_FiniteElementMeshImpl_hpp

#include <string>
#include <vector>

namespace stk {
namespace mesh {
namespace impl {

void verify_spatial_dimension( const char * , unsigned spatial_dimension );

std::vector<std::string>
finite_element_mesh_entity_rank_names( unsigned spatial_dimension );

}
}
}

#endif

