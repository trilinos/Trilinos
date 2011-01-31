/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_Stencils_deprecated_hpp
#define stk_mesh_Stencils_deprecated_hpp

#include <stk_util/util/StaticAssert.hpp>
#include <stk_mesh/base/Types.hpp>
// TO BE REMOVED AFTER FEM refactor is finished:
#include <stk_mesh/fem/EntityRanks.hpp> 

namespace stk {
namespace mesh {
namespace { // To prevent multiple copies for the linker

enum EntityRankEnum3D {
  Node3D                = 0 ,
  Edge3D                = 1 ,
  Face3D                = 2 ,
  Element3D             = 3 
};

//----------------------------------------------------------------------

template< class TopologyTraits >
int element_node_stencil( EntityRank , EntityRank , unsigned );

template<>
int element_node_stencil<void>( EntityRank from_type ,
                                EntityRank to_type ,
                                unsigned identifier )
{
  int ordinal = -1 ;

  if ( Element3D == from_type && Node3D == to_type ) {
    ordinal = static_cast<int>(identifier);
  }

  return ordinal ;
}

template< class TopologyTraits >
int element_node_stencil( EntityRank from_type ,
                          EntityRank to_type ,
                          unsigned   identifier )
{
  enum { number_node = TopologyTraits::node_count };

  int ordinal = -1 ;

  if ( Element3D == from_type &&
       Node3D    == to_type &&
       identifier < number_node ) {
    ordinal = static_cast<int>(identifier);
  }

  return ordinal ;
}

template< class TopologyTraits >
int entity_node_stencil( EntityRank , EntityRank , unsigned );

template<>
int entity_node_stencil<void>( EntityRank from_entity_rank ,
                               EntityRank to_entity_rank ,
                               unsigned identifier )
{
  int ordinal = -1 ;

  if ( Node3D == to_entity_rank ) {
    ordinal = static_cast<int>(identifier);
  }

  return ordinal ;
}

template< class TopologyTraits >
int entity_node_stencil( EntityRank from_entity_rank ,
                         EntityRank to_entity_rank ,
                         unsigned   identifier )
{
  enum { number_node = TopologyTraits::node_count };

  int ordinal = -1 ;

  if ( Node3D    == to_entity_rank &&
       identifier < number_node ) {
    ordinal = static_cast<int>(identifier);
  }

  return ordinal ;
}

//----------------------------------------------------------------------

} // namespace <empty>
} // namespace mesh
} // namespace stk

#endif //  stk_mesh_Stencils_deprecated_hpp
