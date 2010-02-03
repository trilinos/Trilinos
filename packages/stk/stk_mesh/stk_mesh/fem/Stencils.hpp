#ifndef stk_mesh_Stencils_hpp
#define stk_mesh_Stencils_hpp

#include <stk_util/util/StaticAssert.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/fem/EntityTypes.hpp>

namespace stk {
namespace mesh {
namespace { // To prevent multiple copies for the linker

enum { ElementStencils_OK =
       StaticAssert< stk::mesh::Node == 0 &&
       stk::mesh::Edge == 1 &&
       stk::mesh::Face == 2 &&
       stk::mesh::Element == 3 >::OK };

//----------------------------------------------------------------------

template< class TopologyTraits >
int element_node_stencil( unsigned , unsigned , unsigned );

template<>
int element_node_stencil<void>( unsigned from_type ,
                                unsigned to_type ,
                                unsigned identifier )
{
  int ordinal = -1 ;

  if ( Element == from_type && Node == to_type ) {
    ordinal = static_cast<int>(identifier);
  }

  return ordinal ;
}

template< class TopologyTraits >
int element_node_stencil( unsigned from_type ,
                          unsigned to_type ,
                          unsigned   identifier )
{
  enum { number_node = TopologyTraits::node_count };

  int ordinal = -1 ;

  if ( Element == from_type &&
       Node    == to_type &&
       identifier < number_node ) {
    ordinal = static_cast<int>(identifier);
  }

  return ordinal ;
}

template< class TopologyTraits >
int element_node_lock_stencil( unsigned , unsigned , unsigned );

template<>
int element_node_lock_stencil<void>( unsigned from_type ,
                                     unsigned to_type ,
                                     unsigned identifier )
{
  int ordinal = -1 ;

  if ( Element == from_type && Node == to_type ) {
    ordinal = (int) identifier ;
  }

  return ordinal ;
}

template< class TopologyTraits >
int element_node_lock_stencil( unsigned from_type ,
                               unsigned to_type ,
                               unsigned identifier )
{
  enum { number_node = TopologyTraits::node_count };

  int ordinal = -1 ;

  if ( Element == from_type &&
       Node    == to_type &&
       identifier < number_node ) {
    ordinal = (int) identifier ;
  }

  return ordinal ;
}

//----------------------------------------------------------------------

} // namespace <empty>
} // namespace mesh
} // namespace stk

#endif //  stk_mesh_Stencils_hpp
