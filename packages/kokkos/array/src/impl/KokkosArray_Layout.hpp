/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

//----------------------------------------------------------------------------

#ifndef KOKKOSARRAY_LAYOUT_HPP
#define KOKKOSARRAY_LAYOUT_HPP

#include <impl/KokkosArray_StaticAssert.hpp>
#include <impl/KokkosArray_Shape.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

/** \brief  Left-to-right striding of multi-indices (Fortran scheme). */
struct Left { typedef Left array_layout ; };

/** \brief  Right-to-left striding of multi-indices (C or lexigraphical scheme). */
struct Right { typedef Right array_layout ; };

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

/** \brief  Map multi-indices.  Callable from host or device code.
 */
template < class ShapeType , class DeviceSpace >
struct LayoutMap ;

//----------------------------------------------------------------------------

template< class MemorySpace , class T >
inline
size_t stride( const Shape<T,Left> & shape )
{
  enum { Rank = Shape<T,Left>::rank };

  return 0 == Rank ? 0 :
         1 == Rank ? shape.N0 :
         MemorySpace::preferred_stride( shape.value_size , shape.N0 );
}

template< class T >
inline
size_t allocation_size( const Shape<T,Left> & shape )
{
  enum { Rank = Shape<T,Left>::rank };

  return shape.value_size * (
           0 == Rank ? 1 : shape.Stride * (
           1 == Rank ? 1 : shape.N1 * (
           2 == Rank ? 1 : shape.N2 * (
           3 == Rank ? 1 : shape.N3 * (
           4 == Rank ? 1 : shape.N4 * (
           5 == Rank ? 1 : shape.N5 * (
           6 == Rank ? 1 : shape.N6 * (
           7 == Rank ? 1 : shape.N7 ))))))));
}

//----------------------------------------------------------------------------

template< class MemorySpace , class T >
inline
size_t stride( const Shape<T,Right> & shape )
{
  enum { Rank = Shape<T,Right>::rank };

  return
    0 == Rank ? 0 :
    1 == Rank ? 0 :
    MemorySpace::preferred_stride(
      shape.value_size ,
      shape.N1 * ( 2 == Rank ? 1 : shape.N2 * (
                   3 == Rank ? 1 : shape.N3 * (
                   4 == Rank ? 1 : shape.N4 * (
                   5 == Rank ? 1 : shape.N5 * (
                   6 == Rank ? 1 : shape.N6 * (
                   7 == Rank ? 1 : shape.N7 )))))) );
}

template< class T >
inline
size_t allocation_size( const Shape<T,Right> & shape )
{
  enum { Rank = Shape<T,Right>::rank };

  return shape.value_size * (
         0 == Rank ? 1 : shape.N0 * (
         1 == Rank ? 1 : shape.Stride ));
}

} /* namespace Impl */
} /* namespace KokkosArray */

#endif /* #ifndef KOKKOSARRAY_LAYOUT_HPP */


/** \Notes

template< class Type , class MapSpec , class Device = MapSpec >
class Array 
{
  typedef MapSpec::map_scheme    map_scheme ;
  typedef Device ::memory_space  memory_space ;

  typedef ArrayShape< Type > array_shape_type ;
  typedef ArrayIndexMap< array_shape_type , map_scheme , memory_space > index_map_type ;
  typedef Array< Type , MapSpec , Host > HostMirror ;

  index_map_type m_index_map ;

};

template< class Type , class MapSpec , class Device >
class Factory< Array< Type , MapSpec , Host > ,
               Array< Type , MapSpec , Device > , true_type >
{
  typedef Array< Type , MapSpec , Device > input_type ;
  typedef Array< Type , MapSpec , Host >   output_type ;
  typedef ArrayShape< Type >               shape_type ;
  typedef typename MapSpec::map_scheme     map_scheme ;
  typedef typename Device::memory_space    memory_space ;
  typedef ArrayShapeMap< shape_type , map_scheme > shape_map ;

  inline static
  output_type create( const input_type & input ) 
  {
    output_type output ;

    // Retain shape exactly
    output.m_shape = input.m_shape ;

    const size_t alloc_size = shape_map::allocation_size( output.m_shape );

    output.m_ptr_on_device = memory_space::allocate( alloc_size );

    return output ;
  }
};


create_array< Array< Type , MapSpec , Space > >( size_t n0 , size_t n1 )
{
  typedef typename rank_dynamic<Type>::type  dynamic_rank ;
  typedef unsigned_<2>                       input_rank ;

  typedef bool_< is_same< dynamic_rank , input_rank >::value &&
                 other_test::value > ok ;

  return Factory< array_type , input_rank , ok >::create( n0 , n1 );
}

class Factory< Array< Type , MapSpec , Device > , unsigned_<2> , true_type >
{
  typedef ArrayShapeMap< shape_type , map_scheme > shape_map ;

  inline static
  output_type create( n0 , n1 )
  {
    output_type output ;
    output.m_index_map.shape.N0 = n0 ;
    output.m_index_map.shape.N1 = n1 ;
    output.m_index_map.shape.Stride =
      shape_map::template stride<memory_space>( output.m_index_map.shape );

    size_t alloc_size =
      shape_map::allocation_size( output.m_index_map.shape );
  }
};



Array< double[0][0][3] , Host >
Array< double[0][0][3] , Tile<8,8> , Host >


template < class Type , unsigned NX , unsigned NY >
struct IndexMap< Type , Tile<NX,NY> , KOKKOS_MACRO_DEVICE::memory_space >
{
};

*/


