/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_VIEWTILELEFT_HPP
#define KOKKOSARRAY_VIEWTILELEFT_HPP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class Layout >
struct is_LayoutTileLeft ;

template< unsigned ArgN0 , unsigned ArgN1 >
struct is_LayoutTileLeft< LayoutTileLeft<ArgN0,ArgN1> >
{
  enum { value = true };
  enum { N0 = ArgN0 };
  enum { N1 = ArgN1 };
};


template< class DstViewType >
struct ViewAssignment<
  DstViewType ,
  typename DstViewType::memory_space ,
  typename enable_if< (
    is_LayoutTileLeft< typename DstViewType::array_layout >::value
    &&
    is_same< typename DstViewType::memory_traits , MemoryManaged >::value
  ) >::type >
{
  typedef typename DstViewType::shape_type shape_type ;

private:

  static inline
  void allocate( DstViewType & dst , const std::string & label )
  {
    typedef is_LayoutTileLeft< typename DstViewType::array_layout > layout ;
    typedef typename DstViewType::memory_space  memory_space ;

    ViewAssignment< DstViewType >::decrement( dst.m_ptr_on_device );

    dst.m_stride = 1 ;

    const size_t allocation_count = 
       layout::N0 * layout::N1 * 
       ( ( dst.m_shape.N0 + layout::N0 - 1 ) / layout::N0 ) *
       ( ( dst.m_shape.N1 + layout::N1 - 1 ) / layout::N1 );

    dst.m_ptr_on_device = (typename DstViewType::scalar_type *)
      memory_space::allocate( label ,
                              typeid(typename DstViewType::scalar_type) ,
                              sizeof(typename DstViewType::scalar_type) ,
                              allocation_count );
  }

public:

  ViewAssignment( DstViewType & dst , const std::string & label , const shape_type shape )
  {
    dst.m_shape = shape ;

    allocate( dst , label );
  }

  ViewAssignment( DstViewType & dst , const std::string & label ,
                  const size_t n0 = 0 ,
                  const size_t n1 = 0 ,
                  const size_t n2 = 0 ,
                  const size_t n3 = 0 ,
                  const size_t n4 = 0 ,
                  const size_t n5 = 0 ,
                  const size_t n6 = 0 ,
                  const size_t n7 = 0 )
  {
    shape_type::assign( dst.m_shape, n0, n1, n2, n3, n4, n5, n6, n7 );

    allocate( dst , label );
  }
};

} /* namespace Impl */
} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAY_VIEWTILELEFT_HPP */

