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

#ifndef KOKKOSARRAY_SHAPETILELEFT_HPP
#define KOKKOSARRAY_SHAPETILELEFT_HPP

#include <impl/KokkosArray_Shape.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template < unsigned ValueSize , unsigned M, unsigned N, unsigned s0 , unsigned s1 >
struct ShapeMap< Shape<LayoutTileLeft<M,N>,ValueSize,2,s0,s1> >
{
  typedef Shape<LayoutTileLeft<M,N>,ValueSize,2,s0,s1> shape_type ;

  static inline
  size_t allocation_count( const shape_type & shape , unsigned )
  { return M*N*((shape.N0+M-1)/M)*((shape.N1+N-1)/N); }

  static inline
  size_t offset( const shape_type & shape , const size_t i0, const size_t i1 )
  { assert_shape_bounds( shape, i0, i1 ); return shape.offset(i0,i1)  ; }

  template< class MemorySpace >
  static inline
  size_t stride( const shape_type & ) { return 1 ; }
};

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace KokkosArray */

#endif /* #ifndef KOKKOSARRAY_SHAPETILELEFT_HPP */

