/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#ifndef KOKKOS_ARRAYBOUNDS_HPP
#define KOKKOS_ARRAYBOUNDS_HPP

#include <cstddef>

namespace Kokkos {
namespace Impl {

void multivector_require_equal_dimension(
  size_t length_x , size_t count_x ,
  size_t length_y , size_t count_y );

template< class MultiVectorX , class MultiVectorY >
void multivector_require_equal_dimension( const MultiVectorX & x ,
                                          const MultiVectorY & y )
{
  multivector_require_equal_dimension( x.length() , x.count() ,
                                       y.length() , y.count() );
}

void multivectorview_range_error( size_t , size_t , size_t );

size_t mdarray_deduce_rank( size_t , size_t , size_t , size_t ,
                            size_t , size_t , size_t , size_t );

void require_less( size_t , size_t );

void mdarray_require_dimension(
  size_t n_rank ,
  size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
  size_t n4 , size_t n5 , size_t n6 , size_t n7 ,
  size_t i_rank ,
  size_t i0 , size_t i1 , size_t i2 , size_t i3 ,
  size_t i4 , size_t i5 , size_t i6 , size_t i7 );

void mdarray_require_equal_dimension(
  size_t n_rank , const size_t n_dims[] ,
  size_t m_rank , const size_t m_dims[] );

template < class MDArrayType1 , class MDArrayType2 >
inline
void mdarray_require_equal_dimension( const MDArrayType1 & array1 ,
                                      const MDArrayType2 & array2 )
{
  const size_t array1_rank = array1.rank();
  const size_t array2_rank = array2.rank();
  size_t array1_dims[8] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };
  size_t array2_dims[8] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };
  array1.dimensions( array1_dims );
  array2.dimensions( array2_dims );
  mdarray_require_equal_dimension( array1_rank, array1_dims,
                                   array2_rank, array2_dims);
}

} // namespace Impl
} // namespace Kokkos

#endif /* KOKKOS_ARRAYBOUNDS_HPP */

