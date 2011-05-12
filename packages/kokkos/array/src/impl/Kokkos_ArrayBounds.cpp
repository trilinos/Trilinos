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
#include <sstream>
#include <stdexcept>

namespace Kokkos {
namespace Impl {

size_t mdarray_deduce_rank( size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
                            size_t n4 , size_t n5 , size_t n6 , size_t n7 )
{
  enum { MAX_RANK = 8 };

  const size_t dim[ MAX_RANK ] = { n0 , n1 , n2 , n3 , n4 , n5 , n6 , n7 };

  size_t rank = 0 , r = 0 ;

  for ( ; rank < MAX_RANK && 0 < dim[rank] ; ++rank );

  for ( r = rank ; r < MAX_RANK && 0 == dim[r] ; ++r );

  if ( MAX_RANK != r ) {
    std::ostringstream msg ;
    msg << "Kokkos::Impl::mdarray_deduce_rank( " << dim[0] ;
    for ( r = 1 ; r < MAX_RANK ; ++r ) {
       msg << " , " << dim[r] ;
    }
    msg << " ) FAILED : Zero dimension inside non-zero dimension" ;
    throw std::runtime_error( msg.str() );
  }

  return rank ;
}

void require_less( size_t i , size_t j )
{
  if ( j <= i ) {
    std::ostringstream msg ;
    msg << "Kokkos::Impl::require_less( " << i << " , " << j << " ) FAILED" ;
    throw std::runtime_error( msg.str() );
  }
}

void mdarray_require_dimension(
  size_t n_rank ,
  size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
  size_t n4 , size_t n5 , size_t n6 , size_t n7 ,
  size_t i_rank ,
  size_t i0 , size_t i1 , size_t i2 , size_t i3 ,
  size_t i4 , size_t i5 , size_t i6 , size_t i7 )
{
  enum { MAX_RANK = 8 };

  const size_t dim[ MAX_RANK ] = { n0 , n1 , n2 , n3 , n4 , n5 , n6 , n7 };
  const size_t ind[ MAX_RANK ] = { i0 , i1 , i2 , i3 , i4 , i5 , i6 , i7 };

  bool ok = n_rank <= MAX_RANK && n_rank == i_rank ;

  for ( size_t r = 0 ; ok && r < n_rank ; ++r ) {
    ok = ind[r] < dim[r] ;
  }

  if ( ! ok ) {
    std::ostringstream msg ;
    msg << "Kokkos::Impl::require_dimension( dimension( " ;
    for ( size_t r = 0 ; r < MAX_RANK && r < n_rank ; ++r ) {
      if ( r ) msg << " , " ;
      msg << dim[r] ;
    }
    msg << " ) , index( " ;
    for ( size_t r = 0 ; r < MAX_RANK && r < i_rank ; ++r ) {
      if ( r ) msg << " , " ;
      msg << ind[r] ;
    }
    msg << " ) ) FAILED Index is incompatible with dimension" ;

    throw std::runtime_error( msg.str() );
  }
}

} // namespace Impl
} // namespace Kokkos

#endif /* KOKKOS_ARRAYBOUNDS_HPP */

