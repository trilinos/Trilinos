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

#if ! defined(KOKKOSARRAY_MACRO_DEVICE_TEMPLATE_SPECIALIZATION) || \
    ! defined(KOKKOSARRAY_MACRO_DEVICE)                  || \
    ! defined(KOKKOSARRAY_INLINE_FUNCTION)

#error "Including <impl/KokkosArray_SymmetricDiagonalSpec_macros.hpp> without macros defined"

#else

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

template<>
class SymmetricDiagonalSpec< KOKKOSARRAY_MACRO_DEVICE > {
public:
  typedef KOKKOSARRAY_MACRO_DEVICE::size_type size_type ;

  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension() const { return m_dimension ; }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type matrix_offset( const size_type row , const size_type column ) const
    {
      const int diag_count = 1 + ( m_dimension >> 1 );
      const int diag = (int) column - (int) row ;

      size_type offset = 0 ;

      if ( ( 0 <= diag && diag < diag_count ) || ( diag <= - diag_count ) ) {
        offset = row + m_dimension * ( ( m_dimension + diag ) % m_dimension );
      }
      else {
        offset = column + m_dimension * ( ( m_dimension - diag ) % m_dimension );
      }

      return offset ;
    }

  KOKKOSARRAY_INLINE_FUNCTION
  size_type matrix_size() const
    { return ( m_dimension * ( m_dimension + 1 ) ) >> 1 ; }

  SymmetricDiagonalSpec()
    : m_dimension( 0 ) {}

  SymmetricDiagonalSpec( const SymmetricDiagonalSpec & rhs )
    : m_dimension( rhs.m_dimension ) {}

  SymmetricDiagonalSpec & operator =
    ( const SymmetricDiagonalSpec & rhs )
      { m_dimension = rhs.m_dimension ; return *this ; }

  explicit
  SymmetricDiagonalSpec( const size_type dim )
    : m_dimension( dim ) {}

private:
  size_type m_dimension ;
};

//----------------------------------------------------------------------------

} /* namespace KokkosArray */

//----------------------------------------------------------------------------

#endif

