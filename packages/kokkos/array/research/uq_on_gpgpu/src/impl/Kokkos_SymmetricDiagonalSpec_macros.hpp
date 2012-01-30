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

#if ! defined(KOKKOS_MACRO_DEVICE_TEMPLATE_SPECIALIZATION) || \
    ! defined(KOKKOS_MACRO_DEVICE)                  || \
    ! defined(KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION)

#error "Including <impl/Kokkos_SymmetricDiagonalSpec_macros.hpp> without macros defined"

#else

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template<>
class SymmetricDiagonalSpec< KOKKOS_MACRO_DEVICE > {
public:
  typedef KOKKOS_MACRO_DEVICE::size_type size_type ;

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension() const { return m_dimension ; }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type size() const
    { return ( m_dimension * ( m_dimension + 1 ) ) >> 1 ; }

  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type offset( const size_type row , const size_type column ) const
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

  template< typename MatrixValue , typename VectorValue >
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  void multiply( const MatrixValue *       a ,
                 const VectorValue * const x ,
                       VectorValue * const y ) const
  {
    // Number of full diagonals
    // If even number of diagonals then the last diagonal is half length
    const size_type dim_half = ( m_dimension + 1 ) >> 1 ;

    // Multiply the main diagonal (first diagonal)
    for ( size_type j = 0 ; j < m_dimension ; ++j ) {
      y[j] += a[j] * x[j] ; // Contiguous access
    }

    // Multiply remaining full diagionals, each diagonal is accessed twice
    for ( size_type d = 1 ; d < dim_half ; ++d ) {
      size_type kx  = d ;
      size_type kxr = m_dimension - d ;

      a += m_dimension ; // next diagonal

      for ( size_type j = 0 ; j < m_dimension ; ++j ) {
        y[j] += a[j] * x[kx] + a[kxr] * x[kxr]; // Contiguous access
        if ( m_dimension == ++kx )  kx = 0 ;
        if ( m_dimension == ++kxr ) kxr = 0 ;
      }
    }

    // If even number of diagonals then the last diagonal is half-length
    if ( ! ( m_dimension & 01 ) ) {
      size_type kx = dim_half ;

      a += m_dimension ; // next diagonal

      for ( size_type j = 0 ; j < dim_half ; ++j , ++kx ) {
        y[j]  += a[j] * x[kx] ; // Contiguous access
        y[kx] += a[j] * x[j] ;  // Contiguous access
      }
    }
  }

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

} /* namespace Kokkos */

//----------------------------------------------------------------------------

#endif

