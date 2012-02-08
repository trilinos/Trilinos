/*
//@HEADER
// ************************************************************************
// 
//                         Kokkos Array
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
// Questions? Contact H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#if ! defined(KOKKOS_MACRO_DEVICE_TEMPLATE_SPECIALIZATION) || \
    ! defined(KOKKOS_MACRO_DEVICE)                  || \
    ! defined(KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION)

#error "Including <impl/Kokkos_ProductTensor_macros.hpp> without macros defined"

#else

namespace Kokkos {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename ValueType , class PolynomialType >
class ProductTensorFromBases< ValueType, PolynomialType, KOKKOS_MACRO_DEVICE >
{
public:

  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef device_type::size_type  size_type ;
  typedef ValueType               value_type ;

private:

  typedef Kokkos::MDArray<     size_type, device_type>  map_type ;
  typedef Kokkos::MultiVector< value_type,device_type>  vec_type ;

public:

  inline
  ~ProductTensorFromBases() {}

  inline
  ProductTensorFromBases()
    : m_degree_map()
    , m_coord()
    , m_value()
    , m_variable(0)
    , m_dimension(0)
    , m_tensor_count(0)
    {}

  inline
  ProductTensorFromBases( const ProductTensorFromBases & rhs )
    : m_degree_map(   rhs.m_degree_map )
    , m_coord(        rhs.m_coord )
    , m_value(        rhs.m_value )
    , m_variable(     rhs.m_variable )
    , m_dimension(    rhs.m_dimension )
    , m_tensor_count( rhs.m_tensor_count )
    {}

  inline
  ProductTensorFromBases & operator = ( const ProductTensorFromBases & rhs )
  {
    m_degree_map   = rhs.m_degree_map ;
    m_coord        = rhs.m_coord ;
    m_value        = rhs.m_value ;
    m_variable     = rhs.m_variable ;
    m_dimension    = rhs.m_dimension ;
    m_tensor_count = rhs.m_tensor_count ;
    return *this ;
  }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type variable_count() const { return m_variable ; }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension() const { return m_dimension ; }

  template< typename iType >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type variable_degree( const iType & iVariable ) const
    { return m_degree_map( 0 , iVariable ); }

  template< typename iType , typename jType >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type bases_degree( const iType & iBasis , const jType & iVariable ) const
    { return m_degree_map( iBasis + 1 , iVariable ); }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type tensor_count() const
  { return m_tensor_count ; }

  template< typename iType , typename jType >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  size_type tensor_coord( const iType & iEntry , const jType & c ) const
  { return m_coord( iEntry , c ); }

  template< typename iType >
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  const value_type & tensor_value( const iType & entry ) const
  { return m_value[ entry ]; }

  template< typename MatrixValue , typename VectorValue >
  KOKKOS_MACRO_DEVICE_FUNCTION
  void multiply( const MatrixValue *       a ,
                 const VectorValue * const x ,
                       VectorValue * const y ) const
  {
    const size_type nEntry = m_value.length();
    for ( size_type iEntry = 0 ; iEntry < nEntry ; ++iEntry ) {
      const size_type i = m_coord(iEntry,0);
      const size_type j = m_coord(iEntry,1);
      const size_type k = m_coord(iEntry,2);
      const size_type v = m_value(iEntry);

      const bool neq_ij = i != j ;
      const bool neq_jk = j != k ;
      const bool neq_ki = k != i ;

      y[k] += neq_ij ? v * ( a[i] * x[j] + x[i] * a[j] )
                     : v * ( a[i] * x[i] );

      if ( neq_jk ) {
        y[j] += neq_ki ? v * ( a[i] * x[k] + x[i] * a[k] )
                       : v * ( a[i] * x[i] );
      }

      if ( neq_ki && neq_ij ) {
        y[i] += neq_jk ? v * ( a[k] * x[j] + x[k] * a[j] )
                       : v * ( a[j] * x[j] );
      }
    }
  }

private:

  Kokkos::MDArray< size_type , device_type >            m_degree_map ;
  Kokkos::MDArray< size_type , device_type >            m_coord ;
  Kokkos::Impl::MemoryView< value_type , device_type >  m_value ;
  size_type                                             m_variable ;
  size_type                                             m_dimension ;
  size_type                                             m_tensor_count ;

  template< class T , class I >
  friend class Impl::CreateProductTensorFromBases ;
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Kokkos

#endif

