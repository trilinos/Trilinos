/*
//@HEADER
// ************************************************************************
// 
//    KokkosArray: Manycore Performance-Portable Multidimensional Arrays
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

#if ! defined(KOKKOSARRAY_MACRO_DEVICE_TEMPLATE_SPECIALIZATION) || \
    ! defined(KOKKOSARRAY_MACRO_DEVICE)                  || \
    ! defined(KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION)

#error "Including <impl/KokkosArray_ProductTensor_macros.hpp> without macros defined"

#else

namespace KokkosArray {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename ValueType , class PolynomialType ,
          template< unsigned , typename , class > class TensorType >
class StochasticProductTensor< ValueType, PolynomialType, KOKKOSARRAY_MACRO_DEVICE , TensorType >
{
public:

  typedef KOKKOSARRAY_MACRO_DEVICE     device_type ;
  typedef device_type::size_type  size_type ;
  typedef ValueType               value_type ;
  typedef TensorType< 3 , value_type , device_type > tensor_type ;

  inline
  ~StochasticProductTensor() {}

  inline
  StochasticProductTensor()
    : m_tensor()
    , m_degree_map()
    , m_variable(0)
    {}

  inline
  StochasticProductTensor( const StochasticProductTensor & rhs )
    : m_tensor(       rhs.m_tensor )
    , m_degree_map(   rhs.m_degree_map )
    , m_variable(     rhs.m_variable )
    {}

  inline
  StochasticProductTensor & operator = ( const StochasticProductTensor & rhs )
  {
    m_tensor       = rhs.m_tensor ;
    m_degree_map   = rhs.m_degree_map ;
    m_variable     = rhs.m_variable ;
    return *this ;
  }

  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  const tensor_type & tensor() const { return m_tensor ; }

  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension() const { return m_tensor.dimension(); }

  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type variable_count() const { return m_variable ; }

  template< typename iType >
  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type variable_degree( const iType & iVariable ) const
    { return m_degree_map( 0 , iVariable ); }

  template< typename iType , typename jType >
  inline
  KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type bases_degree( const iType & iBasis , const jType & iVariable ) const
    { return m_degree_map( iBasis + 1 , iVariable ); }

private:

  tensor_type                                        m_tensor ;
  KokkosArray::View< size_type** , device_type >  m_degree_map ;
  size_type                                          m_variable ;

  template< class T , class I >
  friend class Impl::CreateSparseProductTensor ;
};

//----------------------------------------------------------------------------

} // namespace KokkosArray

#endif

