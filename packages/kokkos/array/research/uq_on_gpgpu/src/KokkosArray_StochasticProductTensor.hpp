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

#ifndef KOKKOSARRAY_STOCHASTICPRODUCTTENSOR_HPP
#define KOKKOSARRAY_STOCHASTICPRODUCTTENSOR_HPP

#include <KokkosArray_ProductTensor.hpp>

namespace KokkosArray {

//----------------------------------------------------------------------------

/** \brief  Bases defined by combinatorial product of polynomial bases.
 *
 *  Bases: \prod_{j=0}^{N-1} P_k(x) \forall j and k \in M(j)
 *  Where:  P_k      is a polynomial of degree k
 *  Where: <P_a,P_b> is the the integral on [-1,1]
 *  Where: <P_a,P_b> is the Kronecker delta \delta_{a,b}
 *                   thus the polynomials are normalized
 *                   with respect to this inner product.
 *
 *  Where: N = the number of variables expanded via polynomial bases
 *  Where: M(j) = the degree of a particular variable
 *
 *  Where: \psi_I(x) = is one basis function and I is a multi-index
 *                     of rank N, denoting one function from each 
 *                     variable's polynomial bases.
 *
 *  Were: <\psi_I,\psi_J,\psi_K> is the integral on [-1,1]
 *
 *  The bases space is sparse due to orthogonality within the
 *  expansion.
 */
template< typename ValueType , class PolynomialType , class Device ,
          template< unsigned , typename , class >
          class TensorType = SparseProductTensor >
class StochasticProductTensor {
public:

  typedef Device                                      device_type ;
  typedef typename device_type::size_type             size_type ;
  typedef ValueType                                   value_type ;
  typedef TensorType< 3 , value_type , device_type >  tensor_type ;

private:

  tensor_type                                     m_tensor ;
  KokkosArray::View< size_type** , device_type >  m_degree_map ;
  size_type                                       m_variable ;

  template< class T , class I >
  friend class Impl::CreateSparseProductTensor ;

public:

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

  KOKKOSARRAY_INLINE_FUNCTION
  const tensor_type & tensor() const { return m_tensor ; }

  /** \brief  Dimension: number of bases and
   *          length of the vector block (and tensor).
   */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type dimension() const { return m_tensor.dimension(); }

  /** \brief  How many variables are being expanded. */
  KOKKOSARRAY_INLINE_FUNCTION
  size_type variable_count() const { return m_variable ; }

  /** \brief  Polynomial degree of a given variable */
  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  size_type variable_degree( const iType & iVariable ) const
    { return m_degree_map( 0 , iVariable ); }

  /** \brief  Basis function 'iBasis' is the product of
   *          'variable_count()' polynomials.  Return the
   *          polynomial degree of component 'iVariable'.
   */
  template< typename iType , typename jType >
  KOKKOSARRAY_INLINE_FUNCTION
  size_type bases_degree( const iType & iBasis , const jType & iVariable ) const
    { return m_degree_map( iBasis + 1 , iVariable ); }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace KokkosArray

#include <impl/KokkosArray_StochasticProductTensor_create.hpp>

#endif /* #ifndef KOKKOSARRAY_STOCHASTICPRODUCTTENSOR_HPP */


