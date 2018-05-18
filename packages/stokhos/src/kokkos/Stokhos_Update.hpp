// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_UPDATE_HPP
#define STOKHOS_UPDATE_HPP

namespace Stokhos {

template< typename ValueType, typename VectorType >
class Update
{
public:
  typedef VectorType                        vector_type;
  typedef ValueType                         value_type;
  typedef typename vector_type::execution_space execution_space;
  typedef typename execution_space::size_type   size_type;


        vector_type  m_x;
  const vector_type  m_y;
  const value_type   m_alpha;
  const value_type   m_beta;

  Update(const value_type& alpha, vector_type& x,
         const value_type& beta,  const vector_type& y)
  : m_x( x )
  , m_y( y )
  , m_alpha( alpha )
  , m_beta( beta )
  {}

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type iRow ) const
  {
    m_x(iRow) = m_alpha * m_x(iRow) + m_beta * m_y(iRow) ;
  }

  static void apply(const value_type& alpha, vector_type& x,
                    const value_type& beta,  const vector_type& y)
  {
    const size_t row_count = x.extent(0);
    Kokkos::parallel_for( row_count , Update(alpha,x,beta,y) );
  }
};

template <typename ValueType, typename VectorType>
void update(const ValueType& alpha, VectorType& x,
            const ValueType& beta,  const VectorType& y)
{
  Update<ValueType,VectorType>::apply( alpha , x , beta, y );
}

} // namespace Stokhos

#endif
