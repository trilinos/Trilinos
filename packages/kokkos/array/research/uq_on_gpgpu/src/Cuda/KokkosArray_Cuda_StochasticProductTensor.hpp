/*
//@HEADER
// ************************************************************************
// 
//           Kokkos: Node API and Parallel Node Kernels
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
// Questions? Contact H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_CUDA_STOCHASTICPRODUCTTENSOR_HPP
#define KOKKOSARRAY_CUDA_STOCHASTICPRODUCTTENSOR_HPP

#include <utility>
#include <sstream>
#include <stdexcept>
#include <Cuda/KokkosArray_Cuda_Parallel.hpp>
#include <Cuda/KokkosArray_Cuda_ProductTensor.hpp>

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename TensorScalar ,
          class PolynomialType ,
          typename MatrixScalar ,
          typename VectorScalar ,
          template< unsigned , typename , class > class TensorType >
class Multiply<
  BlockCrsMatrix<
    StochasticProductTensor< TensorScalar, PolynomialType, Cuda , TensorType >,
    MatrixScalar , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > ,
  View< VectorScalar** , LayoutLeft , Cuda > >
{
public:

  typedef Cuda                    device_type ;
  typedef device_type::size_type  size_type ;

  typedef StochasticProductTensor< TensorScalar , PolynomialType , Cuda , TensorType > tensor_type ;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, device_type > matrix_type ;
  typedef View< VectorScalar** , LayoutLeft , Cuda >           vector_type ;


  static void apply( const matrix_type & A ,
                     const vector_type & x ,
                     const vector_type & y )
  {
    typedef BlockCrsMatrix< TensorType< 3 , TensorScalar, Cuda >,
                            MatrixScalar , Cuda > base_matrix_type ;

    typedef Multiply< base_matrix_type , vector_type , vector_type >
      base_multiply_type ;

    base_matrix_type base_matrix ;

    base_matrix.values = A.values ;
    base_matrix.graph  = A.graph ;
    base_matrix.block  = A.block.tensor();

    base_multiply_type::apply( base_matrix , x , y );
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOSARRAY_CUDA_STOCHASTICPRODUCTTENSOR_HPP */

