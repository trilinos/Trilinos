// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_CUDA_STOCHASTICPRODUCTTENSOR_HPP
#define STOKHOS_CUDA_STOCHASTICPRODUCTTENSOR_HPP

#include "Kokkos_Core.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_BlockCrsMatrix.hpp"
#include "Stokhos_StochasticProductTensor.hpp"

namespace Stokhos {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename TensorScalar ,
          typename MatrixScalar ,
          typename VectorScalar ,
          typename TensorType >
class Multiply<
  BlockCrsMatrix<
    StochasticProductTensor<TensorScalar, TensorType, Kokkos::Cuda>,
    MatrixScalar, Kokkos::Cuda> ,
  Kokkos::View<VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda>,
  Kokkos::View<VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda> >
{
public:

  typedef Kokkos::Cuda                    execution_space ;
  typedef execution_space::size_type  size_type ;

  typedef StochasticProductTensor< TensorScalar , TensorType , Kokkos::Cuda > tensor_type ;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, execution_space > matrix_type ;
  typedef Kokkos::View< VectorScalar** , Kokkos::LayoutLeft , Kokkos::Cuda >           vector_type ;


  static void apply( const matrix_type & A ,
                     const vector_type & x ,
                     const vector_type & y )
  {
    typedef BlockCrsMatrix< TensorType, MatrixScalar , Kokkos::Cuda > base_matrix_type ;

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

} // namespace Stokhos

#endif /* #ifndef STOKHOS_CUDA_STOCHASTICPRODUCTTENSOR_HPP */
