// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_CUDA_LEXICOGRAPHIC_BLOCK_SPARSE_3_TENSOR_HPP
#define STOKHOS_CUDA_LEXICOGRAPHIC_BLOCK_SPARSE_3_TENSOR_HPP

#include <iostream>

#include "Kokkos_Core.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_BlockCrsMatrix.hpp"
#include "Stokhos_LexicographicBlockSparse3Tensor.hpp"

#include "cuda_profiler_api.h"

namespace Stokhos {

//----------------------------------------------------------------------------

template< typename TensorScalar ,
          typename MatrixScalar ,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix< LexicographicBlockSparse3Tensor< TensorScalar, Kokkos::Cuda >,
                  MatrixScalar, Kokkos::Cuda >,
  Kokkos::View<VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda>,
  Kokkos::View<VectorScalar**, Kokkos::LayoutLeft, Kokkos::Cuda> >
{
public:

  typedef Kokkos::Cuda execution_space ;
  typedef execution_space::size_type size_type ;

  typedef LexicographicBlockSparse3Tensor< TensorScalar , execution_space > tensor_type ;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, execution_space > matrix_type ;
  typedef Kokkos::View< VectorScalar** ,
                             Kokkos::LayoutLeft ,
                             Kokkos::Cuda > vector_type ;



  //------------------------------------

  static void apply( const matrix_type & A ,
                     const vector_type & x ,
                     const vector_type & y )
  {
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_CUDA_FLAT_SPARSE_3_TENSOR_HPP */
