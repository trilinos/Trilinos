// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Mechanics_i_h)
#define MiniTensor_Mechanics_i_h

namespace minitensor {

//
// R^N volumetric part of 2nd-order tensor
// \return \f$ \frac{1}{N} \mathrm{tr}\:(A) I \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
vol(Tensor<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  T const
  theta = (1.0/dimension) * trace(A);

  return theta * eye<T, N>(dimension);
}

//
// R^N deviatoric part of 2nd-order tensor
// \return \f$ A - vol(A) \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
dev(Tensor<T, N> const & A)
{
  return A - vol(A);
}

} // namespace minitensor

#endif // MiniTensor_Mechanics_i_h
