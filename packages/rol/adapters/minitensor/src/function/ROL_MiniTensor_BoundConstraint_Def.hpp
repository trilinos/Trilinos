// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

namespace ROL {

using Index = minitensor::Index;

template<typename T, Index N>
MiniTensor_BoundConstraint<T, N>::
MiniTensor_BoundConstraint(
    MiniTensorVector<T, N> & lo,
    MiniTensorVector<T, N> & hi
) : Bounds<T>(
      Teuchos::rcp<MiniTensorVector<T, N>>(&lo, false),
      Teuchos::rcp<MiniTensorVector<T, N>>(&hi, false))
{
  return;
}

} // namespace ROL
