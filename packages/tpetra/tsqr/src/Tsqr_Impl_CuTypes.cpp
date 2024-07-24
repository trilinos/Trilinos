// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tsqr_Impl_CuTypes.hpp"
#if defined(HAVE_TPETRATSQR_CUBLAS) && defined(HAVE_TPETRATSQR_CUSOLVER)

namespace TSQR {
namespace Impl {

cublasSideMode_t cuBlasSide (const char side)
{
  if (side == 'L' || side == 'l') {
    return CUBLAS_SIDE_LEFT;
  }
  else {
    return CUBLAS_SIDE_RIGHT;
  }
}

cublasOperation_t cuBlasTrans (const char trans)
{
  if (trans == 'C' || trans == 'c') {
    return CUBLAS_OP_C;
  }
  else if (trans == 'T' || trans == 't') {
    return CUBLAS_OP_T;
  }
  else {
    return CUBLAS_OP_N;
  }
}

} // namespace Impl
} // namespace TSQR

#endif // HAVE_TPETRATSQR_CUBLAS && HAVE_TPETRATSQR_CUSOLVER
