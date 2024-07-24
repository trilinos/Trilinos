// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TSQR_IMPL_CUBLASHANDLE_HPP
#define TSQR_IMPL_CUBLASHANDLE_HPP

#include "Tsqr_Impl_CuBlasHandle_fwd.hpp"
#ifdef HAVE_TPETRATSQR_CUBLAS
#include <cublas_v2.h>

namespace TSQR {
namespace Impl {

class CuBlasHandle {
public:
  CuBlasHandle () = delete;
  CuBlasHandle (const CuBlasHandle&) = delete;
  CuBlasHandle& operator= (const CuBlasHandle&) = delete;
  CuBlasHandle (CuBlasHandle&&) = delete;
  CuBlasHandle& operator= (CuBlasHandle&&) = delete;

  CuBlasHandle (cublasHandle_t handle);
  cublasHandle_t getHandle () const;

private:
  // cublasHandle_t is actually a pointer type.
  cublasHandle_t handle_ {nullptr};
};

} // namespace Impl
} // namespace TSQR

#endif // HAVE_TPETRATSQR_CUBLAS

#endif // TSQR_IMPL_CUBLASHANDLE_HPP
