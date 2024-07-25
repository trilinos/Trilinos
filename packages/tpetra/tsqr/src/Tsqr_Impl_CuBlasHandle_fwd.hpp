// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TSQR_IMPL_CUBLASHANDLE_FWD_HPP
#define TSQR_IMPL_CUBLASHANDLE_FWD_HPP

#include "TpetraTSQR_config.h"
#ifdef HAVE_TPETRATSQR_CUBLAS

#include <memory>

namespace TSQR {
namespace Impl {

/// \class CuBlasHandle
/// \brief Opaque wrapper for cublasHandle_t (cuBLAS handle instance)
///
/// \note To developers: Do not expose the declaration of this class
///   to downstream code.  Users should only deal with this class by
///   the forward declaration and functions available in this header
///   file.  Do not expose cuBLAS headers or extern declarations to
///   downstream code.
class CuBlasHandle;

//! Get TSQR's global cuBLAS handle wrapper.
std::shared_ptr<CuBlasHandle> getCuBlasHandleSingleton();

} // namespace Impl
} // namespace TSQR

#endif // HAVE_TPETRATSQR_CUBLAS

#endif // TSQR_IMPL_CUBLASHANDLE_FWD_HPP
