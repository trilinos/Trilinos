// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TSQR_IMPL_CUSOLVERHANDLE_FWD_HPP
#define TSQR_IMPL_CUSOLVERHANDLE_FWD_HPP

#include "TpetraTSQR_config.h"
#ifdef HAVE_TPETRATSQR_CUSOLVER

#include <memory>

namespace TSQR {
namespace Impl {

/// \class CuSolverHandle
/// \brief Opaque wrapper for cusolverDnHandle_t (cuSOLVER dense
///   handle instance)
///
/// \note To developers: Do not expose the declaration of this class
///   to downstream code.  Users should only deal with this class by
///   the forward declaration and functions available in this header
///   file.  Do not expose cuSOLVER headers or extern declarations to
///   downstream code.
class CuSolverHandle;

//! Get TSQR's global cuSOLVER handle wrapper.
std::shared_ptr<CuSolverHandle> getCuSolverHandleSingleton();

} // namespace Impl
} // namespace TSQR

#endif // HAVE_TPETRATSQR_CUSOLVER

#endif // TSQR_IMPL_CUSOLVERHANDLE_FWD_HPP
