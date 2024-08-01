// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TSQR_IMPL_CUSOLVERHANDLE_HPP
#define TSQR_IMPL_CUSOLVERHANDLE_HPP

#include "Tsqr_Impl_CuSolverHandle_fwd.hpp"
#ifdef HAVE_TPETRATSQR_CUSOLVER
#include <cusolverDn.h>

namespace TSQR {
namespace Impl {

class CuSolverHandle {
public:
  CuSolverHandle () = delete;
  CuSolverHandle (const CuSolverHandle&) = delete;
  CuSolverHandle& operator= (const CuSolverHandle&) = delete;
  CuSolverHandle (CuSolverHandle&&) = delete;
  CuSolverHandle& operator= (CuSolverHandle&&) = delete;

  CuSolverHandle (cusolverDnHandle_t handle);
  cusolverDnHandle_t getHandle () const;

private:
  // cusolverDnHandle_t is actually a pointer type.
  cusolverDnHandle_t handle_ {nullptr};
};

} // namespace Impl
} // namespace TSQR

#endif // HAVE_TPETRATSQR_CUSOLVER

#endif // TSQR_IMPL_CUSOLVERHANDLE_HPP
