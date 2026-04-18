// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Ifpack2_UnitTestCudaSafeCall.cpp
    \brief Unit tests for \c IFPACK2_IMPL_CUDA_SAFE_CALL (CUDA builds only).
*/

#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include "Ifpack2_CudaSafeCall.hpp"

namespace {

TEUCHOS_UNIT_TEST(CudaSafeCall, MacroSuccessPath) {
  int nDevices = -1;
  IFPACK2_IMPL_CUDA_SAFE_CALL(cudaGetDeviceCount(&nDevices));
  TEST_INEQUALITY(nDevices, -1);
  TEST_ASSERT(nDevices >= 0);
}

TEUCHOS_UNIT_TEST(CudaSafeCall, DirectCallCudaSuccess) {
  Ifpack2::Impl::cuda_internal_safe_call(cudaSuccess, "cudaSuccess", __FILE__, __LINE__);
}

}  // namespace
