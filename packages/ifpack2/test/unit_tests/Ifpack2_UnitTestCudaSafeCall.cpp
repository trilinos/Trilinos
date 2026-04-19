// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Ifpack2_UnitTestCudaSafeCall.cpp
    \brief Unit tests for \c IFPACK2_IMPL_CUDA_SAFE_CALL when \c KOKKOS_ENABLE_CUDA is defined.
*/

// Ifpack2_ConfigDefs pulls in Kokkos; KOKKOS_ENABLE_CUDA must be visible before the #if below.
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#if defined(KOKKOS_ENABLE_CUDA)

#include "Ifpack2_CudaSafeCall.hpp"

#include <limits>
#include <string>

namespace {

//! Clear any pending CUDA sticky error so later tests see a clean runtime state.
inline void clearCudaLastError() { (void)cudaGetLastError(); }

TEUCHOS_UNIT_TEST(CudaSafeCall, MacroSuccessPath) {
  int nDevices = -1;
  IFPACK2_IMPL_CUDA_SAFE_CALL(cudaGetDeviceCount(&nDevices));
  TEST_INEQUALITY(nDevices, -1);
  TEST_ASSERT(nDevices >= 0);
}

TEUCHOS_UNIT_TEST(CudaSafeCall, DirectCallCudaSuccess) {
  Ifpack2::Impl::cuda_internal_safe_call(cudaSuccess, "cudaSuccess", __FILE__, __LINE__);
}

TEUCHOS_UNIT_TEST(CudaSafeCall, DirectCallInvalidValueThrows) {
  TEST_THROW(Ifpack2::Impl::cuda_internal_safe_call(cudaErrorInvalidValue, "synthetic_InvalidValue",
                                                    __FILE__, __LINE__),
             std::runtime_error);
}

TEUCHOS_UNIT_TEST(CudaSafeCall, DirectCallInvalidDeviceThrows) {
  TEST_THROW(Ifpack2::Impl::cuda_internal_safe_call(cudaErrorInvalidDevice, "synthetic_InvalidDevice",
                                                    __FILE__, __LINE__),
             std::runtime_error);
}

TEUCHOS_UNIT_TEST(CudaSafeCall, DirectCallErrorMessageFormat) {
  bool caught = false;
  try {
    Ifpack2::Impl::cuda_internal_safe_call(cudaErrorInvalidValue, "tag_for_message_test", __FILE__,
                                           __LINE__);
  } catch (const std::runtime_error &err) {
    caught = true;
    const std::string msg(err.what());
    TEST_ASSERT(msg.find("tag_for_message_test") != std::string::npos);
    TEST_ASSERT(msg.find("cudaErrorInvalidValue") != std::string::npos);
    TEST_ASSERT(msg.find("Ifpack2_UnitTestCudaSafeCall.cpp") != std::string::npos);
  }
  TEST_ASSERT(caught);
}

TEUCHOS_UNIT_TEST(CudaSafeCall, MacroInvalidDeviceThrows) {
  const int badDevice = std::numeric_limits<int>::max();
  TEST_THROW(IFPACK2_IMPL_CUDA_SAFE_CALL(cudaSetDevice(badDevice)), std::runtime_error);
  clearCudaLastError();
}

TEUCHOS_UNIT_TEST(CudaSafeCall, MacroErrorMessageContainsStringifiedCall) {
  bool caught         = false;
  const int badDevice = std::numeric_limits<int>::max();
  try {
    IFPACK2_IMPL_CUDA_SAFE_CALL(cudaSetDevice(badDevice));
  } catch (const std::runtime_error &err) {
    caught = true;
    const std::string msg(err.what());
    TEST_ASSERT(msg.find("cudaSetDevice") != std::string::npos);
  }
  TEST_ASSERT(caught);
  clearCudaLastError();
}

TEUCHOS_UNIT_TEST(CudaSafeCall, MacroSuccessAfterClearedError) {
  const int badDevice = std::numeric_limits<int>::max();
  TEST_THROW(IFPACK2_IMPL_CUDA_SAFE_CALL(cudaSetDevice(badDevice)), std::runtime_error);
  clearCudaLastError();
  int nDevices = -1;
  TEST_NOTHROW(IFPACK2_IMPL_CUDA_SAFE_CALL(cudaGetDeviceCount(&nDevices)));
  TEST_ASSERT(nDevices >= 0);
}

}  // namespace

#endif /* KOKKOS_ENABLE_CUDA */
