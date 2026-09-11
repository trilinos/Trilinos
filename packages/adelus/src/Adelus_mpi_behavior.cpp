/*
//@HEADER
// *****************************************************************************
//                        Adelus
//
// Copyright 2020 NTESS and the Adelus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER
*/

#include "Adelus_mpi_behavior.hpp"

namespace Adelus {

constexpr const std::string_view USE_GPU_AWARE_MPI_ENV_VAR = "ADELUS_USE_GPU_AWARE_MPI";

constexpr bool MpiIsGPUAwareDefault() {
#ifdef ADELUS_USE_GPU_AWARE_MPI
    return true;
#else
    return false;
#endif  // ADELUS_USE_GPU_AWARE_MPI
}

bool probeMpiIsGPUAware() {
  constexpr bool defaultValue = MpiIsGPUAwareDefault();

  static bool value_       = defaultValue;
  static bool initialized_ = false;
  return Teuchos::idempotentlyGetEnvironmentVariable(value_, initialized_, Adelus::USE_GPU_AWARE_MPI_ENV_VAR, defaultValue);
}

}

