//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef _KOKKOSKERNELS_TPLS_VERSIONS_HPP
#define _KOKKOSKERNELS_TPLS_VERSIONS_HPP

#include "KokkosKernels_config.h"
#include <sstream>

#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
#include "cublas_v2.h"
#endif

#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
#include "cusparse.h"
#endif

#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER)
#include "cusolver_common.h"
#endif

namespace KokkosKernels {

#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
inline std::string cublas_version_string() {
  // Print version
  std::stringstream ss;

  ss << CUBLAS_VER_MAJOR << "." << CUBLAS_VER_MINOR << "." << CUBLAS_VER_PATCH;

  return ss.str();
}
#endif

#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
inline std::string cusparse_version_string() {
  // Print version
  std::stringstream ss;

  ss << CUSPARSE_VER_MAJOR << "." << CUSPARSE_VER_MINOR << "." << CUSPARSE_VER_PATCH << "." << CUSPARSE_VER_BUILD;

  return ss.str();
}
#endif

#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER)
inline std::string cusolver_version_string() {
  std::stringstream ss;

  ss << CUSOLVER_VER_MAJOR << "." << CUSOLVER_VER_MINOR << "." << CUSOLVER_VER_PATCH << "." << CUSOLVER_VER_BUILD;

  return ss.str();
}
#endif

}  // namespace KokkosKernels
#endif  // _KOKKOSKERNELS_TPLS_VERSIONS_HPP
