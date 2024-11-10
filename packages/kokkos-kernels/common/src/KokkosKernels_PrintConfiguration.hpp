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

#ifndef _KOKKOSKERNELS_PRINT_CONFIGURATION_HPP
#define _KOKKOSKERNELS_PRINT_CONFIGURATION_HPP

#include "KokkosKernels_config.h"
#include "KokkosKernels_TplsVersion.hpp"
#include <iostream>

namespace KokkosKernels {
namespace Impl {

inline void print_cublas_version_if_enabled(std::ostream& os) {
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_CUBLAS: " << cublas_version_string() << "\n";
#else
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_CUBLAS: no\n";
#endif
}

inline void print_cusparse_version_if_enabled(std::ostream& os) {
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_CUSPARSE: " << cusparse_version_string() << "\n";
#else
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_CUSPARSE: no\n";
#endif
}

inline void print_cusolver_version_if_enabled(std::ostream& os) {
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSOLVER
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_CUSOLVER: " << cusolver_version_string() << "\n";
#else
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_CUSOLVER: no\n";
#endif
}

inline void print_enabled_tpls(std::ostream& os) {
#ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACK
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_LAPACK: yes\n";
#else
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_LAPACK: no\n";
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_BLAS: yes\n";
#else
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_BLAS: no\n";
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_CBLAS
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_CBLAS: yes\n";
#else
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_CBLAS: no\n";
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACKE
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_LAPACKE: yes\n";
#else
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_LAPACKE: no\n";
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_SUPERLU
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_SUPERLU: yes\n";
#else
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_SUPERLU: no\n";
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_CHOLMOD
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_CHOLMOD: yes\n";
#else
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_CHOLMOD: no\n";
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_MKL: yes\n";
#else
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_MKL: no\n";
#endif
  print_cublas_version_if_enabled(os);
  print_cusparse_version_if_enabled(os);
  print_cusolver_version_if_enabled(os);
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_ROCBLAS: yes\n";
#else
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_ROCBLAS: no\n";
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE: yes\n";
#else
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE: no\n";
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER: yes\n";
#else
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_ROCOLVER: no\n";
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_METIS
  os << "KOKKOSKERNELS_ENABLE_TPL_METIS: yes\n";
#else
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_METIS: no\n";
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_ARMPL
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_ARMPL: yes\n";
#else
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_ARMPL: no\n";
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_MAGMA: yes\n";
#else
  os << "  "
     << "KOKKOSKERNELS_ENABLE_TPL_MAGMA: no\n";
#endif
}

inline void print_version(std::ostream& os) {
  // KOKKOSKERNELS_VERSION is used because MAJOR, MINOR and PATCH macros
  // are not available in Kernels
  os << "  "
     << "KokkosKernels Version: " << KOKKOSKERNELS_VERSION_MAJOR << "." << KOKKOSKERNELS_VERSION_MINOR << "."
     << KOKKOSKERNELS_VERSION_PATCH << '\n';
}

}  // namespace Impl

inline void print_configuration(std::ostream& os) {
  Impl::print_version(os);

  os << "TPLs: \n";
  Impl::print_enabled_tpls(os);
}

}  // namespace KokkosKernels
#endif  // _KOKKOSKERNELS_PRINT_CONFIGURATION_HPP
