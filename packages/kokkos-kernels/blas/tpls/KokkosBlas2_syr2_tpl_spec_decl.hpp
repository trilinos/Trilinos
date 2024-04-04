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

#ifndef KOKKOSBLAS2_SYR2_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS2_SYR2_TPL_SPEC_DECL_HPP_

// BLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
#include <KokkosBlas2_syr2_tpl_spec_decl_blas.hpp>
#endif

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include <KokkosBlas2_syr2_tpl_spec_decl_cublas.hpp>
#endif

// rocBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
#include <KokkosBlas2_syr2_tpl_spec_decl_rocblas.hpp>
#endif

#endif
