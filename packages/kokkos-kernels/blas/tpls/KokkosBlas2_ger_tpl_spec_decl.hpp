// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS2_GER_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS2_GER_TPL_SPEC_DECL_HPP_

// BLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
#include <KokkosBlas2_ger_tpl_spec_decl_blas.hpp>
#endif

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include <KokkosBlas2_ger_tpl_spec_decl_cublas.hpp>
#endif

// rocBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
#include <KokkosBlas2_ger_tpl_spec_decl_rocblas.hpp>
#endif

#endif
