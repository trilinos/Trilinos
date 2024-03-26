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
#ifndef TEST_BLAS_HPP
#define TEST_BLAS_HPP

// Blas 1
#include "Test_Blas1_abs.hpp"
#include "Test_Blas1_asum.hpp"
#include "Test_Blas1_axpby.hpp"
#include "Test_Blas1_axpy.hpp"
#include "Test_Blas1_axpby_unification.hpp"
#include "Test_Blas1_dot.hpp"
#include "Test_Blas1_iamax.hpp"
#include "Test_Blas1_mult.hpp"
#include "Test_Blas1_nrm1.hpp"
#include "Test_Blas1_nrm2_squared.hpp"
#include "Test_Blas1_nrm2.hpp"
#include "Test_Blas1_nrm2w_squared.hpp"
#include "Test_Blas1_nrm2w.hpp"
#include "Test_Blas1_nrminf.hpp"
#include "Test_Blas1_reciprocal.hpp"
#include "Test_Blas1_rot.hpp"
#include "Test_Blas1_rotg.hpp"
#include "Test_Blas1_rotm.hpp"
#include "Test_Blas1_rotmg.hpp"
#include "Test_Blas1_scal.hpp"
#include "Test_Blas1_sum.hpp"
#include "Test_Blas1_swap.hpp"
#include "Test_Blas1_update.hpp"

// Serial Blas 1
#include "Test_Blas1_serial_setscal.hpp"
#include "Test_Blas_serial_axpy.hpp"
#include "Test_Blas_serial_nrm2.hpp"

// Team Blas 1
#include "Test_Blas1_team_setscal.hpp"
#include "Test_Blas1_team_abs.hpp"
#include "Test_Blas1_team_axpby.hpp"
#include "Test_Blas1_team_axpy.hpp"
#include "Test_Blas1_team_dot.hpp"
#include "Test_Blas1_team_mult.hpp"
#include "Test_Blas1_team_nrm2.hpp"
#include "Test_Blas1_team_scal.hpp"
#include "Test_Blas1_team_update.hpp"

// Blas 2
#include "Test_Blas2_gemv.hpp"
#include "Test_Blas2_ger.hpp"
#include "Test_Blas2_syr.hpp"
#include "Test_Blas2_syr2.hpp"

// Serial Blas 2
#include "Test_Blas2_serial_gemv.hpp"

// Team Blas 2
#include "Test_Blas2_team_gemv.hpp"
#include "Test_Blas2_teamvector_gemv.hpp"

// Blas 3
#include "Test_Blas3_gemm.hpp"
#include "Test_Blas3_trmm.hpp"
#include "Test_Blas3_trsm.hpp"

// TPLs
#include "Test_Blas_rocblas.hpp"

#endif  // TEST_BLAS_HPP
