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
#ifndef TEST_BATCHED_DENSE_HPP
#define TEST_BATCHED_DENSE_HPP

// Serial kernels
#include "Test_Batched_SerialAxpy.hpp"
#include "Test_Batched_SerialAxpy_Real.hpp"
#include "Test_Batched_SerialAxpy_Complex.hpp"
#include "Test_Batched_SerialEigendecomposition.hpp"
#include "Test_Batched_SerialEigendecomposition_Real.hpp"
#include "Test_Batched_SerialGesv.hpp"
#include "Test_Batched_SerialGesv_Real.hpp"
#include "Test_Batched_SerialInverseLU.hpp"
#include "Test_Batched_SerialInverseLU_Real.hpp"
#include "Test_Batched_SerialInverseLU_Complex.hpp"
#include "Test_Batched_SerialLU.hpp"
#include "Test_Batched_SerialLU_Real.hpp"
#include "Test_Batched_SerialLU_Complex.hpp"
#include "Test_Batched_SerialSolveLU.hpp"
#include "Test_Batched_SerialSolveLU_Real.hpp"
#include "Test_Batched_SerialSolveLU_Complex.hpp"
#include "Test_Batched_SerialTrmm.hpp"
#include "Test_Batched_SerialTrmm_Real.hpp"
#include "Test_Batched_SerialTrmm_Complex.hpp"
#include "Test_Batched_SerialTrsm.hpp"
#include "Test_Batched_SerialTrsm_Real.hpp"
#include "Test_Batched_SerialTrsm_Complex.hpp"
#include "Test_Batched_SerialTrsv.hpp"
#include "Test_Batched_SerialTrsv_Real.hpp"
#include "Test_Batched_SerialTrsv_Complex.hpp"
#include "Test_Batched_SerialTbsv.hpp"
#include "Test_Batched_SerialTbsv_Real.hpp"
#include "Test_Batched_SerialTbsv_Complex.hpp"
#include "Test_Batched_SerialTrtri.hpp"
#include "Test_Batched_SerialTrtri_Real.hpp"
#include "Test_Batched_SerialTrtri_Complex.hpp"
#include "Test_Batched_SerialSVD.hpp"
#include "Test_Batched_SerialPttrf.hpp"
#include "Test_Batched_SerialPttrf_Real.hpp"
#include "Test_Batched_SerialPttrf_Complex.hpp"

// Team Kernels
#include "Test_Batched_TeamAxpy.hpp"
#include "Test_Batched_TeamAxpy_Real.hpp"
#include "Test_Batched_TeamAxpy_Complex.hpp"
#include "Test_Batched_TeamGesv.hpp"
#include "Test_Batched_TeamGesv_Real.hpp"
#include "Test_Batched_TeamInverseLU.hpp"
#include "Test_Batched_TeamInverseLU_Real.hpp"
#include "Test_Batched_TeamInverseLU_Complex.hpp"
#include "Test_Batched_TeamLU.hpp"
#include "Test_Batched_TeamLU_Real.hpp"
#include "Test_Batched_TeamLU_Complex.hpp"
#include "Test_Batched_TeamSolveLU.hpp"
#include "Test_Batched_TeamSolveLU_Real.hpp"
#include "Test_Batched_TeamSolveLU_Complex.hpp"
#include "Test_Batched_TeamTrsm.hpp"
#include "Test_Batched_TeamTrsm_Real.hpp"
#include "Test_Batched_TeamTrsm_Complex.hpp"
#include "Test_Batched_TeamTrsv.hpp"
#include "Test_Batched_TeamTrsv_Real.hpp"
#include "Test_Batched_TeamTrsv_Complex.hpp"

// TeamVector Kernels
#include "Test_Batched_TeamVectorAxpy.hpp"
#include "Test_Batched_TeamVectorAxpy_Real.hpp"
#include "Test_Batched_TeamVectorAxpy_Complex.hpp"
#include "Test_Batched_TeamVectorEigendecomposition.hpp"
#include "Test_Batched_TeamVectorEigendecomposition_Real.hpp"
#include "Test_Batched_TeamVectorGesv.hpp"
#include "Test_Batched_TeamVectorGesv_Real.hpp"
#include "Test_Batched_TeamVectorQR.hpp"
#include "Test_Batched_TeamVectorQR_Real.hpp"
#include "Test_Batched_TeamVectorQR_WithColumnPivoting.hpp"
#include "Test_Batched_TeamVectorQR_WithColumnPivoting_Real.hpp"
#include "Test_Batched_TeamVectorSolveUTV.hpp"
#include "Test_Batched_TeamVectorSolveUTV_Real.hpp"
#include "Test_Batched_TeamVectorSolveUTV2.hpp"
#include "Test_Batched_TeamVectorSolveUTV2_Real.hpp"
#include "Test_Batched_TeamVectorUTV.hpp"
#include "Test_Batched_TeamVectorUTV_Real.hpp"

// Vector Kernels
#include "Test_Batched_VectorArithmatic.hpp"
#include "Test_Batched_VectorLogical.hpp"
#include "Test_Batched_VectorMath.hpp"
#include "Test_Batched_VectorMisc.hpp"
#include "Test_Batched_VectorRelation.hpp"
#include "Test_Batched_VectorView.hpp"

#endif  // TEST_BATCHED_DENSE_HPP
