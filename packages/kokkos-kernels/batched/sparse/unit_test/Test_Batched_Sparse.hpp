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
#ifndef TEST_BATCHED_SPARSE_HPP
#define TEST_BATCHED_SPARSE_HPP

// Serial kernels
#include "Test_Batched_SerialGMRES.hpp"
#include "Test_Batched_SerialGMRES_Real.hpp"
#include "Test_Batched_SerialSpmv.hpp"
#include "Test_Batched_SerialSpmv_Real.hpp"

// Team Kernels
#include "Test_Batched_TeamCG.hpp"
#include "Test_Batched_TeamCG_Real.hpp"
#include "Test_Batched_TeamGMRES.hpp"
#include "Test_Batched_TeamGMRES_Real.hpp"
#include "Test_Batched_TeamSpmv.hpp"
#include "Test_Batched_TeamSpmv_Real.hpp"

// TeamVector Kernels
#include "Test_Batched_TeamVectorCG.hpp"
#include "Test_Batched_TeamVectorCG_Real.hpp"
#include "Test_Batched_TeamVectorGMRES.hpp"
#include "Test_Batched_TeamVectorGMRES_Real.hpp"
#include "Test_Batched_TeamVectorSpmv.hpp"
#include "Test_Batched_TeamVectorSpmv_Real.hpp"

// Vector Kernels

#endif  // TEST_BATCHED_SPARSE_HPP
