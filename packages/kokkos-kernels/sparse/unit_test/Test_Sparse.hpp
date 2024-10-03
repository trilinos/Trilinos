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
#ifndef TEST_SPARSE_HPP
#define TEST_SPARSE_HPP

#include "Test_Sparse_coo2crs.hpp"
#include "Test_Sparse_crs2coo.hpp"
#include "Test_Sparse_Controls.hpp"
#include "Test_Sparse_CrsMatrix.hpp"
#include "Test_Sparse_mdf.hpp"
#include "Test_Sparse_findRelOffset.hpp"
#include "Test_Sparse_gauss_seidel.hpp"
#include "Test_Sparse_MergeMatrix.hpp"
#include "Test_Sparse_replaceSumInto.hpp"
#include "Test_Sparse_replaceSumIntoLonger.hpp"
#include "Test_Sparse_spadd.hpp"
#include "Test_Sparse_spgemm_jacobi.hpp"
#include "Test_Sparse_spgemm.hpp"
#include "Test_Sparse_SortCrs.hpp"
#include "Test_Sparse_spiluk.hpp"
#include "Test_Sparse_spmv.hpp"
#include "Test_Sparse_sptrsv.hpp"
#include "Test_Sparse_trsv.hpp"
#include "Test_Sparse_par_ilut.hpp"
#include "Test_Sparse_gmres.hpp"
#include "Test_Sparse_Transpose.hpp"
#include "Test_Sparse_TestUtils_RandCsMat.hpp"
#include "Test_Sparse_IOUtils.hpp"
#include "Test_Sparse_ccs2crs.hpp"
#include "Test_Sparse_crs2ccs.hpp"
#include "Test_Sparse_removeCrsMatrixZeros.hpp"
#include "Test_Sparse_extractCrsDiagonalBlocks.hpp"

// TPL specific tests, these require
// particular pairs of backend and TPL
// to actually define tests.

#include "Test_Sparse_Utils_cusparse.hpp"
#include "Test_Sparse_rocsparse.hpp"

#endif  // TEST_SPARSE_HPP
