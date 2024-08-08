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

#ifndef OPENMP_SMART_STATIC_SPMV_HPP_
#define OPENMP_SMART_STATIC_SPMV_HPP_

#ifdef KOKKOS_ENABLE_OPENMP

#include <omp.h>

#define OMP_BENCH_RESTRICT __restrict__

// Done to prevent conflicting definitions of threadStarts
extern int* OMP_BENCH_RESTRICT threadStarts;

template <typename AType, typename Offset, typename Ordinal, typename Scalar>
void establishSmartSchedule(AType A) {
  const Ordinal rowCount                            = A.numRows();
  const Offset* OMP_BENCH_RESTRICT matrixRowOffsets = &A.graph.row_map(0);

  // Generate a schedule
  Ordinal* rowSizes = NULL;
  posix_memalign((void**)&rowSizes, 64, sizeof(int) * A.numRows());
  posix_memalign((void**)&threadStarts, 128, sizeof(int) * (omp_get_max_threads() + 1));

  for (int i = 0; i < omp_get_max_threads(); ++i) {
    threadStarts[i] = A.numRows();
  }

  unsigned long long int nnz = 0;

#pragma omp parallel for reduction(+ : nnz)
  for (Ordinal row = 0; row < rowCount; ++row) {
    const Ordinal rowElements = matrixRowOffsets[row + 1] - matrixRowOffsets[row];
    rowSizes[row]             = rowElements;
    nnz += rowElements;
  }

  Ordinal nzPerThreadTarget = (int)(nnz / (unsigned long long int)omp_get_max_threads());

  if (nzPerThreadTarget > 128) {
    nzPerThreadTarget &= 0xFFFFFFFC;
  }

  Ordinal nextRow = 0;

  printf("Target NZ Per Thread: %20d\n", nzPerThreadTarget);
  threadStarts[0] = 0;

  for (int thread = 1; thread < omp_get_max_threads(); ++thread) {
    Ordinal nzAccum = 0;

    while (nzAccum < nzPerThreadTarget) {
      if (nextRow >= rowCount) break;

      nzAccum += rowSizes[nextRow];
      nextRow++;
    }

    threadStarts[thread] = nextRow;
  }

  threadStarts[omp_get_max_threads()] = A.numRows();

  // printf("Schedule: Target-per-Row=%20d\n", rowsPerThreadTarget);
  // for(int i = 0; i < omp_get_max_threads(); ++i) {
  //	printf("thread [%5d] start=%20d end=%20d\n", i, threadStarts[i],
  // threadStarts[i+1]);
  //}

  free(rowSizes);
}

template <typename AType, typename XType, typename YType, typename Offset, typename Ordinal, typename Scalar>
void openmp_smart_static_matvec(AType A, XType x, YType y) {
  if (NULL == threadStarts) {
    // printf("Generating Schedule...\n");
    establishSmartSchedule<AType, Offset, Ordinal, Scalar>(A);
  }

  const Scalar s_a = 1.0;
  const Scalar s_b = 0.0;

  // const Ordinal rowCount                           = A.numRows();
  const Scalar* OMP_BENCH_RESTRICT x_ptr            = (Scalar*)x.data();
  Scalar* OMP_BENCH_RESTRICT y_ptr                  = (Scalar*)y.data();
  const Scalar* OMP_BENCH_RESTRICT matrixCoeffs     = A.values.data();
  const Ordinal* OMP_BENCH_RESTRICT matrixCols      = A.graph.entries.data();
  const Offset* OMP_BENCH_RESTRICT matrixRowOffsets = &A.graph.row_map(0);

#ifdef KOKKOS_ENABLE_PROFILING
  uint64_t kpID = 0;
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::beginParallelFor("KokkosSparse::Test_SPMV_raw_openmp", 0, &kpID);
  }
#endif

#pragma omp parallel
  {
#if defined(KOKKOS_COMPILER_INTEL) && !defined(__clang__)
    __assume_aligned(x_ptr, 64);
    __assume_aligned(y_ptr, 64);
#endif

    const int myID        = omp_get_thread_num();
    const Ordinal myStart = threadStarts[myID];
    const Ordinal myEnd   = threadStarts[myID + 1];

    for (int row = myStart; row < myEnd; ++row) {
      const Offset rowStart = matrixRowOffsets[row];
      const Offset rowEnd   = matrixRowOffsets[row + 1];

      Scalar sum = 0.0;

      for (Offset i = rowStart; i < rowEnd; ++i) {
        const Ordinal x_entry = matrixCols[i];
        const Scalar alpha_MC = s_a * matrixCoeffs[i];
        sum += alpha_MC * x_ptr[x_entry];
      }

      if (0.0 == s_b) {
        y_ptr[row] = sum;
      } else {
        y_ptr[row] = s_b * y_ptr[row] + sum;
      }
    }
  }
#ifdef KOKKOS_ENABLE_PROFILING
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::endParallelFor(kpID);
  }
#endif
}

#undef OMP_BENCH_RESTRICT

#endif /* KOKKOS_ENABLE_OPENMP */
#endif /* OPENMP_SMART_STATIC_SPMV_HPP_ */
