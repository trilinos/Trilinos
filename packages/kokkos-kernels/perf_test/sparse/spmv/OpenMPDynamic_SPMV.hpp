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

#ifndef OPENMP_DYNAMIC_SPMV_HPP_
#define OPENMP_DYNAMIC_SPMV_HPP_

template <typename AType, typename XType, typename YType, typename Offset, typename Ordinal, typename Scalar>
void openmp_dynamic_matvec(AType A, XType x, YType y) {
#define OMP_BENCH_RESTRICT __restrict__

  const Scalar s_a(1.0);
  const Scalar s_b(0.0);

  const Ordinal rowCount                            = A.numRows();
  const Scalar* OMP_BENCH_RESTRICT x_ptr            = x.data();
  Scalar* OMP_BENCH_RESTRICT y_ptr                  = y.data();
  const Scalar* OMP_BENCH_RESTRICT matrixCoeffs     = A.values.data();
  const Ordinal* OMP_BENCH_RESTRICT matrixCols      = A.graph.entries.data();
  const Offset* OMP_BENCH_RESTRICT matrixRowOffsets = A.graph.row_map.data();

#pragma omp parallel for schedule(dynamic)
  for (Ordinal row = 0; row < rowCount; ++row) {
    const Offset rowStart = matrixRowOffsets[row];
    const Offset rowEnd   = matrixRowOffsets[row + 1];

    Scalar sum(0.0);

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

#undef OMP_BENCH_RESTRICT
}

#endif /* OPENMP_DYNAMIC_SPMV_HPP_ */
