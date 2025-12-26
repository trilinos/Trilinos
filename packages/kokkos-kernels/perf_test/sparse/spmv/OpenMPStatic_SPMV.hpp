// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef OPENMP_STATIC_SPMV_HPP_
#define OPENMP_STATIC_SPMV_HPP_

template <typename AType, typename XType, typename YType, typename Offset, typename Ordinal, typename Scalar>
void openmp_static_matvec(AType A, XType x, YType y) {
#define OMP_BENCH_RESTRICT __restrict__

  const Scalar s_a(1.0);
  const Scalar s_b(0.0);

  const Ordinal rowCount                            = A.numRows();
  const Scalar* OMP_BENCH_RESTRICT x_ptr            = x.data();
  Scalar* OMP_BENCH_RESTRICT y_ptr                  = y.data();
  const Scalar* OMP_BENCH_RESTRICT matrixCoeffs     = A.values.data();
  const Ordinal* OMP_BENCH_RESTRICT matrixCols      = A.graph.entries.data();
  const Offset* OMP_BENCH_RESTRICT matrixRowOffsets = &A.graph.row_map(0);

#pragma omp parallel for
  for (Ordinal row = 0; row < rowCount; ++row) {
    const Offset rowStart = matrixRowOffsets[row];
    const Offset rowEnd   = matrixRowOffsets[row + 1];

    Scalar sum = 0.0;

    for (Offset i = rowStart; i < rowEnd; ++i) {
      const Ordinal x_entry = matrixCols[i];
      const Scalar alpha_MC = s_a * matrixCoeffs[i];
      sum += alpha_MC * x_ptr[x_entry];
    }

    if (0.0 == s_b)
      y_ptr[row] = sum;
    else
      y_ptr[row] = s_b * y_ptr[row] + sum;
  }

#undef OMP_BENCH_RESTRICT
}

#endif /* OPENMP_STATIC_SPMV_HPP_ */
