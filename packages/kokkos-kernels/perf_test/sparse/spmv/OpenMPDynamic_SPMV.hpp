/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef OPENMP_DYNAMIC_SPMV_HPP_
#define OPENMP_DYNAMIC_SPMV_HPP_

template <typename AType, typename XType, typename YType, typename Offset,
          typename Ordinal, typename Scalar>
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
