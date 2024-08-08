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
namespace KokkosSparse {
namespace Impl {

#ifdef KOKKOS_ENABLE_OPENMP
template <typename AMatrix, typename XVector, typename YVector>
void spmv_raw_openmp_no_transpose(typename YVector::const_value_type& s_a, AMatrix A, XVector x,
                                  typename YVector::const_value_type& s_b, YVector y) {
  typedef typename YVector::non_const_value_type value_type;
  typedef typename AMatrix::ordinal_type ordinal_type;
  typedef typename AMatrix::non_const_size_type size_type;

  typename XVector::const_value_type* KOKKOS_RESTRICT x_ptr     = x.data();
  typename YVector::non_const_value_type* KOKKOS_RESTRICT y_ptr = y.data();

  const typename AMatrix::value_type* KOKKOS_RESTRICT matrixCoeffs = A.values.data();
  const ordinal_type* KOKKOS_RESTRICT matrixCols                   = A.graph.entries.data();
  const size_type* KOKKOS_RESTRICT matrixRowOffsets                = A.graph.row_map.data();
  const size_type* KOKKOS_RESTRICT threadStarts                    = A.graph.row_block_offsets.data();

#if defined(KOKKOS_ENABLE_PROFILING)
  uint64_t kpID = 0;
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::beginParallelFor("KokkosSparse::spmv<RawOpenMP,NoTranspose>", 0, &kpID);
  }
#endif

  typename YVector::const_value_type zero = 0;
#pragma omp parallel
  {
#if defined(KOKKOS_COMPILER_INTEL) && !defined(__clang__)
    __assume_aligned(x_ptr, 64);
    __assume_aligned(y_ptr, 64);
#endif

    const int myID          = omp_get_thread_num();
    const size_type myStart = threadStarts[myID];
    const size_type myEnd   = threadStarts[myID + 1];

    for (size_type row = myStart; row < myEnd; ++row) {
      const size_type rowStart = matrixRowOffsets[row];
      const size_type rowEnd   = matrixRowOffsets[row + 1];

      value_type sum = 0.0;

      for (size_type i = rowStart; i < rowEnd; ++i) {
        const ordinal_type x_entry = matrixCols[i];
        const value_type alpha_MC  = s_a * matrixCoeffs[i];
        sum += alpha_MC * x_ptr[x_entry];
      }

      if (zero == s_b) {
        y_ptr[row] = sum;
      } else {
        y_ptr[row] = s_b * y_ptr[row] + sum;
      }
    }
  }
#if defined(KOKKOS_ENABLE_PROFILING)
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::endParallelFor(kpID);
  }
#endif
}

#endif
}  // namespace Impl
}  // namespace KokkosSparse
