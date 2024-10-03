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
//
// Created by Poliakoff, David Zoeller on 4/26/21.
//

#ifndef KOKKOSKERNELS_KOKKOSSPARSE_SPMV_TEST_HPP
#define KOKKOSKERNELS_KOKKOSSPARSE_SPMV_TEST_HPP

#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosKernels_IOUtils.hpp>
#include <KokkosSparse_spmv.hpp>
#include "KokkosKernels_default_types.hpp"
#include <spmv/Kokkos_SPMV.hpp>
#include <spmv/Kokkos_SPMV_Inspector.hpp>

#include <spmv/KokkosKernels_spmv_data.hpp>

#ifdef KOKKOSKERNELS_ENABLE_TESTS_AND_PERFSUITE
#include <PerfTestUtilities.hpp>
#endif

#ifdef KOKKOS_ENABLE_OPENMP
#include <spmv/OpenMPStatic_SPMV.hpp>
#include <spmv/OpenMPDynamic_SPMV.hpp>
#include <spmv/OpenMPSmartStatic_SPMV.hpp>
#endif

#ifdef HAVE_CUSPARSE
#include <spmv/CuSparse_SPMV.hpp>
#endif

#ifdef HAVE_MKL
#include <spmv/MKL_SPMV.hpp>
#endif

template <typename AType, typename XType, typename YType>
void armpl_matvec(AType /*A*/, XType x, YType y, spmv_additional_data* data);

enum { KOKKOS, MKL, ARMPL, CUSPARSE, KK_KERNELS, KK_KERNELS_INSP, KK_INSP, OMP_STATIC, OMP_DYNAMIC, OMP_INSP };
enum { AUTO, DYNAMIC, STATIC };

using Scalar  = default_scalar;
using Ordinal = default_lno_t;
using Offset  = default_size_type;
using Layout  = default_layout;

#ifdef KOKKOSKERNELS_ENABLE_TESTS_AND_PERFSUITE
std::vector<rajaperf::KernelBase*> make_spmv_kernel_base(const rajaperf::RunParams& params);

test_list construct_kernel_base(const rajaperf::RunParams& run_params, Ordinal numRows, Ordinal numCols,
                                spmv_additional_data* data, Ordinal rows_per_thread, int team_size, int vector_length,
                                int schedule, int loop);

#endif

struct SPMVTestData {
  using matrix_type   = KokkosSparse::CrsMatrix<Scalar, Ordinal, Kokkos::DefaultExecutionSpace, void, Offset>;
  using mv_type       = Kokkos::View<Scalar*, Layout>;
  using h_mv_type     = mv_type::HostMirror;
  using h_graph_type  = matrix_type::StaticCrsGraphType::HostMirror;
  using h_values_type = matrix_type::values_type::HostMirror;

  Ordinal numRows;
  Ordinal numCols;
  matrix_type A;
  mv_type x1;
  mv_type y1;
  Ordinal rows_per_thread;
  Ordinal team_size;
  Ordinal vector_length;
  Ordinal nnz;
  h_mv_type h_y_compare;
  h_mv_type h_x;
  h_mv_type h_y;
  // int test;
  spmv_additional_data* data;
  int schedule;
  int num_errors;
  int num_error_checks;
  double total_error;
  double ave_time;
  double max_time;
  double min_time;
  inline void check_errors() {
    double error = 0.0;
    double sum   = 0.0;

    for (int i = 0; i < numRows; i++) {
      error += (h_y_compare(i) - h_y(i)) * (h_y_compare(i) - h_y(i));
      sum += h_y_compare(i) * h_y_compare(i);
    }
    if (sum == 0.0) {
      sum = 1.0;  // protect against div by zero
    }
    if ((error / sum) > 1e-5) {  // if above tolerance
      ++num_errors;
    }
    total_error += error;
    ++num_error_checks;
  }
  inline void generate_gold_standard(h_graph_type h_graph, h_values_type h_values) {
    for (int i = 0; i < numRows; i++) {
      int start = h_graph.row_map(i);
      int end   = h_graph.row_map(i + 1);
      for (int j = start; j < end; j++) {
        h_values(j) = h_graph.entries(j) + i;
      }

      h_y_compare(i) = 0;
      for (int j = start; j < end; j++) {
        Scalar tmp_val = h_graph.entries(j) + i;
        int idx        = h_graph.entries(j);
        h_y_compare(i) += tmp_val * h_x(idx);
      }
    }
  }
};

SPMVTestData setup_test(spmv_additional_data* data, SPMVTestData::matrix_type A, Ordinal rows_per_thread, int team_size,
                        int vector_length, int schedule, int loop);

template <typename AType, typename XType, typename YType>
void matvec(AType& A, XType x, YType y, Ordinal rows_per_thread, int team_size, int vector_length,
            spmv_additional_data* data, int schedule) {
  switch (data->test) {
    case KOKKOS:
      if (schedule == AUTO) schedule = A.nnz() > 10000000 ? DYNAMIC : STATIC;
      if (schedule == STATIC)
        kokkos_matvec<AType, XType, YType, Kokkos::Static>(A, x, y, rows_per_thread, team_size, vector_length);
      if (schedule == DYNAMIC)
        kokkos_matvec<AType, XType, YType, Kokkos::Dynamic>(A, x, y, rows_per_thread, team_size, vector_length);
      break;
    case KK_INSP:
      if (schedule == AUTO) schedule = A.nnz() > 10000000 ? DYNAMIC : STATIC;
      if (schedule == STATIC)
        kk_inspector_matvec<AType, XType, YType, Kokkos::Static>(A, x, y, team_size, vector_length);
      if (schedule == DYNAMIC)
        kk_inspector_matvec<AType, XType, YType, Kokkos::Dynamic>(A, x, y, team_size, vector_length);
      break;

#ifdef KOKKOS_ENABLE_OPENMP
    case OMP_STATIC: openmp_static_matvec<AType, XType, YType, Offset, Ordinal, Scalar>(A, x, y); break;
    case OMP_DYNAMIC: openmp_dynamic_matvec<AType, XType, YType, Offset, Ordinal, Scalar>(A, x, y); break;
    case OMP_INSP: openmp_smart_static_matvec<AType, XType, YType, Offset, Ordinal, Scalar>(A, x, y); break;
#endif

#ifdef HAVE_MKL
    case MKL: mkl_matvec(A, x, y); break;
#endif
#ifdef HAVE_CUSPARSE
    case CUSPARSE: cusparse_matvec(A, x, y); break;
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_ARMPL
    case ARMPL: armpl_matvec(A, x, y, data); break;
#endif
    case KK_KERNELS: KokkosSparse::spmv(KokkosSparse::NoTranspose, 1.0, A, x, 0.0, y); break;
    case KK_KERNELS_INSP:
      if (A.graph.row_block_offsets.data() == NULL) {
        printf("PTR: %p\n", static_cast<const void*>(A.graph.row_block_offsets.data()));
        A.graph.create_block_partitioning(typename AType::execution_space().concurrency());
        printf("PTR2: %p\n", static_cast<const void*>(A.graph.row_block_offsets.data()));
      }
      KokkosSparse::spmv(KokkosSparse::NoTranspose, 1.0, A, x, 0.0, y);
      break;
    default: fprintf(stderr, "Selected test is not available.\n");
  }
}

void run_benchmark(SPMVTestData& data);

#endif  // KOKKOSKERNELS_KOKKOSSPARSE_SPMV_HPP
