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

#include <cstdio>

#include <ctime>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <limits.h>
#include <cmath>
#include <unordered_map>
#include <KokkosSparse_spmv_test.hpp>
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

#ifdef KOKKOSKERNELS_ENABLE_TPL_ARMPL
#include <spmv/ArmPL_SPMV.hpp>
#endif

// return std::make_tuple(newnumRows, newnumCols, A, x1, y1,
//    rows_per_thread, team_size, vector_length,
//    test, schedule, ave_time, max_time, min_time);

SPMVTestData setup_test(spmv_additional_data* data, SPMVTestData::matrix_type A, Ordinal rows_per_thread, int team_size,
                        int vector_length, int schedule, int) {
  SPMVTestData test_data;
  using mv_type         = SPMVTestData::mv_type;
  using h_graph_type    = SPMVTestData::h_graph_type;
  using h_values_type   = SPMVTestData::h_values_type;
  test_data.A           = A;
  test_data.numRows     = A.numRows();
  test_data.numCols     = A.numCols();
  test_data.num_errors  = 0;
  test_data.total_error = 0;
  test_data.nnz         = A.nnz();
  mv_type x("X", test_data.numCols);
  mv_type y("Y", test_data.numRows);
  test_data.h_x         = Kokkos::create_mirror_view(x);
  test_data.h_y         = Kokkos::create_mirror_view(y);
  test_data.h_y_compare = Kokkos::create_mirror(y);

  h_graph_type h_graph   = Kokkos::create_mirror(test_data.A.graph);
  h_values_type h_values = Kokkos::create_mirror_view(test_data.A.values);

  for (int i = 0; i < test_data.numCols; i++) {
    test_data.h_x(i) = (Scalar)(1.0 * (rand() % 40) - 20.);
  }
  for (int i = 0; i < test_data.numRows; i++) {
    test_data.h_y(i) = (Scalar)(1.0 * (rand() % 40) - 20.);
  }

  test_data.generate_gold_standard(h_graph, h_values);

  Kokkos::deep_copy(x, test_data.h_x);
  Kokkos::deep_copy(y, test_data.h_y);
  Kokkos::deep_copy(test_data.A.graph.entries, h_graph.entries);
  Kokkos::deep_copy(test_data.A.values, h_values);
  test_data.x1 = mv_type("X1", test_data.numCols);
  Kokkos::deep_copy(test_data.x1, test_data.h_x);
  test_data.y1 = mv_type("Y1", test_data.numRows);

  // int nnz_per_row = A.nnz()/A.numRows(); // TODO: relocate
  matvec(A, test_data.x1, test_data.y1, rows_per_thread, team_size, vector_length, data, schedule);

  // Error Check
  Kokkos::deep_copy(test_data.h_y, test_data.y1);

  test_data.check_errors();
  test_data.min_time        = 1.0e32;
  test_data.ave_time        = 0.0;
  test_data.max_time        = 0.0;
  test_data.rows_per_thread = rows_per_thread;
  test_data.team_size       = team_size;
  test_data.vector_length   = vector_length;
  test_data.data            = data;
  test_data.schedule        = schedule;

  return test_data;
}
void run_benchmark(SPMVTestData& data) {
  Kokkos::Timer timer;
  matvec(data.A, data.x1, data.y1, data.rows_per_thread, data.team_size, data.vector_length, data.data, data.schedule);
  Kokkos::fence();
  double time = timer.seconds();
  data.ave_time += time;
  if (time > data.max_time) data.max_time = time;
  if (time < data.min_time) data.min_time = time;
}

struct SPMVConfiguration {
  int test;
  Ordinal rows_per_thread;
  Ordinal team_size;
  Ordinal vector_length;
  int schedule;
  int loop;
};

#ifdef KOKKOSKERNELS_ENABLE_TESTS_AND_PERFSUITE

namespace readers {
template <>
struct test_reader<SPMVConfiguration> {
  static SPMVConfiguration read(const std::string& filename) {
    std::ifstream input(filename);
    SPMVConfiguration config;
    input >> config.test >> config.rows_per_thread >> config.team_size >> config.vector_length >> config.schedule >>
        config.loop;
    return config;
  }
};

}  // namespace readers
test_list construct_kernel_base(const rajaperf::RunParams& run_params) {
  using matrix_type = SPMVTestData::matrix_type;
  srand(17312837);
  data_retriever<matrix_type, SPMVConfiguration> reader("sparse/spmv/", "sample.mtx", "config.cfg");
  std::vector<rajaperf::KernelBase*> test_cases;
  for (auto test_case : reader.test_cases) {
    auto& config = std::get<1>(test_case.test_data);
    test_cases.push_back(rajaperf::make_kernel_base(
        "Sparse_SPMV:" + test_case.filename, run_params,
        [=](const int, const int) {
          spmv_additional_data data(config.test);
          return std::make_tuple(setup_test(&data, std::get<0>(test_case.test_data), config.rows_per_thread,
                                            config.team_size, config.vector_length, config.schedule, config.loop));
        },
        [&](const int, const int, SPMVTestData& data) { run_benchmark(data); }));
  }
  return test_cases;
}

std::vector<rajaperf::KernelBase*> make_spmv_kernel_base(const rajaperf::RunParams& params) {
  return construct_kernel_base(params);
}

#endif  // KOKKOSKERNELS_ENABLE_TESTS_AND_PERFSUITE
