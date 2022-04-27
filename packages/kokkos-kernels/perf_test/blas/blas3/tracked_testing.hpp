//
// Created by Poliakoff, David Zoeller on 4/26/21.
//
#ifndef KOKKOSKERNELS_BLAS3_TRACKED_TESTING_HPP
#define KOKKOSKERNELS_BLAS3_TRACKED_TESTING_HPP

#include <common/RAJAPerfSuite.hpp>
#include <common/Executor.hpp>

#include "KokkosBlas3_gemm_tracked_perf_test.hpp"

/*
 *Three cases to test:
 *
 * 1) m = n = k
 * 2) one case for m, n, k all pretty large,
 *  3) and another for m, k small but n  large
 *
 *  You could use m = k = 5, n = 1 million or something like that for dot based
 *  gemm
 *
 */

namespace test {
namespace blas3 {

// Change n and k values in the context of the backend
template <class ExecSpace>
std::vector<gemmConfig> create_m_n_k_vect() {
  std::string exec_space_name = ExecSpace::name();

  return {

      // m = n = k; one case for m, n, k all pretty large,
      {1000, 1000, 1000},
      // and another for m, k small but n  large
      {5, 1000000, 5}};
}

// Register kernels for a specific test
void build_gemm_executor(rajaperf::Executor& exec, int argc, char* argv[],
                         const rajaperf::RunParams& params) {
  for (auto* kernel : construct_gemm_kernel_base(
           params, create_m_n_k_vect<Kokkos::DefaultExecutionSpace>())) {
    exec.registerKernel("BLAS3", kernel);
  }
}

void build_blas3_executor(rajaperf::Executor& exec, int argc, char* argv[],
                          const rajaperf::RunParams& params) {
  exec.registerGroup("BLAS3");
  build_gemm_executor(exec, argc, argv, params);
}

}  // namespace blas3
}  // namespace test

#endif  // KOKKOSKERNELS_TRACKED_TESTING_HPP
