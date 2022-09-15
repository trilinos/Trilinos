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

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

#include "KokkosSparse_Utils.hpp"
#include "KokkosSparse_SortCrs.hpp"
#include "KokkosSparse_spgemm.hpp"
#include "KokkosSparse_BsrMatrix.hpp"
#include "KokkosSparse_IOUtils.hpp"

using namespace KokkosSparse;

namespace Test {

template <typename bsrMat_t>
int run_block_spgemm(const bsrMat_t A, const bsrMat_t B, bsrMat_t &C,
                     // parameters
                     KokkosSparse::SPGEMMAlgorithm spgemm_algorithm,
                     bool use_dynamic_scheduling = true,
                     size_t shmem_size           = 0) {
  typedef typename bsrMat_t::size_type size_type;
  typedef typename bsrMat_t::ordinal_type lno_t;
  typedef typename bsrMat_t::value_type scalar_t;
  typedef typename bsrMat_t::device_type device;
  typedef typename bsrMat_t::memory_space memory_space;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle<
      size_type, lno_t, scalar_t, typename device::execution_space,
      memory_space, memory_space>
      KernelHandle;

  KernelHandle kh;
  kh.set_team_work_size(16);
  kh.set_dynamic_scheduling(use_dynamic_scheduling);

  kh.create_spgemm_handle(spgemm_algorithm);

  if (shmem_size > 0) {
    kh.set_shmem_size(shmem_size);
  }
  KokkosSparse::block_spgemm_symbolic(kh, A, false, B, false, C);
  KokkosSparse::block_spgemm_numeric(kh, A, false, B, false, C);
  kh.destroy_spgemm_handle();

  return 0;
}

template <typename bsrMat_t>
bool is_same_block_matrix(bsrMat_t output_mat_actual,
                          bsrMat_t output_mat_reference) {
  using device         = typename bsrMat_t::device_type;
  using graph_t        = typename bsrMat_t::StaticCrsGraphType;
  using lno_view_t     = typename graph_t::row_map_type::non_const_type;
  using lno_nnz_view_t = typename graph_t::entries_type::non_const_type;
  using scalar_view_t  = typename bsrMat_t::values_type::non_const_type;

  size_t nrows_actual    = output_mat_actual.numRows();
  size_t nentries_actual = output_mat_actual.graph.entries.extent(0);
  size_t nvals_actual    = output_mat_actual.values.extent(0);

  size_t nrows_reference    = output_mat_reference.numRows();
  size_t nentries_reference = output_mat_reference.graph.entries.extent(0);
  size_t nvals_reference    = output_mat_reference.values.extent(0);

  if (nrows_actual != nrows_reference) {
    std::cout << "nrows_actual:" << nrows_actual
              << " nrows_reference:" << nrows_reference << std::endl;
    return false;
  }
  if (nentries_actual != nentries_reference) {
    std::cout << "nentries_actual:" << nentries_actual
              << " nentries_reference:" << nentries_reference << std::endl;
    return false;
  }
  if (nvals_actual != nvals_reference) {
    std::cout << "nvals_actual:" << nvals_actual
              << " nvals_reference:" << nvals_reference << std::endl;
    return false;
  }

  KokkosSparse::sort_bsr_matrix(output_mat_actual);
  KokkosSparse::sort_bsr_matrix(output_mat_reference);

  bool is_identical = true;
  is_identical      = KokkosKernels::Impl::kk_is_identical_view<
      typename graph_t::row_map_type, typename graph_t::row_map_type,
      typename lno_view_t::value_type, typename device::execution_space>(
      output_mat_actual.graph.row_map, output_mat_reference.graph.row_map, 0);

  if (!is_identical) {
    std::cout << "rowmaps are different." << std::endl;
    std::cout << "Actual rowmap:\n";
    KokkosKernels::Impl::kk_print_1Dview(output_mat_actual.graph.row_map);
    std::cout << "Correct rowmap (SPGEMM_DEBUG):\n";
    KokkosKernels::Impl::kk_print_1Dview(output_mat_reference.graph.row_map);
    return false;
  }

  is_identical = KokkosKernels::Impl::kk_is_identical_view<
      lno_nnz_view_t, lno_nnz_view_t, typename lno_nnz_view_t::value_type,
      typename device::execution_space>(output_mat_actual.graph.entries,
                                        output_mat_reference.graph.entries, 0);

  if (!is_identical) {
    std::cout << "entries are different." << std::endl;
    KokkosKernels::Impl::kk_print_1Dview(output_mat_actual.graph.entries);
    KokkosKernels::Impl::kk_print_1Dview(output_mat_reference.graph.entries);
    return false;
  }

  typedef typename Kokkos::Details::ArithTraits<
      typename scalar_view_t::non_const_value_type>::mag_type eps_type;
  eps_type eps = std::is_same<eps_type, float>::value ? 3.7e-3 : 1e-7;

  is_identical = KokkosKernels::Impl::kk_is_relatively_identical_view<
      scalar_view_t, scalar_view_t, eps_type, typename device::execution_space>(
      output_mat_actual.values, output_mat_reference.values, eps);

  if (!is_identical) {
    std::cout << "values are different." << std::endl;
    KokkosKernels::Impl::kk_print_1Dview(output_mat_actual.values);
    KokkosKernels::Impl::kk_print_1Dview(output_mat_reference.values);

    return false;
  }
  return true;
}
}  // namespace Test

// Generate matrices and test all supported spgemm algorithms.
// C := AB, where A is m*k, B is k*n, and C is m*n.
template <typename scalar_t, typename lno_t, typename size_type,
          typename device>
void test_bspgemm(lno_t blkDim, lno_t m, lno_t k, lno_t n, size_type nnz,
                  lno_t bandwidth, lno_t row_size_variance,
                  const bool use_dynamic_scheduling = true,
                  const size_t shared_memory_size   = 0) {
  using namespace Test;
  // device::execution_space::initialize();
  // device::execution_space::print_configuration(std::cout);

  using bsrMat_t =
      KokkosSparse::Experimental::BsrMatrix<scalar_t, lno_t, device, void,
                                            size_type>;

  // Generate random compressed sparse row matrix. Randomly generated (non-zero)
  // values are stored in a 1-D (1 rank) array.
  bsrMat_t A = KokkosSparse::Impl::kk_generate_sparse_matrix<bsrMat_t>(
      blkDim, m, k, nnz, row_size_variance, bandwidth);
  bsrMat_t B = KokkosSparse::Impl::kk_generate_sparse_matrix<bsrMat_t>(
      blkDim, k, n, nnz, row_size_variance, bandwidth);

  const bool is_empy_case = m < 1 || n < 1 || k < 1 || nnz < 1;

  bsrMat_t output_mat2;
  run_block_spgemm(A, B, output_mat2, SPGEMM_DEBUG, use_dynamic_scheduling,
                   shared_memory_size);

  std::vector<SPGEMMAlgorithm> algorithms = {
      SPGEMM_KK,
      SPGEMM_KK_MEMORY /* alias SPGEMM_KK_MEMSPEED */,
      SPGEMM_KK_SPEED /* alias SPGEMM_KK_DENSE */,
      SPGEMM_MKL /* verify failure in case of missing build */,
  };

  if (!KokkosKernels::Impl::kk_is_gpu_exec_space<
          typename device::execution_space>()) {
    // SPGEMM_KK_LP is useful on CPU to cover MultiCoreTag4 functor
    // (otherwise skipped) but on GPU it's same as SPGEMM_KK, so we can skip it.
    algorithms.push_back(SPGEMM_KK_LP);
  }

  for (auto spgemm_algorithm : algorithms) {
    const uint64_t max_integer = Kokkos::ArithTraits<int>::max();
    std::string algo           = "UNKNOWN";
    bool is_expected_to_fail   = false;

    switch (spgemm_algorithm) {
      case SPGEMM_CUSPARSE:
        // TODO: add these test failure cases for cusparse too.
        algo = "SPGEMM_CUSPARSE";
#ifndef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
        is_expected_to_fail = true;
#endif
        break;

      case SPGEMM_MKL:
        algo                = "SPGEMM_MKL";
        is_expected_to_fail = !is_empy_case;  // TODO: add block MKL impl
#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
        if (!KokkosSparse::Impl::mkl_is_supported_value_type<scalar_t>::value) {
          is_expected_to_fail = true;
        }
#else
        is_expected_to_fail = true;  // fail: MKL not enabled in build
#endif
        // MKL requires local ordinals to be int.
        // Note: empty-array special case will NOT fail on this.
        if (!std::is_same<int, lno_t>::value && !is_empy_case) {
          is_expected_to_fail = true;
        }
        // if size_type is larger than int, mkl casts it to int.
        // it will fail if casting cause overflow.
        if (A.values.extent(0) > max_integer) {
          is_expected_to_fail = true;
        }
        break;

      case SPGEMM_KK: algo = "SPGEMM_KK"; break;
      case SPGEMM_KK_LP: algo = "SPGEMM_KK_LP"; break;
      case SPGEMM_KK_MEMSPEED: algo = "SPGEMM_KK_MEMSPEED"; break;
      case SPGEMM_KK_SPEED: algo = "SPGEMM_KK_SPEED"; break;
      case SPGEMM_KK_MEMORY: algo = "SPGEMM_KK_MEMORY"; break;
      default: algo = "!!! UNKNOWN ALGO !!!";
    }

    Kokkos::Timer timer1;
    bsrMat_t output_mat;

    bool failed = false;
    int res     = 0;
    try {
      res = run_block_spgemm(A, B, output_mat, spgemm_algorithm,
                             use_dynamic_scheduling, shared_memory_size);
    } catch (const char *message) {
      EXPECT_TRUE(is_expected_to_fail) << algo << ": " << message;
      failed = true;
    } catch (std::string message) {
      EXPECT_TRUE(is_expected_to_fail) << algo << ": " << message;
      failed = true;
    } catch (std::exception &e) {
      EXPECT_TRUE(is_expected_to_fail) << algo << ": " << e.what();
      failed = true;
    }
    EXPECT_EQ(is_expected_to_fail, failed);

    // double spgemm_time = timer1.seconds();

    timer1.reset();
    if (!is_expected_to_fail) {
      EXPECT_TRUE((res == 0)) << algo;
      bool is_identical = is_same_block_matrix(output_mat, output_mat2);
      EXPECT_TRUE(is_identical) << algo;
      // EXPECT_TRUE( equal) << algo;
    }
    // std::cout << "algo:" << algo << " spgemm_time:" << spgemm_time << "
    // output_check_time:" << timer1.seconds() << std::endl;
  }
  // device::execution_space::finalize();
}

// Note: Tests with shared memory specified aim to trigger specific GPU functors
//       dispatched by matrix size and the available shared memory.
#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)        \
  TEST_F(TestCategory,                                                     \
         sparse_block_spgemm_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    auto const SHMEM_AUTO = 0;                                             \
    auto test_case        = test_bspgemm<SCALAR, ORDINAL, OFFSET, DEVICE>; \
    /* Trigger SPGEMM_KK_MEMORY_SPREADTEAM on GPU */                       \
    test_case(2, 50, 50, 50, 2000, 50, 5, true, 16 * 1024);                \
    /* Trigger SPGEMM_KK -> SPGEMM_KK_MEMORY on GPU */                     \
    test_case(2, 50, 50, 50, 1000, 50, 5, false, 16 * 1024);               \
    /* Trigger SPGEMM_KK_MEMORY_BIGSPREADTEAM on GPU */                    \
    test_case(2, 500, 500, 500, 32000, 500, 500, true, 16 * 1024);         \
    /* trigger dense dispatch in hash method */                            \
    test_case(2, 2, 3, 4, 2, 2, 0, true, 16 * 1024);                       \
    /* zero-size handling */                                               \
    test_case(2, 0, 0, 0, 0, 10, 10, true, SHMEM_AUTO);                    \
    test_case(2, 0, 12, 5, 0, 10, 0, true, SHMEM_AUTO);                    \
    test_case(2, 10, 10, 0, 0, 10, 10, true, SHMEM_AUTO);                  \
  }

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST
