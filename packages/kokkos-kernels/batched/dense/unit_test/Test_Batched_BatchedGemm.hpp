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
#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

#include "KokkosBatched_HostLevel_Gemm.hpp"
#include "KokkosBatched_HostLevel_Gemm_DblBuf_Impl.hpp"

#include "KokkosKernels_TestUtils.hpp"

using namespace KokkosBatched;

namespace Test {
template <typename DeviceType, typename ViewType, typename ScalarType, typename ParamTagType>
void impl_test_batched_gemm_with_handle(BatchedGemmHandle* batchedGemmHandle, const int N, const int matAdim1,
                                        const int matAdim2, const int matBdim1, const int matBdim2, const int matCdim1,
                                        const int matCdim2, ScalarType alpha, ScalarType beta) {
  using execution_space = typename DeviceType::execution_space;
  using transA          = typename ParamTagType::transA;
  using transB          = typename ParamTagType::transB;
  using batchLayout     = typename ParamTagType::batchLayout;
  using ats             = Kokkos::ArithTraits<ScalarType>;

  int ret        = 0;
  auto algo_type = batchedGemmHandle->get_kernel_algo_type();
  ViewType a_expected, a_actual, b_expected, b_actual, c_expected, c_actual;
  std::string fmsg;
  std::string fmsg_rhs = "algo_type:" + batchedGemmHandle->get_kernel_algo_type_str() + ", ";
  fmsg_rhs += ("N:" + std::to_string(N) + ", ");
  fmsg_rhs += ("A:" + std::to_string(matAdim1) + "x" + std::to_string(matAdim2) + ", ");
  fmsg_rhs += ("B:" + std::to_string(matBdim1) + "x" + std::to_string(matBdim2) + ", ");
  fmsg_rhs += ("C:" + std::to_string(matCdim1) + "x" + std::to_string(matCdim2) + "\n");

  if (std::is_same<batchLayout, BatchLayout::Left>::value) {
    a_expected = ViewType("a_expected", N, matAdim1, matAdim2);
    a_actual   = ViewType("a_actual", N, matAdim1, matAdim2);
    b_expected = ViewType("b_expected", N, matBdim1, matBdim2);
    b_actual   = ViewType("b_actual", N, matBdim1, matBdim2);
    c_expected = ViewType("c_expected", N, matCdim1, matCdim2);
    c_actual   = ViewType("c_actual", N, matCdim1, matCdim2);
  } else if (std::is_same<batchLayout, BatchLayout::Right>::value) {
    a_expected = ViewType("a_expected", matAdim1, matAdim2, N);
    a_actual   = ViewType("a_actual", matAdim1, matAdim2, N);
    b_expected = ViewType("b_expected", matBdim1, matBdim2, N);
    b_actual   = ViewType("b_actual", matBdim1, matBdim2, N);
    c_expected = ViewType("c_expected", matCdim1, matCdim2, N);
    c_actual   = ViewType("c_actual", matCdim1, matCdim2, N);
  }

  Kokkos::Random_XorShift64_Pool<execution_space> random(13718);

  Kokkos::fill_random(a_expected, random, ScalarType(1.0));
  Kokkos::fill_random(b_expected, random, ScalarType(1.0));
  Kokkos::fill_random(c_expected, random, ScalarType(1.0));

  Kokkos::fence();

  Kokkos::deep_copy(a_actual, a_expected);
  Kokkos::deep_copy(b_actual, b_expected);
  Kokkos::deep_copy(c_actual, c_expected);

  // Check for expected BatchedDblBufGemm runtime errors
  if (algo_type == GemmKokkosBatchedAlgos::KK_DBLBUF) {
    // Check for DblBuf runtime errors related to team_size
    try {
      fmsg = kk_failure_str(__FILE__, __FUNCTION__, __LINE__);
      Impl::BatchedDblBufGemm<transA, transB, batchLayout, BatchedGemmHandle, ScalarType, decltype(a_actual),
                              decltype(b_actual), decltype(c_actual), BoundsCheck::Yes, AlphaTag::No, 65536, 1, 65536>(
          batchedGemmHandle, alpha, a_actual, b_actual, beta, c_actual)
          .invoke();
      FAIL() << (fmsg + fmsg_rhs);
    } catch (const std::runtime_error& error) {
      ;
    }

    // Check for DblBuf runtime errors related to vector_len
    try {
      fmsg = kk_failure_str(__FILE__, __FUNCTION__, __LINE__);
      Impl::BatchedDblBufGemm<transA, transB, batchLayout, BatchedGemmHandle, ScalarType, decltype(a_actual),
                              decltype(b_actual), decltype(c_actual), BoundsCheck::No, AlphaTag::No, 65536, 65536 * 2,
                              65536>(batchedGemmHandle, alpha, a_actual, b_actual, beta, c_actual)
          .invoke();
      FAIL() << (fmsg + fmsg_rhs);
    } catch (const std::runtime_error& error) {
      ;
    }
  }

  // Check for expected BatchedGemm runtime errors
  try {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ARMPL)
    if (algo_type == BaseTplAlgos::ARMPL && N % 2 == 0) {
      auto ninter = batchedGemmHandle->get_tpl_params();
      // Set ninter parameter for underlying armpl_dgemm_interleave_batch call
      *ninter = N / 2;
    }
#endif

    fmsg = kk_failure_str(__FILE__, __FUNCTION__, __LINE__);
    ret  = BatchedGemm<transA, transB, batchLayout>(batchedGemmHandle, alpha, a_actual, b_actual, beta,
                                                   c_actual);  // Compute c_actual
  } catch (const std::runtime_error& error) {
    std::string error_msg = error.what();
    if (algo_type == BaseHeuristicAlgos::SQUARE && matCdim1 != matCdim2) {
      ;
    } else if (algo_type == BaseTplAlgos::ARMPL) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ARMPL) && ARMPL_BUILD >= 1058
      auto ninter = batchedGemmHandle->get_tpl_params()[0];
      // No runtime errors expected since layout is valid, double is a supported
      // type, and ninter != 0
      if (std::is_same<typename ViewType::value_type, double>::value && ninter != 0) {
        FAIL() << (error_msg + fmsg + fmsg_rhs);
      }
#else
      ;  // We expect a runtime error if the ARMPL TPL is not enabled
#endif
    } else {
      FAIL() << (error_msg + fmsg + fmsg_rhs);
    }
    return;
  }
  ASSERT_EQ(ret, 0) << (fmsg + fmsg_rhs);

  Functor_BatchedVanillaGEMM<ViewType, ViewType, ViewType, execution_space> vgemm;
  vgemm.A_t                 = std::is_same<transA, Trans::Transpose>::value;
  vgemm.B_t                 = std::is_same<transB, Trans::Transpose>::value;
  vgemm.batch_size_last_dim = std::is_same<batchLayout, BatchLayout::Right>::value;
  vgemm.A_c = vgemm.B_c = false;
  vgemm.A               = a_expected;
  vgemm.B               = b_expected;
  vgemm.C               = c_expected;
  vgemm.alpha           = alpha;
  vgemm.beta            = beta;
  vgemm.run();  // Compute c_expected

  Kokkos::fence();

  typename ViewType::HostMirror c_expected_host = Kokkos::create_mirror_view(c_expected);
  typename ViewType::HostMirror c_actual_host   = Kokkos::create_mirror_view(c_actual);

  // Copy to host
  Kokkos::deep_copy(c_expected_host, c_expected);
  Kokkos::deep_copy(c_actual_host, c_actual);

  Kokkos::fence();

  // check c_expected = c_actual ; this eps is about 2^-9
  // Set mag_type to host_value_type, we may not have half precision on host
  using mag_type = float;
  mag_type sum(1), diff(0);

  auto eps = static_cast<mag_type>(ats::epsilon());

  eps *= std::is_same<ScalarType, Kokkos::Experimental::half_t>::value ||
                 std::is_same<ScalarType, Kokkos::Experimental::bhalf_t>::value
             ? 4
             : 1e3;

  for (int k = 0; k < N; ++k) {
    for (int i = 0; i < matCdim1; ++i) {
      for (int j = 0; j < matCdim2; ++j) {
        if (std::is_same<batchLayout, BatchLayout::Right>::value) {
          sum += ats::abs(c_expected_host(i, j, k));
          diff += ats::abs(c_expected_host(i, j, k) - c_actual_host(i, j, k));
        } else {
          sum += ats::abs(c_expected_host(k, i, j));
          diff += ats::abs(c_expected_host(k, i, j) - c_actual_host(k, i, j));
        }
      }
    }
  }

  EXPECT_NEAR_KK(diff / sum, 0, eps, fmsg + fmsg_rhs);
}

template <typename DeviceType, typename ViewType, typename ScalarType, typename ParamTagType>
void impl_test_batched_gemm(const int N, const int matAdim1, const int matAdim2, const int matBdim1, const int matBdim2,
                            const int matCdim1, const int matCdim2) {
  {
    BatchedGemmHandle batchedGemmHandle;

    ASSERT_EQ(batchedGemmHandle.get_kernel_algo_type(), BaseHeuristicAlgos::SQUARE);
    ASSERT_EQ(batchedGemmHandle.teamSz, 0);
    ASSERT_EQ(batchedGemmHandle.vecLen, 0);

#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
    cublasHandle_t cublas_handle;
    BatchedGemmHandle batchedGemmHandleCublas(cublas_handle, GemmTplAlgos::CUBLAS, 0, 0);
    ASSERT_EQ(&cublas_handle, batchedGemmHandleCublas.get_tpl_params());
    ASSERT_EQ(batchedGemmHandleCublas.get_kernel_algo_type(), (int)GemmTplAlgos::CUBLAS);
    ASSERT_EQ(batchedGemmHandleCublas.teamSz, 0);
    ASSERT_EQ(batchedGemmHandleCublas.vecLen, 0);
#endif

    // FIXME temporary workaround to run this magma test only if cublas is not
    // enabled the design of the BatchedGemmHandle currently does not allow
    // simultanous testing in this way. See issue #2177
#if defined(KOKKOSKERNELS_ENABLE_TPL_MAGMA) && !defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
    magma_queue_t magma_queue;
    BatchedGemmHandle batchedGemmHandleMagma(magma_queue, GemmTplAlgos::MAGMA, 0, 0);
    ASSERT_EQ(&magma_queue, batchedGemmHandleMagma.get_tpl_params());
    ASSERT_EQ(batchedGemmHandleMagma.get_kernel_algo_type(), (int)GemmTplAlgos::MAGMA);
    ASSERT_EQ(batchedGemmHandleMagma.teamSz, 0);
    ASSERT_EQ(batchedGemmHandleMagma.vecLen, 0);
#endif
  }

  for (int algo_type = BaseHeuristicAlgos::SQUARE; algo_type < GemmKokkosBatchedAlgos::N; ++algo_type) {
    {
      try {
        BatchedGemmHandle batchedGemmHandle(algo_type);

        ASSERT_EQ(batchedGemmHandle.get_kernel_algo_type(), algo_type);

        if (algo_type == BaseTplAlgos::ARMPL || algo_type == BaseKokkosBatchedAlgos::KK_SERIAL ||
            algo_type == GemmKokkosBatchedAlgos::KK_SERIAL_RANK0 || algo_type == GemmKokkosBatchedAlgos::KK_DBLBUF) {
          impl_test_batched_gemm_with_handle<DeviceType, ViewType, ScalarType, ParamTagType>(
              &batchedGemmHandle, N, matAdim1, matAdim2, matBdim1, matBdim2, matCdim1, matCdim2, 1.5, 3.0);
        } else if (algo_type == BaseHeuristicAlgos::SQUARE) {
          // Invoke 4 times to ensure we cover all paths for alpha and beta
          impl_test_batched_gemm_with_handle<DeviceType, ViewType, ScalarType, ParamTagType>(
              &batchedGemmHandle, N, matAdim1, matAdim2, matBdim1, matBdim2, matCdim1, matCdim2, 0.0, 0.0);
          impl_test_batched_gemm_with_handle<DeviceType, ViewType, ScalarType, ParamTagType>(
              &batchedGemmHandle, N, matAdim1, matAdim2, matBdim1, matBdim2, matCdim1, matCdim2, 1.0, 0.0);
          impl_test_batched_gemm_with_handle<DeviceType, ViewType, ScalarType, ParamTagType>(
              &batchedGemmHandle, N, matAdim1, matAdim2, matBdim1, matBdim2, matCdim1, matCdim2, 0.0, 1.0);
          impl_test_batched_gemm_with_handle<DeviceType, ViewType, ScalarType, ParamTagType>(
              &batchedGemmHandle, N, matAdim1, matAdim2, matBdim1, matBdim2, matCdim1, matCdim2, 1.5, 3.0);
        } else {
          try {
            // Allocate these views to invoke BatchedGemm with an unsupported
            // algo type
            ViewType a_actual("a_actual", N, matAdim1, matAdim2);
            ViewType b_actual("b_actual", N, matBdim1, matBdim2);
            ViewType c_actual("c_actual", N, matCdim1, matCdim2);
            using ta         = typename ParamTagType::transA;
            using tb         = typename ParamTagType::transB;
            using bl         = typename ParamTagType::batchLayout;
            ScalarType alpha = 0.34;
            ScalarType beta  = 0.43;
            BatchedGemm<ta, tb, bl>(&batchedGemmHandle, alpha, a_actual, b_actual, beta, c_actual);
            std::string fmsg = kk_failure_str(__FILE__, __FUNCTION__, __LINE__);
            FAIL() << fmsg;
          } catch (const std::runtime_error& error) {
            ;
          }
        }
      } catch (const std::runtime_error& error) {
#if !defined(KOKKOSKERNELS_ENABLE_TPL_ARMPL) || (ARMPL_BUILD < 1058)
        if (algo_type == BaseTplAlgos::ARMPL) {
          ;
        } else {
          std::string fmsg = kk_failure_str(__FILE__, __FUNCTION__, __LINE__);
          FAIL() << fmsg;
        }
#endif  // KOKKOSKERNELS_ENABLE_TPL_ARMPL
      }
    }
  }
}
}  // namespace Test

template <typename ViewType, typename DeviceType, typename ValueType, typename ScalarType, typename ParamTagType>
void test_batched_gemm_with_layout(int N) {
  // Square cases
  {
    int i = 0;
    Test::impl_test_batched_gemm<DeviceType, ViewType, ScalarType, ParamTagType>(N, i, i, i, i, i, i);

    i = 10;
    Test::impl_test_batched_gemm<DeviceType, ViewType, ScalarType, ParamTagType>(N, i, i, i, i, i, i);

    i = 25;
    Test::impl_test_batched_gemm<DeviceType, ViewType, ScalarType, ParamTagType>(N, i, i, i, i, i, i);

    i = 32;
    Test::impl_test_batched_gemm<DeviceType, ViewType, ScalarType, ParamTagType>(N, i, i, i, i, i, i);
  }

  // Non-square cases
  for (int i = 1; i < 5; ++i) {
    int dimM = 1 * i;
    int dimN = 2 * i;
    int dimK = 3 * i;
    if ((std::is_same<typename ParamTagType::transA, KokkosBatched::Trans::NoTranspose>::value) &&
        (std::is_same<typename ParamTagType::transB, KokkosBatched::Trans::NoTranspose>::value)) {
      Test::impl_test_batched_gemm<DeviceType, ViewType, ScalarType, ParamTagType>(N, dimM, dimK, dimK, dimN, dimM,
                                                                                   dimN);
    }
    if ((std::is_same<typename ParamTagType::transA, KokkosBatched::Trans::NoTranspose>::value) &&
        (std::is_same<typename ParamTagType::transB, KokkosBatched::Trans::Transpose>::value)) {
      Test::impl_test_batched_gemm<DeviceType, ViewType, ScalarType, ParamTagType>(N, dimM, dimK, dimN, dimK, dimM,
                                                                                   dimN);
    }
    if ((std::is_same<typename ParamTagType::transA, KokkosBatched::Trans::Transpose>::value) &&
        (std::is_same<typename ParamTagType::transB, KokkosBatched::Trans::NoTranspose>::value)) {
      Test::impl_test_batched_gemm<DeviceType, ViewType, ScalarType, ParamTagType>(N, dimK, dimM, dimK, dimN, dimM,
                                                                                   dimN);
    }
    if ((std::is_same<typename ParamTagType::transA, KokkosBatched::Trans::Transpose>::value) &&
        (std::is_same<typename ParamTagType::transB, KokkosBatched::Trans::Transpose>::value)) {
      Test::impl_test_batched_gemm<DeviceType, ViewType, ScalarType, ParamTagType>(N, dimK, dimM, dimN, dimK, dimM,
                                                                                   dimN);
    }
  }
}

template <typename DeviceType, typename ValueType, typename ScalarType, typename ParamTagType>
int test_batched_gemm() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  if constexpr (std::is_same_v<typename ParamTagType::batchLayout, typename BatchLayout::Right>) {
    using param_tag_type =
        ::Test::SharedParamTag<typename ParamTagType::transA, typename ParamTagType::transB, BatchLayout::Right>;
    typedef Kokkos::View<ValueType***, Kokkos::LayoutLeft, DeviceType> llVt;
    test_batched_gemm_with_layout<llVt, DeviceType, ValueType, ScalarType, param_tag_type>(0);
    test_batched_gemm_with_layout<llVt, DeviceType, ValueType, ScalarType, param_tag_type>(1);
    test_batched_gemm_with_layout<llVt, DeviceType, ValueType, ScalarType, param_tag_type>(4);
    test_batched_gemm_with_layout<llVt, DeviceType, ValueType, ScalarType, param_tag_type>(8);
    test_batched_gemm_with_layout<llVt, DeviceType, ValueType, ScalarType, param_tag_type>(16);
  } else {
    std::cerr << "TEST SKIPPED since BatchLayout is not Right." << std::endl;
  }
#else
  std::cerr << "TEST SKIPPED since LayoutLeft is not ETI'd." << std::endl;
#endif  // KOKKOSKERNELS_INST_LAYOUTLEFT

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  if constexpr (std::is_same_v<typename ParamTagType::batchLayout, typename BatchLayout::Left>) {
    using param_tag_type =
        ::Test::SharedParamTag<typename ParamTagType::transA, typename ParamTagType::transB, BatchLayout::Left>;
    typedef Kokkos::View<ValueType***, Kokkos::LayoutRight, DeviceType> lrVt;
    test_batched_gemm_with_layout<lrVt, DeviceType, ValueType, ScalarType, param_tag_type>(0);
    test_batched_gemm_with_layout<lrVt, DeviceType, ValueType, ScalarType, param_tag_type>(1);
    test_batched_gemm_with_layout<lrVt, DeviceType, ValueType, ScalarType, param_tag_type>(4);
    test_batched_gemm_with_layout<lrVt, DeviceType, ValueType, ScalarType, param_tag_type>(8);
    test_batched_gemm_with_layout<lrVt, DeviceType, ValueType, ScalarType, param_tag_type>(16);
  } else {
    std::cerr << "TEST SKIPPED since BatchLayout is not Left." << std::endl;
  }
#else
  std::cerr << "TEST SKIPPED since LayoutRight is not ETI'd." << std::endl;
#endif  // KOKKOSKERNELS_INST_LAYOUTRIGHT
  return 0;
}
