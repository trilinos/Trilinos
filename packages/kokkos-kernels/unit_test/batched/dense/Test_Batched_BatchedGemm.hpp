#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

#include "KokkosBatched_Gemm_Decl.hpp"

#include "KokkosKernels_TestUtils.hpp"

using namespace KokkosBatched;

namespace Test {
template <typename DeviceType, typename ViewType, typename ScalarType,
          typename ParamTagType>
void impl_test_batched_gemm_with_handle(BatchedGemmHandle* batchedGemmHandle,
                                        const int N, const int matAdim1,
                                        const int matAdim2, const int matBdim1,
                                        const int matBdim2, const int matCdim1,
                                        const int matCdim2, ScalarType alpha,
                                        ScalarType beta) {
  using execution_space = typename DeviceType::execution_space;
  using transA          = typename ParamTagType::transA;
  using transB          = typename ParamTagType::transB;
  using batchLayout     = typename ParamTagType::batchLayout;
  using view_layout     = typename ViewType::array_layout;
  using ats             = Kokkos::Details::ArithTraits<ScalarType>;

  int ret        = 0;
  auto algo_type = batchedGemmHandle->get_kernel_algo_type();

  ViewType a_expected, a_actual, b_expected, b_actual, c_expected, c_actual;
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
      Impl::BatchedDblBufGemm<transA, transB, batchLayout, BatchedGemmHandle,
                              ScalarType, decltype(a_actual),
                              decltype(b_actual), decltype(c_actual),
                              BoundsCheck::Yes, 65536, 1, 65536>(
          batchedGemmHandle, alpha, a_actual, b_actual, beta, c_actual)
          .invoke();
      FAIL();
    } catch (const std::runtime_error& error) {
      ;
    }

    // Check for DblBuf runtime errors related to vector_len
    try {
      Impl::BatchedDblBufGemm<transA, transB, batchLayout, BatchedGemmHandle,
                              ScalarType, decltype(a_actual),
                              decltype(b_actual), decltype(c_actual),
                              BoundsCheck::No, 65536, 65536 * 2, 65536>(
          batchedGemmHandle, alpha, a_actual, b_actual, beta, c_actual)
          .invoke();
      FAIL();
    } catch (const std::runtime_error& error) {
      ;
    }
  }

  // Check for expected BatchedGemm runtime errors
  try {
    ret = BatchedGemm<transA, transB, batchLayout>(
        batchedGemmHandle, alpha, a_actual, b_actual, beta,
        c_actual);  // Compute c_actual
  } catch (const std::runtime_error& error) {
    // std::cout << "Caught expected runtime error" << std::endl;
    if (algo_type == BaseHeuristicAlgos::SQUARE && matCdim1 != matCdim2)
      ;
    else if (!((std::is_same<view_layout, Kokkos::LayoutLeft>::value &&
                !std::is_same<batchLayout, BatchLayout::Right>::value) ||
               (std::is_same<view_layout, Kokkos::LayoutRight>::value &&
                !std::is_same<batchLayout, BatchLayout::Left>::value))) {
      FAIL();
    }
    return;
  }
  ASSERT_EQ(ret, 0);

  Functor_BatchedVanillaGEMM<ViewType, ViewType, ViewType, execution_space>
      vgemm;
  vgemm.A_t = std::is_same<transA, Trans::Transpose>::value;
  vgemm.B_t = std::is_same<transB, Trans::Transpose>::value;
  vgemm.batch_size_last_dim =
      std::is_same<batchLayout, BatchLayout::Right>::value;
  vgemm.A_c = vgemm.B_c = false;
  vgemm.A               = a_expected;
  vgemm.B               = b_expected;
  vgemm.C               = c_expected;
  vgemm.alpha           = alpha;
  vgemm.beta            = beta;
  vgemm.run();  // Compute c_expected

  Kokkos::fence();

  typename ViewType::HostMirror c_expected_host =
      Kokkos::create_mirror_view(c_expected);
  typename ViewType::HostMirror c1_host = Kokkos::create_mirror_view(c_actual);

  // Copy to host
  Kokkos::deep_copy(c_expected_host, c_expected);
  Kokkos::deep_copy(c1_host, c_actual);

  Kokkos::fence();

  // check c_expected = c_actual ; this eps is about 2^-9
  // Set mag_type to host_value_type, we may not have half precision on host
  using mag_type = float;
  mag_type sum(1), diff(0);

  mag_type eps = (mag_type)(1 << 1) * KOKKOSKERNELS_IMPL_FP16_EPSILON;

  for (int k = 0; k < N; ++k) {
    for (int i = 0; i < matCdim1; ++i) {
      for (int j = 0; j < matCdim2; ++j) {
        if (std::is_same<batchLayout, BatchLayout::Right>::value) {
          sum += ats::abs(c_expected_host(i, j, k));
          diff += ats::abs(c_expected_host(i, j, k) - c1_host(i, j, k));
        } else {
          sum += ats::abs(c_expected_host(k, i, j));
          diff += ats::abs(c_expected_host(k, i, j) - c1_host(k, i, j));
        }
      }
    }
  }
  EXPECT_NEAR_KK(diff / sum, 0, eps);
}

template <typename DeviceType, typename ViewType, typename ScalarType,
          typename ParamTagType>
void impl_test_batched_gemm(const int N, const int matAdim1, const int matAdim2,
                            const int matBdim1, const int matBdim2,
                            const int matCdim1, const int matCdim2) {
  {
    BatchedGemmHandle batchedGemmHandle;

    ASSERT_EQ(batchedGemmHandle.get_kernel_algo_type(),
              BaseHeuristicAlgos::SQUARE);
    ASSERT_EQ(batchedGemmHandle.teamSz, 0);
    ASSERT_EQ(batchedGemmHandle.vecLen, 0);

#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
    cublasHandle_t cublas_handle;
    BatchedGemmHandle batchedGemmHandleCublas(cublas_handle,
                                              GemmTplAlgos::CUBLAS, 0, 0);
    ASSERT_EQ(&cublas_handle, batchedGemmHandleCublas.get_tpl_params());
    ASSERT_EQ(batchedGemmHandleCublas.get_kernel_algo_type(),
              (int)GemmTplAlgos::CUBLAS);
    ASSERT_EQ(batchedGemmHandleCublas.teamSz, 0);
    ASSERT_EQ(batchedGemmHandleCublas.vecLen, 0);
#endif

#if defined(KOKKOSKERNELS_ENABLE_TPL_MAGMA)
    magma_queue_t magma_queue;
    BatchedGemmHandle batchedGemmHandleMagma(magma_queue, GemmTplAlgos::MAGMA,
                                             0, 0);
    ASSERT_EQ(&magma_queue, batchedGemmHandleMagma.get_tpl_params());
    ASSERT_EQ(batchedGemmHandleMagma.get_kernel_algo_type(),
              (int)GemmTplAlgos::MAGMA);
    ASSERT_EQ(batchedGemmHandleMagma.teamSz, 0);
    ASSERT_EQ(batchedGemmHandleMagma.vecLen, 0);
#endif
  }

  for (int algo_type = BaseHeuristicAlgos::SQUARE;
       algo_type < GemmKokkosBatchedAlgos::N; ++algo_type) {
    {
      BatchedGemmHandle batchedGemmHandle(algo_type);

      ASSERT_EQ(batchedGemmHandle.get_kernel_algo_type(), algo_type);

      if (algo_type == BaseKokkosBatchedAlgos::KK_SERIAL ||
          algo_type == BaseHeuristicAlgos::SQUARE ||
          algo_type == GemmKokkosBatchedAlgos::KK_DBLBUF ||
          algo_type == GemmKokkosBatchedAlgos::KK_SERIAL_RANK0) {
        // Invoke 4 times to ensure we cover all paths for alpha and beta
        impl_test_batched_gemm_with_handle<DeviceType, ViewType, ScalarType,
                                           ParamTagType>(
            &batchedGemmHandle, N, matAdim1, matAdim2, matBdim1, matBdim2,
            matCdim1, matCdim2, 0.0, 0.0);
        impl_test_batched_gemm_with_handle<DeviceType, ViewType, ScalarType,
                                           ParamTagType>(
            &batchedGemmHandle, N, matAdim1, matAdim2, matBdim1, matBdim2,
            matCdim1, matCdim2, 1.0, 0.0);
        impl_test_batched_gemm_with_handle<DeviceType, ViewType, ScalarType,
                                           ParamTagType>(
            &batchedGemmHandle, N, matAdim1, matAdim2, matBdim1, matBdim2,
            matCdim1, matCdim2, 0.0, 1.0);
        impl_test_batched_gemm_with_handle<DeviceType, ViewType, ScalarType,
                                           ParamTagType>(
            &batchedGemmHandle, N, matAdim1, matAdim2, matBdim1, matBdim2,
            matCdim1, matCdim2, 1.5, 3.0);
      } else {
        try {
          // Allocate these views to invoke BatchedGemm with an unsupported algo
          // type
          ViewType a_actual("a_actual", N, matAdim1, matAdim2);
          ViewType b_actual("b_actual", N, matBdim1, matBdim2);
          ViewType c_actual("c_actual", N, matCdim1, matCdim2);
          using ta = typename ParamTagType::transA;
          using tb = typename ParamTagType::transB;
          using bl = typename ParamTagType::batchLayout;
          BatchedGemm<ta, tb, bl>(&batchedGemmHandle, 0.34, a_actual, b_actual,
                                  0.43, c_actual);
          FAIL();
        } catch (const std::runtime_error& error) {
          ;
        }
      }
    }
  }
}
}  // namespace Test

template <typename ViewType, typename DeviceType, typename ValueType,
          typename ScalarType, typename ParamTagType>
void test_batched_gemm_with_layout() {
  // Square cases
  for (int i = 0; i < 5; ++i) {
    Test::impl_test_batched_gemm<DeviceType, ViewType, ScalarType,
                                 ParamTagType>(4, i, i, i, i, i, i);
  }

  {
    int i = 10;
    Test::impl_test_batched_gemm<DeviceType, ViewType, ScalarType,
                                 ParamTagType>(0, i, i, i, i, i, i);

    i = 25;
    Test::impl_test_batched_gemm<DeviceType, ViewType, ScalarType,
                                 ParamTagType>(8, i, i, i, i, i, i);

    i = 32;
    Test::impl_test_batched_gemm<DeviceType, ViewType, ScalarType,
                                 ParamTagType>(8, i, i, i, i, i, i);
  }

  // Non-square cases
  for (int i = 0; i < 5; ++i) {
    int dimM = 1 * i;
    int dimN = 2 * i;
    int dimK = 3 * i;
    if ((std::is_same<typename ParamTagType::transA,
                      KokkosBatched::Trans::NoTranspose>::value) &&
        (std::is_same<typename ParamTagType::transB,
                      KokkosBatched::Trans::NoTranspose>::value)) {
      Test::impl_test_batched_gemm<DeviceType, ViewType, ScalarType,
                                   ParamTagType>(16, dimM, dimK, dimK, dimN,
                                                 dimM, dimN);
    }
    if ((std::is_same<typename ParamTagType::transA,
                      KokkosBatched::Trans::NoTranspose>::value) &&
        (std::is_same<typename ParamTagType::transB,
                      KokkosBatched::Trans::Transpose>::value)) {
      Test::impl_test_batched_gemm<DeviceType, ViewType, ScalarType,
                                   ParamTagType>(16, dimM, dimK, dimN, dimK,
                                                 dimM, dimN);
    }
    if ((std::is_same<typename ParamTagType::transA,
                      KokkosBatched::Trans::Transpose>::value) &&
        (std::is_same<typename ParamTagType::transB,
                      KokkosBatched::Trans::NoTranspose>::value)) {
      Test::impl_test_batched_gemm<DeviceType, ViewType, ScalarType,
                                   ParamTagType>(16, dimK, dimM, dimK, dimN,
                                                 dimM, dimN);
    }
    if ((std::is_same<typename ParamTagType::transA,
                      KokkosBatched::Trans::Transpose>::value) &&
        (std::is_same<typename ParamTagType::transB,
                      KokkosBatched::Trans::Transpose>::value)) {
      Test::impl_test_batched_gemm<DeviceType, ViewType, ScalarType,
                                   ParamTagType>(16, dimK, dimM, dimN, dimK,
                                                 dimM, dimN);
    }
  }
}

template <typename DeviceType, typename ValueType, typename ScalarType,
          typename ParamTagType>
int test_batched_gemm() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  typedef Kokkos::View<ValueType***, Kokkos::LayoutLeft, DeviceType> llVt;
  test_batched_gemm_with_layout<llVt, DeviceType, ValueType, ScalarType,
                                ParamTagType>();
#endif  // KOKKOSKERNELS_INST_LAYOUTLEFT

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  typedef Kokkos::View<ValueType***, Kokkos::LayoutRight, DeviceType> lrVt;
  test_batched_gemm_with_layout<lrVt, DeviceType, ValueType, ScalarType,
                                ParamTagType>();
#endif  // KOKKOSKERNELS_INST_LAYOUTRIGHT
  return 0;
}
