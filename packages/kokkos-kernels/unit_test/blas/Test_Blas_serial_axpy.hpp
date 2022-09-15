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
// Questions? Contact Luc Berger-Vergiat (lberge@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef TEST_BLAS_SERIAL_AXPY_HPP_
#define TEST_BLAS_SERIAL_AXPY_HPP_

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"
#include "KokkosKernels_TestUtils.hpp"

#include "KokkosBlas1_axpby.hpp"

namespace Test {

struct KokkosKernelAxpyTag {};
struct NaiveAxpyTag {};

template <typename DeviceType, typename ViewType, typename ScalarType,
          typename AlgoTagType>
struct Functor_TestBlasSerialAxpy {
  ScalarType _alpha;
  ViewType _x;
  ViewType _y;

  KOKKOS_INLINE_FUNCTION
  Functor_TestBlasSerialAxpy(const ScalarType alpha, const ViewType &x,
                             const ViewType &y)
      : _alpha(alpha), _x(x), _y(y) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const KokkosKernelAxpyTag &, const int i) const {
    auto X = Kokkos::subview(_x, i, Kokkos::ALL(), Kokkos::ALL());
    auto Y = Kokkos::subview(_y, i, Kokkos::ALL(), Kokkos::ALL());
    KokkosBlas::serial_axpy(_alpha, X, Y);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const NaiveAxpyTag &, const int k) const {
    auto X      = Kokkos::subview(_x, k, Kokkos::ALL(), Kokkos::ALL());
    auto Y      = Kokkos::subview(_y, k, Kokkos::ALL(), Kokkos::ALL());
    const int m = X.extent(0), n = X.extent(1);
    for (int i = 0; i < m; ++i)
      for (int j = 0; j < n; ++j) Y(i, j) += _alpha * X(i, j);
  }

  inline void run() {
    using value_type = typename ViewType::value_type;
    std::string name_region("KokkosBlas::Test::SerialAxpy");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name_work_tag =
        (std::is_same<AlgoTagType, KokkosKernelAxpyTag>::value
             ? "::KokkosBlas"
             : std::is_same<AlgoTagType, NaiveAxpyTag>::value
                   ? "::Naive"
                   : "::UnknownWorkTag");
    std::string name_test_id = "Axpy";
    std::string name =
        name_region + name_value_type + name_work_tag + name_test_id;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<DeviceType, AlgoTagType> policy(0, _x.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
    return;
  }
};

template <typename DeviceType, typename ViewType, typename ScalarType>
void impl_test_blas_serial_axpy(const int N, const int BlkSize) {
  /// typedefs
  using value_type = typename ViewType::value_type;
  using ats        = Kokkos::ArithTraits<value_type>;

  /// radomized input testing views
  const ScalarType alpha = 11.1;
  ViewType X("X", N, BlkSize, BlkSize);
  ViewType Y("Y", N, BlkSize, BlkSize);
  ViewType Yref("Yref", N, BlkSize, BlkSize);

  Kokkos::Random_XorShift64_Pool<typename DeviceType::execution_space> random(
      13718);
  Kokkos::fill_random(X, random, ats::one());
  Kokkos::fill_random(Y, random, ats::one());
  Kokkos::fence();
  Kokkos::deep_copy(Yref, Y);

  /// test body
  Functor_TestBlasSerialAxpy<DeviceType, ViewType, ScalarType, NaiveAxpyTag>(
      alpha, X, Yref)
      .run();
  Functor_TestBlasSerialAxpy<DeviceType, ViewType, ScalarType,
                             KokkosKernelAxpyTag>(alpha, X, Y)
      .run();

  Kokkos::fence();

  /// for comparison send it to host
  typename ViewType::HostMirror Y_host    = Kokkos::create_mirror_view(Y);
  typename ViewType::HostMirror Yref_host = Kokkos::create_mirror_view(Yref);

  Kokkos::deep_copy(Y_host, Y);
  Kokkos::deep_copy(Yref_host, Yref);

  /// check a = b
  typename ats::mag_type eps =
      100 * std::numeric_limits<typename ats::mag_type>::epsilon();
  for (int k = 0; k < N; ++k)
    for (int i = 0; i < BlkSize; ++i)
      for (int j = 0; j < BlkSize; ++j)
        EXPECT_NEAR_KK(Y_host(k, i, j), Yref_host(k, i, j), eps);
}

}  // namespace Test

template <typename DeviceType, typename ValueType, typename ScalarType>
int test_blas_serial_axpy() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    typedef Kokkos::View<ValueType ***, Kokkos::LayoutLeft, DeviceType>
        ViewType;
    Test::impl_test_blas_serial_axpy<DeviceType, ViewType, ScalarType>(0, 10);
    Test::impl_test_blas_serial_axpy<DeviceType, ViewType, ScalarType>(10, 15);
    Test::impl_test_blas_serial_axpy<DeviceType, ViewType, ScalarType>(1024, 9);
    Test::impl_test_blas_serial_axpy<DeviceType, ViewType, ScalarType>(132231,
                                                                       3);
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    typedef Kokkos::View<ValueType ***, Kokkos::LayoutRight, DeviceType>
        ViewType;
    Test::impl_test_blas_serial_axpy<DeviceType, ViewType, ScalarType>(0, 10);
    Test::impl_test_blas_serial_axpy<DeviceType, ViewType, ScalarType>(10, 15);
    Test::impl_test_blas_serial_axpy<DeviceType, ViewType, ScalarType>(1024, 9);
    Test::impl_test_blas_serial_axpy<DeviceType, ViewType, ScalarType>(132231,
                                                                       3);
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, serial_axpy_float_float) {
  test_blas_serial_axpy<TestExecSpace, float, float>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, serial_axpy_double_double) {
  test_blas_serial_axpy<TestExecSpace, double, double>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, serial_axpy_dcomplex_dcomplex) {
  test_blas_serial_axpy<TestExecSpace, Kokkos::complex<double>,
                        Kokkos::complex<double> >();
}

TEST_F(TestCategory, serial_axpy_dcomplex_double) {
  test_blas_serial_axpy<TestExecSpace, Kokkos::complex<double>, double>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, serial_axpy_fcomplex_fcomplex) {
  test_blas_serial_axpy<TestExecSpace, Kokkos::complex<float>,
                        Kokkos::complex<double> >();
}

TEST_F(TestCategory, serial_axpy_fcomplex_float) {
  test_blas_serial_axpy<TestExecSpace, Kokkos::complex<float>, float>();
}
#endif

#endif  // TEST_BLAS_SERIAL_AXPY_HPP_
