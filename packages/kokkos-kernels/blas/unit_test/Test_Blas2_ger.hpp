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

// **********************************************************************
// The tests executed by the code below cover many combinations for
// the operation A += alpha * x * y^{T,H}.
// 01) Type of 'x' components: float, double, complex, ...
// 02) Type of 'y' components: float, double, complex, ...
// 03) Type of 'A' components: float, double, complex, ...
// 04) Execution space: serial, threads, OpenMP, Cuda, ...
// 05) Layout of 'x'
// 06) Layout of 'y'
// 07) Layout of 'A'
// 08) Dimension of 'A'
// 09) Options 'const' or 'non const' for x view, when calling ger()
// 10) Options 'const' or 'non const' for y view, when calling ger()
// 11) Usage of analytical results in the tests
// 12) Options 'T' or 'H' when calling ger()
//
// Choices (01)-(04) are selected in the routines TEST_F() at the
// very bottom of the file, when calling test_ger<...>().
//
// Choices (05)-(12) are selected in routine test_gerr<...>(),
// when calling the method test() of class Test::GerTester<...>.
//
// The class Test::GerTester<...> represents the "core" of the test
// logic, where all calculations, comparisons, and success/failure
// decisions are performed.
//
// A high level explanation of method Test::GerTester<...>::test()
// is given by the 9 steps named "Step 1 of 9" to "Step 9 of 9"
// in the code below.
// **********************************************************************

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBlas2_ger.hpp>
#include <Kokkos_MathematicalConstants.hpp>

namespace Test {

template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
class GerTester {
 public:
  GerTester();

  ~GerTester();

  void test(const int M, const int N, const int nonConstConstCombinations, const bool useAnalyticalResults = false,
            const bool useHermitianOption = false);

 private:
  using _ViewTypeX = Kokkos::View<ScalarX*, tLayoutX, Device>;
  using _ViewTypeY = Kokkos::View<ScalarY*, tLayoutY, Device>;
  using _ViewTypeA = Kokkos::View<ScalarA**, tLayoutA, Device>;

  using _HostViewTypeX    = typename _ViewTypeX::HostMirror;
  using _HostViewTypeY    = typename _ViewTypeY::HostMirror;
  using _HostViewTypeA    = typename _ViewTypeA::HostMirror;
  using _ViewTypeExpected = Kokkos::View<ScalarA**, tLayoutA, Kokkos::HostSpace>;

  using _KAT_A   = Kokkos::ArithTraits<ScalarA>;
  using _AuxType = typename _KAT_A::mag_type;

  void populateVariables(ScalarA& alpha, view_stride_adapter<_ViewTypeX, false>& x,
                         view_stride_adapter<_ViewTypeY, false>& y, view_stride_adapter<_ViewTypeA, false>& A,
                         _ViewTypeExpected& h_expected, bool& expectedResultIsKnown);

  template <class T>
  typename std::enable_if<
      std::is_same<T, Kokkos::complex<float>>::value || std::is_same<T, Kokkos::complex<double>>::value, void>::type
  populateAnalyticalValues(T& alpha, _HostViewTypeX& h_x, _HostViewTypeY& h_y, _HostViewTypeA& h_A,
                           _ViewTypeExpected& h_expected);

  template <class T>
  typename std::enable_if<
      !std::is_same<T, Kokkos::complex<float>>::value && !std::is_same<T, Kokkos::complex<double>>::value, void>::type
  populateAnalyticalValues(T& alpha, _HostViewTypeX& h_x, _HostViewTypeY& h_y, _HostViewTypeA& h_A,
                           _ViewTypeExpected& h_expected);

  template <class T>
  typename std::enable_if<
      std::is_same<T, Kokkos::complex<float>>::value || std::is_same<T, Kokkos::complex<double>>::value, void>::type
  populateVanillaValues(const T& alpha, const _HostViewTypeX& h_x, const _HostViewTypeY& h_y, const _HostViewTypeA& h_A,
                        _ViewTypeExpected& h_vanilla);

  template <class T>
  typename std::enable_if<
      !std::is_same<T, Kokkos::complex<float>>::value && !std::is_same<T, Kokkos::complex<double>>::value, void>::type
  populateVanillaValues(const T& alpha, const _HostViewTypeX& h_x, const _HostViewTypeY& h_y, const _HostViewTypeA& h_A,
                        _ViewTypeExpected& h_vanilla);

  template <class T>
  typename std::enable_if<
      std::is_same<T, Kokkos::complex<float>>::value || std::is_same<T, Kokkos::complex<double>>::value, void>::type
  compareVanillaAgainstExpected(const T& alpha, const _ViewTypeExpected& h_vanilla,
                                const _ViewTypeExpected& h_expected);

  template <class T>
  typename std::enable_if<
      !std::is_same<T, Kokkos::complex<float>>::value && !std::is_same<T, Kokkos::complex<double>>::value, void>::type
  compareVanillaAgainstExpected(const T& alpha, const _ViewTypeExpected& h_vanilla,
                                const _ViewTypeExpected& h_expected);

  template <class T>
  typename std::enable_if<
      std::is_same<T, Kokkos::complex<float>>::value || std::is_same<T, Kokkos::complex<double>>::value, void>::type
  compareKkGerAgainstExpected(const T& alpha, const _HostViewTypeA& h_A, const _ViewTypeExpected& h_expected);

  template <class T>
  typename std::enable_if<
      !std::is_same<T, Kokkos::complex<float>>::value && !std::is_same<T, Kokkos::complex<double>>::value, void>::type
  compareKkGerAgainstExpected(const T& alpha, const _HostViewTypeA& h_A, const _ViewTypeExpected& h_expected);

  template <class T>
  T shrinkAngleToZeroTwoPiRange(const T input);

  template <class TX, class TY>
  void callKkGerAndCompareAgainstExpected(const ScalarA& alpha, TX& x, TY& y, view_stride_adapter<_ViewTypeA, false>& A,
                                          const _ViewTypeExpected& h_expected, const std::string& situation);

  const bool _A_is_complex;
  const bool _A_is_lr;
  const bool _A_is_ll;
  const bool _testIsGpu;
  const bool _vanillaUsesDifferentOrderOfOps;
  const _AuxType _absTol;
  const _AuxType _relTol;
  int _M;
  int _N;
  bool _useAnalyticalResults;
  bool _useHermitianOption;
  bool _kkGerShouldThrowException;
};

template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
GerTester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::GerTester()
    : _A_is_complex(std::is_same<ScalarA, Kokkos::complex<float>>::value ||
                    std::is_same<ScalarA, Kokkos::complex<double>>::value),
      _A_is_lr(std::is_same<tLayoutA, Kokkos::LayoutRight>::value),
      _A_is_ll(std::is_same<tLayoutA, Kokkos::LayoutLeft>::value),
      _testIsGpu(KokkosKernels::Impl::kk_is_gpu_exec_space<typename Device::execution_space>())
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
      ,
      _vanillaUsesDifferentOrderOfOps(_A_is_lr && _testIsGpu)
#else
      ,
      _vanillaUsesDifferentOrderOfOps(false)
#endif
      ,
      // ****************************************************************
      // Tolerances for double can be tighter than tolerances for float.
      //
      // In the case of calculations with float, a small amount of
      // discrepancies between reference results and CUDA results are
      // large enough to require 'relTol' to value 5.0e-3. The same
      // calculations show no discrepancies for calculations with double.
      // ****************************************************************
      _absTol(std::is_same<_AuxType, float>::value ? 1.0e-6 : (std::is_same<_AuxType, double>::value ? 1.0e-9 : 0)),
      _relTol(std::is_same<_AuxType, float>::value ? 5.0e-3 : (std::is_same<_AuxType, double>::value ? 1.0e-6 : 0)),
      _M(-1),
      _N(-1),
      _useAnalyticalResults(false),
      _useHermitianOption(false),
      _kkGerShouldThrowException(false) {
}

template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
GerTester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::~GerTester() {
  // Nothing to do
}

template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
void GerTester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::test(
    const int M, const int N, const int nonConstConstCombinations, const bool useAnalyticalResults,
    const bool useHermitianOption) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Entering GerTester::test()... - - - - - - - - - - - - - - - - "
               "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
               "- - - - - - - - - "
            << std::endl;

  std::cout << "_A_is_complex = " << _A_is_complex << ", _A_is_lr = " << _A_is_lr << ", _A_is_ll = " << _A_is_ll
            << ", _testIsGpu = " << _testIsGpu
            << ", _vanillaUsesDifferentOrderOfOps = " << _vanillaUsesDifferentOrderOfOps << ", _absTol = " << _absTol
            << ", _relTol = " << _relTol << std::endl;
#endif
  // ********************************************************************
  // Step 1 of 9: declare main types and variables
  // ********************************************************************
  _M                    = M;
  _N                    = N;
  _useAnalyticalResults = useAnalyticalResults;
  _useHermitianOption   = useHermitianOption;

#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
  _kkGerShouldThrowException = false;
  if (_A_is_complex && _useHermitianOption) {
    if ((_testIsGpu == false) && (_A_is_ll == false)) {
      _kkGerShouldThrowException = true;
    } else if ((_testIsGpu == true) && (_A_is_ll == false)) {
      _kkGerShouldThrowException = true;
    }
  }
#endif

  bool test_x_y(false);
  bool test_cx_y(false);
  bool test_x_cy(false);
  bool test_cx_cy(false);
  if (nonConstConstCombinations == 0) {
    test_x_y = true;
  } else if (nonConstConstCombinations == 1) {
    test_cx_y = true;
  } else if (nonConstConstCombinations == 2) {
    test_x_cy = true;
  } else if (nonConstConstCombinations == 3) {
    test_cx_cy = true;
  } else {
    test_x_y   = true;
    test_cx_y  = true;
    test_x_cy  = true;
    test_cx_cy = true;
  }

  view_stride_adapter<_ViewTypeX, false> x("X", _M);
  view_stride_adapter<_ViewTypeY, false> y("Y", _N);
  view_stride_adapter<_ViewTypeA, false> A("A", _M, _N);

  view_stride_adapter<_ViewTypeExpected, true> h_expected("expected A += alpha * x * y^{t,h}", _M, _N);
  bool expectedResultIsKnown = false;

  ScalarA alpha(0.);

  // ********************************************************************
  // Step 2 of 9: populate alpha, h_x, h_y, h_A, h_expected, x, y, A
  // ********************************************************************
  this->populateVariables(alpha, x, y, A, h_expected.d_view, expectedResultIsKnown);

  // ********************************************************************
  // Step 3 of 9: populate h_vanilla
  // ********************************************************************
  view_stride_adapter<_ViewTypeExpected, true> h_vanilla("vanilla = A + alpha * x * y^{t,h}", _M, _N);
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  Kokkos::printf("In Test_Blas2_ger.hpp, computing vanilla A with alpha type = %s\n", typeid(alpha).name());
#endif
  this->populateVanillaValues(alpha, x.h_view, y.h_view, A.h_view, h_vanilla.d_view);

  // ********************************************************************
  // Step 4 of 9: use h_vanilla and h_expected as appropriate
  // ********************************************************************
  if (expectedResultIsKnown) {
    // ******************************************************************
    // Compare h_vanilla against h_expected
    // ******************************************************************
    this->compareVanillaAgainstExpected(alpha, h_vanilla.d_view, h_expected.d_view);
  } else {
    // ******************************************************************
    // Copy h_vanilla to h_expected
    // ******************************************************************
    Kokkos::deep_copy(h_expected.d_base, h_vanilla.d_base);
  }

  // ********************************************************************
  // Step 5 of 9: test with 'non const x' and 'non const y'
  // ********************************************************************
  view_stride_adapter<_ViewTypeA, false> org_A("Org_A", _M, _N);
  Kokkos::deep_copy(org_A.d_base, A.d_base);

  if (test_x_y) {
    this->callKkGerAndCompareAgainstExpected(alpha, x.d_view, y.d_view, A, h_expected.d_view, "non const {x,y}");
  }

  // ********************************************************************
  // Step 6 of 9: test with const x
  // ********************************************************************
  if (test_cx_y) {
    Kokkos::deep_copy(A.d_base, org_A.d_base);

    this->callKkGerAndCompareAgainstExpected(alpha, x.d_view_const, y.d_view, A, h_expected.d_view, "const x");
  }

  // ********************************************************************
  // Step 7 of 9: test with const y
  // ********************************************************************
  if (test_x_cy) {
    Kokkos::deep_copy(A.d_base, org_A.d_base);

    this->callKkGerAndCompareAgainstExpected(alpha, x.d_view, y.d_view_const, A, h_expected.d_view, "const y");
  }

  // ********************************************************************
  // Step 8 of 9: test with const x and const y
  // ********************************************************************
  if (test_cx_cy) {
    Kokkos::deep_copy(A.d_base, org_A.d_base);

    this->callKkGerAndCompareAgainstExpected(alpha, x.d_view_const, y.d_view_const, A, h_expected.d_view,
                                             "const {x,y}");
  }

  // ********************************************************************
  // Step 9 of 9: tests with invalid values on the first input parameter
  // ********************************************************************
  EXPECT_ANY_THROW(KokkosBlas::ger(".", alpha, x.d_view, y.d_view, A.d_view))
      << "Failed test: kk ger should have thrown an exception for mode '.'";
  EXPECT_ANY_THROW(KokkosBlas::ger("", alpha, x.d_view, y.d_view, A.d_view))
      << "Failed test: kk ger should have thrown an exception for mode ''";

#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Leaving GerTester::test() - - - - - - - - - - - - - - - - - - "
               "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
               "- - - - - - - "
            << std::endl;
#endif
}

template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
void GerTester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::populateVariables(
    ScalarA& alpha, view_stride_adapter<_ViewTypeX, false>& x, view_stride_adapter<_ViewTypeY, false>& y,
    view_stride_adapter<_ViewTypeA, false>& A, _ViewTypeExpected& h_expected, bool& expectedResultIsKnown) {
  expectedResultIsKnown = false;

  if (_useAnalyticalResults) {
    this->populateAnalyticalValues(alpha, x.h_view, y.h_view, A.h_view, h_expected);
    Kokkos::deep_copy(x.d_base, x.h_base);
    Kokkos::deep_copy(y.d_base, y.h_base);
    Kokkos::deep_copy(A.d_base, A.h_base);

    expectedResultIsKnown = true;
  } else if ((_M == 1) && (_N == 1)) {
    alpha = 3;

    x.h_view[0] = 2;

    y.h_view[0] = 3;

    A.h_view(0, 0) = 7;

    Kokkos::deep_copy(x.d_base, x.h_base);
    Kokkos::deep_copy(y.d_base, y.h_base);
    Kokkos::deep_copy(A.d_base, A.h_base);

    h_expected(0, 0)      = 25;
    expectedResultIsKnown = true;
  } else if ((_M == 1) && (_N == 2)) {
    alpha = 3;

    x.h_view[0] = 2;

    y.h_view[0] = 3;
    y.h_view[1] = 4;

    A.h_view(0, 0) = 7;
    A.h_view(0, 1) = -6;

    Kokkos::deep_copy(x.d_base, x.h_base);
    Kokkos::deep_copy(y.d_base, y.h_base);
    Kokkos::deep_copy(A.d_base, A.h_base);

    h_expected(0, 0)      = 25;
    h_expected(0, 1)      = 18;
    expectedResultIsKnown = true;
  } else if ((_M == 2) && (_N == 2)) {
    alpha = 3;

    x.h_view[0] = 2;
    x.h_view[1] = 9;

    y.h_view[0] = -3;
    y.h_view[1] = 7;

    A.h_view(0, 0) = 17;
    A.h_view(0, 1) = -43;
    A.h_view(1, 0) = 29;
    A.h_view(1, 1) = 101;

    Kokkos::deep_copy(x.d_base, x.h_base);
    Kokkos::deep_copy(y.d_base, y.h_base);
    Kokkos::deep_copy(A.d_base, A.h_base);

    h_expected(0, 0)      = -1;
    h_expected(0, 1)      = -1;
    h_expected(1, 0)      = -52;
    h_expected(1, 1)      = 290;
    expectedResultIsKnown = true;
  } else {
    alpha = 3;

    Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

    {
      ScalarX randStart, randEnd;
      Test::getRandomBounds(1.0, randStart, randEnd);
      Kokkos::fill_random(x.d_view, rand_pool, randStart, randEnd);
    }

    {
      ScalarY randStart, randEnd;
      Test::getRandomBounds(1.0, randStart, randEnd);
      Kokkos::fill_random(y.d_view, rand_pool, randStart, randEnd);
    }

    {
      ScalarA randStart, randEnd;
      Test::getRandomBounds(1.0, randStart, randEnd);
      Kokkos::fill_random(A.d_view, rand_pool, randStart, randEnd);
    }

    Kokkos::deep_copy(x.h_base, x.d_base);
    Kokkos::deep_copy(y.h_base, y.d_base);
    Kokkos::deep_copy(A.h_base, A.d_base);
  }
}

// Code for complex values
template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
template <class T>
typename std::enable_if<
    std::is_same<T, Kokkos::complex<float>>::value || std::is_same<T, Kokkos::complex<double>>::value, void>::type
GerTester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::populateAnalyticalValues(
    T& alpha, _HostViewTypeX& h_x, _HostViewTypeY& h_y, _HostViewTypeA& h_A, _ViewTypeExpected& h_expected) {
  _AuxType auxI(0.);
  _AuxType auxJ(0.);
  _AuxType auxIpJ(0.);
  _AuxType auxImJ(0.);

  alpha.real() = 1.;
  alpha.imag() = -1.;

  for (int i = 0; i < _M; ++i) {
    auxI          = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i));
    h_x[i].real() = sin(auxI);
    h_x[i].imag() = cos(auxI);
  }

  for (int j = 0; j < _N; ++j) {
    auxJ          = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(j));
    h_y[j].real() = cos(auxJ);
    h_y[j].imag() = sin(auxJ);
  }

  if (_useHermitianOption) {
    for (int i = 0; i < _M; ++i) {
      auxI = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i));
      for (int j = 0; j < _N; ++j) {
        auxJ             = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(j));
        auxIpJ           = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i + j));
        h_A(i, j).real() = -sin(auxIpJ) - sin(auxI) * sin(auxJ) - cos(auxI) * cos(auxJ);
        h_A(i, j).imag() = -sin(auxIpJ) - sin(auxI) * sin(auxJ) + cos(auxI) * cos(auxJ);
      }
    }
  } else {
    for (int i = 0; i < _M; ++i) {
      auxI = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i));
      for (int j = 0; j < _N; ++j) {
        auxJ             = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(j));
        auxImJ           = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i - j));
        h_A(i, j).real() = -sin(auxImJ) - sin(auxI) * sin(auxJ) + cos(auxI) * cos(auxJ);
        h_A(i, j).imag() = -sin(auxImJ) - sin(auxI) * sin(auxJ) - cos(auxI) * cos(auxJ);
      }
    }
  }

  if (_useHermitianOption) {
    for (int i = 0; i < _M; ++i) {
      auxI = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i));
      for (int j = 0; j < _N; ++j) {
        auxJ                    = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(j));
        auxIpJ                  = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i + j));
        h_expected(i, j).real() = -2. * sin(auxI) * sin(auxJ);
        h_expected(i, j).imag() = 2. * (cos(auxIpJ) - sin(auxIpJ));
      }
    }
  } else {
    for (int i = 0; i < _M; ++i) {
      auxI = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i));
      for (int j = 0; j < _N; ++j) {
        auxJ                    = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(j));
        auxImJ                  = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i - j));
        h_expected(i, j).real() = 2. * cos(auxI) * cos(auxJ);
        h_expected(i, j).imag() = -2. * sin(auxImJ);
      }
    }
  }
}

// Code for non-complex values
template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
template <class T>
typename std::enable_if<
    !std::is_same<T, Kokkos::complex<float>>::value && !std::is_same<T, Kokkos::complex<double>>::value, void>::type
GerTester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::populateAnalyticalValues(
    T& alpha, _HostViewTypeX& h_x, _HostViewTypeY& h_y, _HostViewTypeA& h_A, _ViewTypeExpected& h_expected) {
  _AuxType auxI(0.);
  _AuxType auxJ(0.);
  _AuxType auxIpJ(0.);

  alpha = 3;

  for (int i = 0; i < _M; ++i) {
    auxI   = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i));
    h_x[i] = sin(auxI);
  }

  for (int j = 0; j < _N; ++j) {
    auxJ   = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(j));
    h_y[j] = cos(auxJ);
  }

  for (int i = 0; i < _M; ++i) {
    auxI = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i));
    for (int j = 0; j < _N; ++j) {
      auxJ      = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(j));
      h_A(i, j) = 3 * cos(auxI) * sin(auxJ);
    }
  }

  for (int i = 0; i < _M; ++i) {
    for (int j = 0; j < _N; ++j) {
      auxIpJ           = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i + j));
      h_expected(i, j) = 3 * sin(auxIpJ);
    }
  }
}

// Code for complex values
template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
template <class T>
typename std::enable_if<
    std::is_same<T, Kokkos::complex<float>>::value || std::is_same<T, Kokkos::complex<double>>::value, void>::type
GerTester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::populateVanillaValues(
    const T& alpha, const _HostViewTypeX& h_x, const _HostViewTypeY& h_y, const _HostViewTypeA& h_A,
    _ViewTypeExpected& h_vanilla) {
  if (_vanillaUsesDifferentOrderOfOps) {
    if (_useHermitianOption) {
      for (int i = 0; i < _M; ++i) {
        for (int j = 0; j < _N; ++j) {
          h_vanilla(i, j) = h_A(i, j) + alpha * _KAT_A::conj(h_y(j)) * h_x(i);
        }
      }
    } else {
      for (int i = 0; i < _M; ++i) {
        for (int j = 0; j < _N; ++j) {
          h_vanilla(i, j) = h_A(i, j) + alpha * h_y(j) * h_x(i);
        }
      }
    }
  } else {
    if (_useHermitianOption) {
      for (int i = 0; i < _M; ++i) {
        for (int j = 0; j < _N; ++j) {
          h_vanilla(i, j) = h_A(i, j) + alpha * h_x(i) * _KAT_A::conj(h_y(j));
        }
      }
    } else {
      for (int i = 0; i < _M; ++i) {
        for (int j = 0; j < _N; ++j) {
          h_vanilla(i, j) = h_A(i, j) + alpha * h_x(i) * h_y(j);
        }
      }
    }
  }
}

// Code for non-complex values
template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
template <class T>
typename std::enable_if<
    !std::is_same<T, Kokkos::complex<float>>::value && !std::is_same<T, Kokkos::complex<double>>::value, void>::type
GerTester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::populateVanillaValues(
    const T& alpha, const _HostViewTypeX& h_x, const _HostViewTypeY& h_y, const _HostViewTypeA& h_A,
    _ViewTypeExpected& h_vanilla) {
  if (_vanillaUsesDifferentOrderOfOps) {
    for (int i = 0; i < _M; ++i) {
      for (int j = 0; j < _N; ++j) {
        h_vanilla(i, j) = h_A(i, j) + alpha * h_y(j) * h_x(i);
      }
    }
  } else {
    for (int i = 0; i < _M; ++i) {
      for (int j = 0; j < _N; ++j) {
        h_vanilla(i, j) = h_A(i, j) + alpha * h_x(i) * h_y(j);
      }
    }
  }
}

template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
template <class T>
T GerTester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::shrinkAngleToZeroTwoPiRange(
    const T input) {
  T output(input);
#if 0
  T twoPi( 2. * Kokkos::numbers::pi );
  if (input > 0.) {
    output -= std::floor( input / twoPi ) * twoPi;
  }
  else if (input < 0.) {
    output += std::floor( -input / twoPi ) * twoPi;
  }
#endif
  return output;
}

// Code for complex values
template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
template <class T>
typename std::enable_if<
    std::is_same<T, Kokkos::complex<float>>::value || std::is_same<T, Kokkos::complex<double>>::value, void>::type
GerTester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::compareVanillaAgainstExpected(
    const T& alpha, const _ViewTypeExpected& h_vanilla, const _ViewTypeExpected& h_expected) {
  int maxNumErrorsAllowed(static_cast<double>(_M) * static_cast<double>(_N) * 1.e-3);

  if (_useAnalyticalResults) {
    int numErrorsRealAbs(0);
    int numErrorsRealRel(0);
    int numErrorsImagAbs(0);
    int numErrorsImagRel(0);
    _AuxType diff(0.);
    _AuxType diffThreshold(0.);
    bool errorHappened(false);
    _AuxType maxErrorRealRel(0.);
    int iForMaxErrorRealRel(0);
    int jForMaxErrorRealRel(0);
    _AuxType maxErrorImagRel(0.);
    int iForMaxErrorImagRel(0);
    int jForMaxErrorImagRel(0);

    for (int i(0); i < _M; ++i) {
      for (int j(0); j < _N; ++j) {
        diff          = _KAT_A::abs(h_expected(i, j).real() - h_vanilla(i, j).real());
        errorHappened = false;
        if (h_expected(i, j).real() == 0.) {
          diffThreshold = _KAT_A::abs(_absTol);
          if (diff > diffThreshold) {
            errorHappened = true;
            numErrorsRealAbs++;
          }
        } else {
          _AuxType aux = diff / _KAT_A::abs(h_expected(i, j).real());
          if (maxErrorRealRel < aux) {
            maxErrorRealRel     = aux;
            iForMaxErrorRealRel = i;
            jForMaxErrorRealRel = j;
          }

          diffThreshold = _KAT_A::abs(_relTol * h_expected(i, j).real());
          if (diff > diffThreshold) {
            errorHappened = true;
            numErrorsRealRel++;
          }
        }
        if (errorHappened && (numErrorsRealAbs + numErrorsRealRel == 1)) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
          std::cout << "ERROR, i = " << i << ", j = " << j << ": h_expected(i,j).real() = " << h_expected(i, j).real()
                    << ", h_vanilla(i,j).real() = " << h_vanilla(i, j).real()
                    << ", _KAT_A::abs(h_expected(i,j).real() - "
                       "h_vanilla(i,j).real()) = "
                    << diff << ", diffThreshold = " << diffThreshold << std::endl;
#endif
        }

        diff          = _KAT_A::abs(h_expected(i, j).imag() - h_vanilla(i, j).imag());
        errorHappened = false;
        if (h_expected(i, j).imag() == 0.) {
          diffThreshold = _KAT_A::abs(_absTol);
          if (diff > diffThreshold) {
            errorHappened = true;
            numErrorsImagAbs++;
          }
        } else {
          _AuxType aux = diff / _KAT_A::abs(h_expected(i, j).imag());
          if (maxErrorImagRel < aux) {
            maxErrorImagRel     = aux;
            iForMaxErrorImagRel = i;
            jForMaxErrorImagRel = j;
          }

          diffThreshold = _KAT_A::abs(_relTol * h_expected(i, j).imag());
          if (diff > diffThreshold) {
            errorHappened = true;
            numErrorsImagRel++;
          }
        }
        if (errorHappened && (numErrorsImagAbs + numErrorsImagRel == 1)) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
          std::cout << "ERROR, i = " << i << ", j = " << j << ": h_expected(i,j).imag() = " << h_expected(i, j).imag()
                    << ", h_vanilla(i,j).imag() = " << h_vanilla(i, j).imag()
                    << ", _KAT_A::abs(h_expected(i,j).imag() - "
                       "h_vanilla(i,j).imag()) = "
                    << diff << ", diffThreshold = " << diffThreshold << std::endl;
#endif
        }
      }  // for j
    }    // for i
    {
      std::ostringstream msg;
      msg << ", A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr << ", _A_is_ll = " << _A_is_ll
          << ", alpha type = " << typeid(alpha).name() << ", _useHermitianOption = " << _useHermitianOption
          << ": vanilla differs too much from analytical on real components"
          << ", numErrorsRealAbs = " << numErrorsRealAbs << ", numErrorsRealRel = " << numErrorsRealRel
          << ", maxErrorRealRel = " << maxErrorRealRel << ", iForMaxErrorRealRel = " << iForMaxErrorRealRel
          << ", jForMaxErrorRealRel = " << jForMaxErrorRealRel << ", h_expected(i,j).real() = "
          << (((_M > 0) && (_N > 0)) ? h_expected(iForMaxErrorRealRel, jForMaxErrorRealRel).real() : 9.999e+99)
          << ", h_vanilla(i,j).real() = "
          << (((_M > 0) && (_N > 0)) ? h_vanilla(iForMaxErrorRealRel, jForMaxErrorRealRel).real() : 9.999e+99)
          << ", maxNumErrorsAllowed = " << maxNumErrorsAllowed;

      int numErrorsReal(numErrorsRealAbs + numErrorsRealRel);
      if (numErrorsReal > 0) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
        std::cout << "WARNING" << msg.str() << std::endl;
#endif
      }
      EXPECT_LE(numErrorsReal, maxNumErrorsAllowed) << "Failed test" << msg.str();
    }
    {
      std::ostringstream msg;
      msg << ", A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr << ", _A_is_ll = " << _A_is_ll
          << ", alpha type = " << typeid(alpha).name() << ", _useHermitianOption = " << _useHermitianOption
          << ": vanilla differs too much from analytical on imag components"
          << ", numErrorsImagAbs = " << numErrorsImagAbs << ", numErrorsImagRel = " << numErrorsImagRel
          << ", maxErrorImagRel = " << maxErrorImagRel << ", iForMaxErrorImagRel = " << iForMaxErrorImagRel
          << ", jForMaxErrorImagRel = " << jForMaxErrorImagRel << ", h_expected(i,j).imag() = "
          << (((_M > 0) && (_N > 0)) ? h_expected(iForMaxErrorImagRel, jForMaxErrorImagRel).imag() : 9.999e+99)
          << ", h_vanilla(i,j).imag() = "
          << (((_M > 0) && (_N > 0)) ? h_vanilla(iForMaxErrorImagRel, jForMaxErrorImagRel).imag() : 9.999e+99)
          << ", maxNumErrorsAllowed = " << maxNumErrorsAllowed;

      int numErrorsImag(numErrorsImagAbs + numErrorsImagRel);
      if (numErrorsImag > 0) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
        std::cout << "WARNING" << msg.str() << std::endl;
#endif
      }
      EXPECT_LE(numErrorsImag, maxNumErrorsAllowed) << "Failed test" << msg.str();
    }
  } else {
    int numErrorsReal(0);
    int numErrorsImag(0);

    for (int i(0); i < _M; ++i) {
      for (int j(0); j < _N; ++j) {
        if (h_expected(i, j).real() != h_vanilla(i, j).real()) {
          if (numErrorsReal == 0) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
            std::cout << "ERROR, i = " << i << ", j = " << j << ": h_expected(i,j).real() = " << h_expected(i, j).real()
                      << ", h_vanilla(i,j).real() = " << h_vanilla(i, j).real() << std::endl;
#endif
          }
          numErrorsReal++;
        }

        if (h_expected(i, j).imag() != h_vanilla(i, j).imag()) {
          if (numErrorsImag == 0) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
            std::cout << "ERROR, i = " << i << ", j = " << j << ": h_expected(i,j).imag() = " << h_expected(i, j).imag()
                      << ", h_vanilla(i,j).imag() = " << h_vanilla(i, j).imag() << std::endl;
#endif
          }
          numErrorsImag++;
        }
      }  // for j
    }    // for i
    EXPECT_EQ(numErrorsReal, 0) << "Failed test"
                                << ", A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr
                                << ", _A_is_ll = " << _A_is_ll << ", alpha type = " << typeid(alpha).name()
                                << ", _useHermitianOption = " << _useHermitianOption
                                << ": vanilla result is incorrect on real components"
                                << ", numErrorsReal = " << numErrorsReal;
    EXPECT_EQ(numErrorsImag, 0) << "Failed test"
                                << ", A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr
                                << ", _A_is_ll = " << _A_is_ll << ", alpha type = " << typeid(alpha).name()
                                << ", _useHermitianOption = " << _useHermitianOption
                                << ": vanilla result is incorrect on imag components"
                                << ", numErrorsImag = " << numErrorsImag;
  }
}

// Code for non-complex values
template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
template <class T>
typename std::enable_if<
    !std::is_same<T, Kokkos::complex<float>>::value && !std::is_same<T, Kokkos::complex<double>>::value, void>::type
GerTester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::compareVanillaAgainstExpected(
    const T& alpha, const _ViewTypeExpected& h_vanilla, const _ViewTypeExpected& h_expected) {
  int maxNumErrorsAllowed(static_cast<double>(_M) * static_cast<double>(_N) * 1.e-3);

  if (_useAnalyticalResults) {
    int numErrorsAbs(0);
    int numErrorsRel(0);
    _AuxType diff(0.);
    _AuxType diffThreshold(0.);
    bool errorHappened(false);
    _AuxType maxErrorRel(0.);
    int iForMaxErrorRel(0);
    int jForMaxErrorRel(0);

    for (int i(0); i < _M; ++i) {
      for (int j(0); j < _N; ++j) {
        diff          = _KAT_A::abs(h_expected(i, j) - h_vanilla(i, j));
        errorHappened = false;
        if (h_expected(i, j) == 0.) {
          diffThreshold = _KAT_A::abs(_absTol);
          if (diff > diffThreshold) {
            errorHappened = true;
            numErrorsAbs++;
          }
        } else {
          _AuxType aux = diff / _KAT_A::abs(h_expected(i, j));
          if (maxErrorRel < aux) {
            maxErrorRel     = aux;
            iForMaxErrorRel = i;
            jForMaxErrorRel = j;
          }

          diffThreshold = _KAT_A::abs(_relTol * h_expected(i, j));
          if (diff > diffThreshold) {
            errorHappened = true;
            numErrorsRel++;
          }
        }
        if (errorHappened && (numErrorsAbs + numErrorsRel == 1)) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
          std::cout << "ERROR, i = " << i << ", j = " << j << ": h_expected(i,j) = " << h_expected(i, j)
                    << ", h_vanilla(i,j) = " << h_vanilla(i, j)
                    << ", _KAT_A::abs(h_expected(i,j) - h_vanilla(i,j)) = " << diff
                    << ", diffThreshold = " << diffThreshold << std::endl;
#endif
        }
      }  // for j
    }    // for i
    {
      std::ostringstream msg;
      msg << ", A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr << ", _A_is_ll = " << _A_is_ll
          << ", alpha type = " << typeid(alpha).name() << ", _useHermitianOption = " << _useHermitianOption
          << ": vanilla differs too much from expected"
          << ", numErrorsAbs = " << numErrorsAbs << ", numErrorsRel = " << numErrorsRel
          << ", maxErrorRel = " << maxErrorRel << ", iForMaxErrorRel = " << iForMaxErrorRel
          << ", jForMaxErrorRel = " << jForMaxErrorRel << ", h_expected(i,j) = "
          << (((_M > 0) && (_N > 0)) ? h_expected(iForMaxErrorRel, jForMaxErrorRel) : 9.999e+99)
          << ", h_vanilla(i,j) = " << (((_M > 0) && (_N > 0)) ? h_vanilla(iForMaxErrorRel, jForMaxErrorRel) : 9.999e+99)
          << ", maxNumErrorsAllowed = " << maxNumErrorsAllowed;

      int numErrors(numErrorsAbs + numErrorsRel);
      if (numErrors > 0) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
        std::cout << "WARNING" << msg.str() << std::endl;
#endif
      }
      EXPECT_LE(numErrors, maxNumErrorsAllowed) << "Failed test" << msg.str();
    }
  } else {
    int numErrors(0);

    for (int i(0); i < _M; ++i) {
      for (int j(0); j < _N; ++j) {
        if (h_expected(i, j) != h_vanilla(i, j)) {
          if (numErrors == 0) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
            std::cout << "ERROR, i = " << i << ", j = " << j << ": h_expected(i,j) = " << h_expected(i, j)
                      << ", h_vanilla(i,j) = " << h_vanilla(i, j) << std::endl;
#endif
          }
          numErrors++;
        }
      }  // for j
    }    // for i
    EXPECT_EQ(numErrors, 0) << "Failed test"
                            << ", A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr
                            << ", _A_is_ll = " << _A_is_ll << ", alpha type = " << typeid(alpha).name()
                            << ", _useHermitianOption = " << _useHermitianOption << ": vanilla result is incorrect"
                            << ", numErrors = " << numErrors;
  }
}

// Code for complex values
template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
template <class T>
typename std::enable_if<
    std::is_same<T, Kokkos::complex<float>>::value || std::is_same<T, Kokkos::complex<double>>::value, void>::type
GerTester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::compareKkGerAgainstExpected(
    const T& alpha, const _HostViewTypeA& h_A, const _ViewTypeExpected& h_expected) {
  int maxNumErrorsAllowed(static_cast<double>(_M) * static_cast<double>(_N) * 1.e-3);

  int numErrorsRealAbs(0);
  int numErrorsRealRel(0);
  int numErrorsImagAbs(0);
  int numErrorsImagRel(0);
  _AuxType diff(0.);
  _AuxType diffThreshold(0.);
  bool errorHappened(false);
  _AuxType maxErrorRealRel(0.);
  int iForMaxErrorRealRel(0);
  int jForMaxErrorRealRel(0);
  _AuxType maxErrorImagRel(0.);
  int iForMaxErrorImagRel(0);
  int jForMaxErrorImagRel(0);
  for (int i(0); i < _M; ++i) {
    for (int j(0); j < _N; ++j) {
      diff          = _KAT_A::abs(h_expected(i, j).real() - h_A(i, j).real());
      errorHappened = false;
      if (h_expected(i, j).real() == 0.) {
        diffThreshold = _KAT_A::abs(_absTol);
        if (diff > diffThreshold) {
          errorHappened = true;
          numErrorsRealAbs++;
        }
      } else {
        _AuxType aux = diff / _KAT_A::abs(h_expected(i, j).real());
        if (maxErrorRealRel < aux) {
          maxErrorRealRel     = aux;
          iForMaxErrorRealRel = i;
          jForMaxErrorRealRel = j;
        }

        diffThreshold = _KAT_A::abs(_relTol * h_expected(i, j).real());
        if (diff > diffThreshold) {
          errorHappened = true;
          numErrorsRealRel++;
        }
      }
      if (errorHappened && (numErrorsRealAbs + numErrorsRealRel == 1)) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
        std::cout << "ERROR, i = " << i << ", j = " << j << ": h_expected(i,j).real() = " << h_expected(i, j).real()
                  << ", h_A(i,j).real() = " << h_A(i, j).real()
                  << ", _KAT_A::abs(h_expected(i,j).real() - h_A(i,j).real()) = " << diff
                  << ", diffThreshold = " << diffThreshold << std::endl;
#endif
      }

      diff          = _KAT_A::abs(h_expected(i, j).imag() - h_A(i, j).imag());
      errorHappened = false;
      if (h_expected(i, j).imag() == 0.) {
        diffThreshold = _KAT_A::abs(_absTol);
        if (diff > diffThreshold) {
          errorHappened = true;
          numErrorsImagAbs++;
        }
      } else {
        _AuxType aux = diff / _KAT_A::abs(h_expected(i, j).imag());
        if (maxErrorImagRel < aux) {
          maxErrorImagRel     = aux;
          iForMaxErrorImagRel = i;
          jForMaxErrorImagRel = j;
        }

        diffThreshold = _KAT_A::abs(_relTol * h_expected(i, j).imag());
        if (diff > diffThreshold) {
          errorHappened = true;
          numErrorsImagRel++;
        }
      }
      if (errorHappened && (numErrorsImagAbs + numErrorsImagRel == 1)) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
        std::cout << "ERROR, i = " << i << ", j = " << j << ": h_expected(i,j).imag() = " << h_expected(i, j).imag()
                  << ", h_A(i,j).imag() = " << h_A(i, j).imag()
                  << ", _KAT_A::abs(h_expected(i,j).imag() - h_A(i,j).imag()) = " << diff
                  << ", diffThreshold = " << diffThreshold << std::endl;
#endif
      }
    }  // for j
  }    // for i
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr << ", _A_is_ll = " << _A_is_ll
            << ", alpha type = " << typeid(alpha).name() << ", _useHermitianOption = " << _useHermitianOption
            << ", numErrorsRealAbs = " << numErrorsRealAbs << ", numErrorsRealRel = " << numErrorsRealRel
            << ", maxErrorRealRel = " << maxErrorRealRel << ", iForMaxErrorRealRel = " << iForMaxErrorRealRel
            << ", jForMaxErrorRealRel = " << jForMaxErrorRealRel << ", h_expected(i,j).real() = "
            << (((_M > 0) && (_N > 0)) ? h_expected(iForMaxErrorRealRel, jForMaxErrorRealRel).real() : 9.999e+99)
            << ", h_A(i,j).real() = "
            << (((_M > 0) && (_N > 0)) ? h_A(iForMaxErrorRealRel, jForMaxErrorRealRel).real() : 9.999e+99)
            << ", numErrorsImagAbs = " << numErrorsImagAbs << ", numErrorsImagRel = " << numErrorsImagRel
            << ", maxErrorImagRel = " << maxErrorImagRel << ", iForMaxErrorImagRel = " << iForMaxErrorImagRel
            << ", jForMaxErrorImagRel = " << jForMaxErrorImagRel << ", h_expected(i,j).imag() = "
            << (((_M > 0) && (_N > 0)) ? h_expected(iForMaxErrorImagRel, jForMaxErrorImagRel).imag() : 9.999e+99)
            << ", h_A(i,j).imag() = "
            << (((_M > 0) && (_N > 0)) ? h_A(iForMaxErrorImagRel, jForMaxErrorImagRel).imag() : 9.999e+99)
            << ", maxNumErrorsAllowed = " << maxNumErrorsAllowed << std::endl;
  if ((_M == 2131) && (_N == 2131)) {
    std::cout << "Information"
              << ": A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr << ", _A_is_ll = " << _A_is_ll
              << ", alpha type = " << typeid(alpha).name() << ", _useHermitianOption = " << _useHermitianOption
              << ", h_expected(11, 2119) = (" << h_expected(11, 2119).real() << ", " << h_expected(11, 2119).imag()
              << ")"
              << ", h_A(11, 2119) = (" << h_A(11, 2119).real() << ", " << h_A(11, 2119).imag() << ")" << std::endl;
    std::cout << "Information"
              << ": A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr << ", _A_is_ll = " << _A_is_ll
              << ", alpha type = " << typeid(alpha).name() << ", _useHermitianOption = " << _useHermitianOption
              << ", h_expected(710, 1065) = (" << h_expected(710, 1065).real() << ", " << h_expected(710, 1065).imag()
              << ")"
              << ", h_A(710, 1065) = (" << h_A(710, 1065).real() << ", " << h_A(710, 1065).imag() << ")" << std::endl;
  }
#endif
  {
    std::ostringstream msg;
    msg << ", A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr << ", _A_is_ll = " << _A_is_ll
        << ", alpha type = " << typeid(alpha).name() << ", _useHermitianOption = " << _useHermitianOption
        << ": ger result is incorrect on real components"
        << ", numErrorsRealAbs = " << numErrorsRealAbs << ", numErrorsRealRel = " << numErrorsRealRel
        << ", maxErrorRealRel = " << maxErrorRealRel << ", iForMaxErrorRealRel = " << iForMaxErrorRealRel
        << ", jForMaxErrorRealRel = " << jForMaxErrorRealRel << ", h_expected(i,j).real() = "
        << (((_M > 0) && (_N > 0)) ? h_expected(iForMaxErrorRealRel, jForMaxErrorRealRel).real() : 9.999e+99)
        << ", h_A(i,j).real() = "
        << (((_M > 0) && (_N > 0)) ? h_A(iForMaxErrorRealRel, jForMaxErrorRealRel).real() : 9.999e+99)
        << ", maxNumErrorsAllowed = " << maxNumErrorsAllowed;

    int numErrorsReal(numErrorsRealAbs + numErrorsRealRel);
    if (numErrorsReal > 0) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
      std::cout << "WARNING" << msg.str() << std::endl;
#endif
    }
    EXPECT_LE(numErrorsReal, maxNumErrorsAllowed) << "Failed test" << msg.str();
  }
  {
    std::ostringstream msg;
    msg << ", A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr << ", _A_is_ll = " << _A_is_ll
        << ", alpha type = " << typeid(alpha).name() << ", _useHermitianOption = " << _useHermitianOption
        << ": ger result is incorrect on imag components"
        << ", numErrorsImagAbs = " << numErrorsImagAbs << ", numErrorsImagRel = " << numErrorsImagRel
        << ", maxErrorImagRel = " << maxErrorImagRel << ", iForMaxErrorImagRel = " << iForMaxErrorImagRel
        << ", jForMaxErrorImagRel = " << jForMaxErrorImagRel << ", h_expected(i,j).imag() = "
        << (((_M > 0) && (_N > 0)) ? h_expected(iForMaxErrorImagRel, jForMaxErrorImagRel).imag() : 9.999e+99)
        << ", h_A(i,j).imag() = "
        << (((_M > 0) && (_N > 0)) ? h_A(iForMaxErrorImagRel, jForMaxErrorImagRel).imag() : 9.999e+99)
        << ", maxNumErrorsAllowed = " << maxNumErrorsAllowed;

    int numErrorsImag(numErrorsImagAbs + numErrorsImagRel);
    if (numErrorsImag > 0) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
      std::cout << "WARNING" << msg.str() << std::endl;
#endif
    }
    EXPECT_LE(numErrorsImag, maxNumErrorsAllowed) << "Failed test" << msg.str();
  }
}

// Code for non-complex values
template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
template <class T>
typename std::enable_if<
    !std::is_same<T, Kokkos::complex<float>>::value && !std::is_same<T, Kokkos::complex<double>>::value, void>::type
GerTester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::compareKkGerAgainstExpected(
    const T& alpha, const _HostViewTypeA& h_A, const _ViewTypeExpected& h_expected) {
  int maxNumErrorsAllowed(static_cast<double>(_M) * static_cast<double>(_N) * 1.e-3);

  int numErrorsAbs(0);
  int numErrorsRel(0);
  _AuxType diff(0.);
  _AuxType diffThreshold(0.);
  bool errorHappened(false);
  _AuxType maxErrorRel(0.);
  int iForMaxErrorRel(0);
  int jForMaxErrorRel(0);
  for (int i(0); i < _M; ++i) {
    for (int j(0); j < _N; ++j) {
      diff          = _KAT_A::abs(h_expected(i, j) - h_A(i, j));
      errorHappened = false;
      if (h_expected(i, j) == 0.) {
        diffThreshold = _KAT_A::abs(_absTol);
        if (diff > diffThreshold) {
          errorHappened = true;
          numErrorsAbs++;
        }
      } else {
        _AuxType aux = diff / _KAT_A::abs(h_expected(i, j));
        if (maxErrorRel < aux) {
          maxErrorRel     = aux;
          iForMaxErrorRel = i;
          jForMaxErrorRel = j;
        }

        diffThreshold = _KAT_A::abs(_relTol * h_expected(i, j));
        if (diff > diffThreshold) {
          errorHappened = true;
          numErrorsRel++;
        }
      }
      if (errorHappened && (numErrorsAbs + numErrorsRel == 1)) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
        std::cout << "ERROR, i = " << i << ", j = " << j << ": h_expected(i,j) = " << h_expected(i, j)
                  << ", h_A(i,j) = " << h_A(i, j) << ", _KAT_A::abs(h_expected(i,j) - h_A(i,j)) = " << diff
                  << ", diffThreshold = " << diffThreshold << std::endl;
#endif
      }
    }  // for j
  }    // for i
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr << ", _A_is_ll = " << _A_is_ll
            << ", alpha type = " << typeid(alpha).name() << ", _useHermitianOption = " << _useHermitianOption
            << ", numErrorsAbs = " << numErrorsAbs << ", numErrorsRel = " << numErrorsRel
            << ", maxErrorRel = " << maxErrorRel << ", iForMaxErrorRel = " << iForMaxErrorRel
            << ", jForMaxErrorRel = " << jForMaxErrorRel << ", h_expected(i,j) = "
            << (((_M > 0) && (_N > 0)) ? h_expected(iForMaxErrorRel, jForMaxErrorRel) : 9.999e+99)
            << ", h_A(i,j) = " << (((_M > 0) && (_N > 0)) ? h_A(iForMaxErrorRel, jForMaxErrorRel) : 9.999e+99)
            << ", maxNumErrorsAllowed = " << maxNumErrorsAllowed << std::endl;
#endif
  {
    std::ostringstream msg;
    msg << ", A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr << ", _A_is_ll = " << _A_is_ll
        << ", alpha type = " << typeid(alpha).name() << ", _useHermitianOption = " << _useHermitianOption
        << ": ger result is incorrect"
        << ", numErrorsAbs = " << numErrorsAbs << ", numErrorsRel = " << numErrorsRel
        << ", maxErrorRel = " << maxErrorRel << ", iForMaxErrorRel = " << iForMaxErrorRel
        << ", jForMaxErrorRel = " << jForMaxErrorRel
        << ", h_expected(i,j) = " << (((_M > 0) && (_N > 0)) ? h_expected(iForMaxErrorRel, jForMaxErrorRel) : 9.999e+99)
        << ", h_A(i,j) = " << (((_M > 0) && (_N > 0)) ? h_A(iForMaxErrorRel, jForMaxErrorRel) : 9.999e+99)
        << ", maxNumErrorsAllowed = " << maxNumErrorsAllowed;

    int numErrors(numErrorsAbs + numErrorsRel);
    if (numErrors > 0) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
      std::cout << "WARNING" << msg.str() << std::endl;
#endif
    }
    EXPECT_LE(numErrors, maxNumErrorsAllowed) << "Failed test" << msg.str();
  }
}

template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
template <class TX, class TY>
void GerTester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::callKkGerAndCompareAgainstExpected(
    const ScalarA& alpha, TX& x, TY& y, view_stride_adapter<_ViewTypeA, false>& A, const _ViewTypeExpected& h_expected,
    const std::string& situation) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  Kokkos::printf(
      "In Test_Blas2_ger.hpp, right before calling KokkosBlas::ger(): "
      "ViewTypeA = %s, _kkGerShouldThrowException=%d\n",
      typeid(_ViewTypeA).name(), _kkGerShouldThrowException);
#endif
  std::string mode = _useHermitianOption ? "H" : "T";
  bool gotStdException(false);
  bool gotUnknownException(false);
  try {
    KokkosBlas::ger(mode.c_str(), alpha, x, y, A.d_view);
  } catch (const std::exception& e) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
    std::cout << "In Test_Blas2_ger, '" << situation << "': caught exception, e.what() = " << e.what() << std::endl;
#endif
    gotStdException = true;
  } catch (...) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
    std::cout << "In Test_Blas2_ger, '" << situation << "': caught unknown exception" << std::endl;
#endif
    gotUnknownException = true;
  }

  EXPECT_EQ(gotUnknownException, false) << "Failed test, '" << situation
                                        << "': unknown exception should not have happened";

  EXPECT_EQ(gotStdException, _kkGerShouldThrowException)
      << "Failed test, '" << situation << "': kk ger() should" << (_kkGerShouldThrowException ? " " : " not ")
      << "have thrown a std::exception";

  if ((gotStdException == false) && (gotUnknownException == false)) {
    Kokkos::deep_copy(A.h_base, A.d_base);

    this->compareKkGerAgainstExpected(alpha, A.h_view, h_expected);
  }
}

}  // namespace Test

template <class ScalarX, class ScalarY, class ScalarA, class Device>
#ifdef HAVE_KOKKOSKERNELS_DEBUG
int test_ger(const std::string& caseName) {
  Kokkos::printf(
      "+======================================================================="
      "===\n");
  Kokkos::printf("Starting %s, device = %s ...\n", caseName.c_str(), typeid(Device).name());
#else
int test_ger(const std::string& /*caseName*/) {
#endif
  bool xBool = std::is_same<ScalarX, float>::value || std::is_same<ScalarX, double>::value ||
               std::is_same<ScalarX, Kokkos::complex<float>>::value ||
               std::is_same<ScalarX, Kokkos::complex<double>>::value;
  bool yBool = std::is_same<ScalarY, float>::value || std::is_same<ScalarY, double>::value ||
               std::is_same<ScalarY, Kokkos::complex<float>>::value ||
               std::is_same<ScalarY, Kokkos::complex<double>>::value;
  bool aBool = std::is_same<ScalarA, float>::value || std::is_same<ScalarA, double>::value ||
               std::is_same<ScalarA, Kokkos::complex<float>>::value ||
               std::is_same<ScalarA, Kokkos::complex<double>>::value;
  bool useAnalyticalResults = xBool && yBool && aBool;

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  Kokkos::printf(
      "+-----------------------------------------------------------------------"
      "---\n");
  Kokkos::printf("Starting %s for LAYOUTLEFT ...\n", caseName.c_str());
#endif
  if (true) {
    Test::GerTester<ScalarX, Kokkos::LayoutLeft, ScalarY, Kokkos::LayoutLeft, ScalarA, Kokkos::LayoutLeft, Device>
        tester;
    tester.test(0, 13, 0);
    tester.test(1024, 0, 0);
    tester.test(1, 1, 0);
    tester.test(2, 2, 0);
    tester.test(1, 2, 0);
    tester.test(13, 13, 0);
    tester.test(13, 1024, 0);
    if (useAnalyticalResults) {
      tester.test(13, 1024, 0, true, false);
      tester.test(13, 1024, 0, true, true);
    } else {
      tester.test(13, 1024, 0, false, true);
    }
    tester.test(50, 40, 4);
    tester.test(1024, 1024, 0);
    tester.test(2131, 2131, 0);
    if (useAnalyticalResults) {
      tester.test(2131, 2131, 0, true, false);
      tester.test(2131, 2131, 0, true, true);
    } else {
      tester.test(2131, 2131, 0, false, true);
    }
  }

#ifdef HAVE_KOKKOSKERNELS_DEBUG
  Kokkos::printf("Finished %s for LAYOUTLEFT\n", caseName.c_str());
  Kokkos::printf(
      "+-----------------------------------------------------------------------"
      "---\n");
#endif
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  Kokkos::printf(
      "+-----------------------------------------------------------------------"
      "---\n");
  Kokkos::printf("Starting %s for LAYOUTRIGHT ...\n", caseName.c_str());
#endif
  if (true) {
    Test::GerTester<ScalarX, Kokkos::LayoutRight, ScalarY, Kokkos::LayoutRight, ScalarA, Kokkos::LayoutRight, Device>
        tester;
    tester.test(0, 13, 0);
    tester.test(1024, 0, 0);
    tester.test(1, 1, 0);
    tester.test(2, 2, 0);
    tester.test(1, 2, 0);
    tester.test(13, 13, 0);
    tester.test(13, 1024, 0);
    if (useAnalyticalResults) {
      tester.test(13, 1024, 0, true, false);
      tester.test(13, 1024, 0, true, true);
    } else {
      tester.test(13, 1024, 0, false, true);
    }
    tester.test(50, 40, 4);
    tester.test(1024, 1024, 0);
    tester.test(2131, 2131, 0);
    if (useAnalyticalResults) {
      tester.test(2131, 2131, 0, true, false);
      tester.test(2131, 2131, 0, true, true);
    } else {
      tester.test(2131, 2131, 0, false, true);
    }
  }

#ifdef HAVE_KOKKOSKERNELS_DEBUG
  Kokkos::printf("Finished %s for LAYOUTRIGHT\n", caseName.c_str());
  Kokkos::printf(
      "+-----------------------------------------------------------------------"
      "---\n");
#endif
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  Kokkos::printf(
      "+-----------------------------------------------------------------------"
      "---\n");
  Kokkos::printf("Starting %s for LAYOUTSTRIDE ...\n", caseName.c_str());
#endif
  if (true) {
    Test::GerTester<ScalarX, Kokkos::LayoutStride, ScalarY, Kokkos::LayoutStride, ScalarA, Kokkos::LayoutStride, Device>
        tester;
    tester.test(0, 13, 0);
    tester.test(1024, 0, 0);
    tester.test(13, 13, 0);
    tester.test(13, 1024, 0);
    if (useAnalyticalResults) {
      tester.test(13, 1024, 0, true, false);
      tester.test(13, 1024, 0, true, true);
    } else {
      tester.test(13, 1024, 0, false, true);
    }
    tester.test(50, 40, 4);
    tester.test(1024, 1024, 0);
    tester.test(2131, 2131, 0);
    if (useAnalyticalResults) {
      tester.test(2131, 2131, 0, true, false);
      tester.test(2131, 2131, 0, true, true);
    } else {
      tester.test(2131, 2131, 0, false, true);
    }
  }

#ifdef HAVE_KOKKOSKERNELS_DEBUG
  Kokkos::printf("Finished %s for LAYOUTSTRIDE\n", caseName.c_str());
  Kokkos::printf(
      "+-----------------------------------------------------------------------"
      "---\n");
#endif
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  Kokkos::printf(
      "+-----------------------------------------------------------------------"
      "---\n");
  Kokkos::printf("Starting %s for MIXED LAYOUTS ...\n", caseName.c_str());
#endif
  if (true) {
    Test::GerTester<ScalarX, Kokkos::LayoutStride, ScalarY, Kokkos::LayoutLeft, ScalarA, Kokkos::LayoutRight, Device>
        tester;
    tester.test(1024, 1024, 0);
    if (useAnalyticalResults) {
      tester.test(1024, 1024, 0, true, false);
      tester.test(1024, 1024, 0, true, true);
    } else {
      tester.test(1024, 1024, 0, false, true);
    }
  }

  if (true) {
    Test::GerTester<ScalarX, Kokkos::LayoutLeft, ScalarY, Kokkos::LayoutStride, ScalarA, Kokkos::LayoutRight, Device>
        tester;
    tester.test(1024, 1024, 0);
  }

#ifdef HAVE_KOKKOSKERNELS_DEBUG
  Kokkos::printf("Finished %s for MIXED LAYOUTS\n", caseName.c_str());
  Kokkos::printf(
      "+-----------------------------------------------------------------------"
      "---\n");
#endif
#endif

#ifdef HAVE_KOKKOSKERNELS_DEBUG
  Kokkos::printf("Finished %s\n", caseName.c_str());
  Kokkos::printf(
      "+======================================================================="
      "===\n");
#endif
  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, ger_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::ger_float");
  test_ger<float, float, float, TestDevice>("test case ger_float");
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, ger_complex_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::ger_complex_float");
  test_ger<Kokkos::complex<float>, Kokkos::complex<float>, Kokkos::complex<float>, TestDevice>(
      "test case ger_complex_float");
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, ger_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::ger_double");
  test_ger<double, double, double, TestDevice>("test case ger_double");
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, ger_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::ger_complex_double");
  test_ger<Kokkos::complex<double>, Kokkos::complex<double>, Kokkos::complex<double>, TestDevice>(
      "test case ger_complex_double");
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, ger_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::ger_int");
  test_ger<int, int, int, TestDevice>("test case ger_int");
  Kokkos::Profiling::popRegion();
}
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
TEST_F(TestCategory, ger_double_int_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::ger_double_int_float");
  test_ger<double, int, float, TestDevice>("test case ger_double_int_float");
  Kokkos::Profiling::popRegion();
}
#endif
