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
// the operations:
// --> A += alpha * x * y^T + alpha * y * x^T, or
// --> A += alpha * x * y^H + conj(alpha) * y * x^H
// 01) Type of 'x' components: float, double, complex, ...
// 02) Type of 'y' components: float, double, complex, ...
// 03) Type of 'A' components: float, double, complex, ...
// 04) Execution space: serial, threads, OpenMP, Cuda, ...
// 05) Layout of 'x'
// 06) Layout of 'y'
// 07) Layout of 'A'
// 08) Dimension of 'A'
// 09) Options 'const' or 'non const' for x view, when calling syr2()
// 10) Options 'const' or 'non const' for y view, when calling syr2()
// 11) Usage of analytical results in the tests
// 12) Options 'T' or 'H' when calling syr2()
// 13) Options 'U' or 'L' when calling syr2()
//
// Choices (01)-(05) are selected in the routines TEST_F() at the
// very bottom of the file, when calling test_syr2<...>().
//
// Choices (06)-(13) are selected in routine test_syr2<...>(),
// when calling the method test() of class Test::Syr2Tester<...>.
//
// The class Test::Syr2Tester<...> represents the "core" of the test
// logic, where all calculations, comparisons, and success/failure
// decisions are performed.
//
// A high level explanation of method Test::SyrTester<...>::test()
// is given by the 7 steps named "Step 1 of 7" to "Step 7 of 7"
// in the code below.
// **********************************************************************

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBlas2_syr2.hpp>
#include <Kokkos_MathematicalConstants.hpp>

namespace Test {

template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
class Syr2Tester {
 public:
  Syr2Tester();

  ~Syr2Tester();

  void test(const int N, const int nonConstConstCombinations, const bool useAnalyticalResults = false,
            const bool useHermitianOption = false, const bool useUpOption = false);

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
  compareKkSyr2AgainstReference(const T& alpha, const _HostViewTypeA& h_A, const _ViewTypeExpected& h_reference);

  template <class T>
  typename std::enable_if<
      !std::is_same<T, Kokkos::complex<float>>::value && !std::is_same<T, Kokkos::complex<double>>::value, void>::type
  compareKkSyr2AgainstReference(const T& alpha, const _HostViewTypeA& h_A, const _ViewTypeExpected& h_reference);

  template <class T>
  T shrinkAngleToZeroTwoPiRange(const T input);

  template <class TX, class TY>
  void callKkSyr2AndCompareAgainstExpected(const ScalarA& alpha, TX& x, TY& y,
                                           view_stride_adapter<_ViewTypeA, false>& A,
                                           const _ViewTypeExpected& h_expected, const std::string& situation);

  template <class TX, class TY>
  void callKkGerAndCompareKkSyr2AgainstIt(const ScalarA& alpha, TX& x, TY& y,
                                          view_stride_adapter<_ViewTypeA, false>& org_A, const _HostViewTypeA& h_A_syr2,
                                          const std::string& situation);

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
  bool _useUpOption;
  bool _kkSyr2ShouldThrowException;
  bool _kkGerShouldThrowException;
};

template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
Syr2Tester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::Syr2Tester()
    : _A_is_complex(std::is_same<ScalarA, Kokkos::complex<float>>::value ||
                    std::is_same<ScalarA, Kokkos::complex<double>>::value),
      _A_is_lr(std::is_same<tLayoutA, Kokkos::LayoutRight>::value),
      _A_is_ll(std::is_same<tLayoutA, Kokkos::LayoutLeft>::value),
      _testIsGpu(KokkosKernels::Impl::kk_is_gpu_exec_space<typename Device::execution_space>())
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
      ,
      _vanillaUsesDifferentOrderOfOps(_A_is_lr)
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
      _useUpOption(false),
      _kkSyr2ShouldThrowException(false),
      _kkGerShouldThrowException(false) {
}

template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
Syr2Tester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::~Syr2Tester() {
  // Nothing to do
}

template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
void Syr2Tester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::test(
    const int N, const int nonConstConstCombinations, const bool useAnalyticalResults, const bool useHermitianOption,
    const bool useUpOption) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Entering Syr2Tester::test()... - - - - - - - - - - - - - - - - "
               "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
               "- - - - - - - - - "
            << std::endl;

  std::cout << "_A_is_complex = " << _A_is_complex << ", _A_is_lr = " << _A_is_lr << ", _A_is_ll = " << _A_is_ll
            << ", _testIsGpu = " << _testIsGpu
            << ", _vanillaUsesDifferentOrderOfOps = " << _vanillaUsesDifferentOrderOfOps << ", _absTol = " << _absTol
            << ", _relTol = " << _relTol << ", nonConstConstCombinations = " << nonConstConstCombinations
            << ", useAnalyticalResults = " << useAnalyticalResults << ", useHermitianOption = " << useHermitianOption
            << ", useUpOption = " << useUpOption << std::endl;
#endif
  // ********************************************************************
  // Step 1 of 7: declare main types and variables
  // ********************************************************************
  _M                    = N;
  _N                    = N;
  _useAnalyticalResults = useAnalyticalResults;
  _useHermitianOption   = useHermitianOption;
  _useUpOption          = useUpOption;

#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
  _kkSyr2ShouldThrowException = false;

  _kkGerShouldThrowException = false;
  if (_A_is_complex && _useHermitianOption) {
    _kkGerShouldThrowException = !_A_is_ll;
  }
#endif

  bool test_x(false);
  bool test_cx(false);
  if (nonConstConstCombinations == 0) {
    test_x = true;
  } else if (nonConstConstCombinations == 1) {
    test_cx = true;
  } else {
    test_x  = true;
    test_cx = true;
  }

  view_stride_adapter<_ViewTypeX, false> x("X", _M);
  view_stride_adapter<_ViewTypeY, false> y("Y", _N);
  view_stride_adapter<_ViewTypeA, false> A("A", _M, _N);

  view_stride_adapter<_ViewTypeExpected, true> h_expected("expected A += alpha * x * x^{t,h}", _M, _N);
  bool expectedResultIsKnown = false;

  using AlphaCoeffType = typename _ViewTypeA::non_const_value_type;
  ScalarA alpha(Kokkos::ArithTraits<AlphaCoeffType>::zero());

  // ********************************************************************
  // Step 2 of 7: populate alpha, h_x, h_A, h_expected, x, A
  // ********************************************************************
  this->populateVariables(alpha, x, y, A, h_expected.d_view, expectedResultIsKnown);

  // ********************************************************************
  // Step 3 of 7: populate h_vanilla
  // ********************************************************************
  view_stride_adapter<_ViewTypeExpected, true> h_vanilla("vanilla = A + alpha * x * x^{t,h}", _M, _N);
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "In Test_Blas2_syr2.hpp, computing vanilla A with alpha type = " << typeid(alpha).name() << std::endl;
#endif
  this->populateVanillaValues(alpha, x.h_view, y.h_view, A.h_view, h_vanilla.d_view);

  // ********************************************************************
  // Step 4 of 7: use h_vanilla and h_expected as appropriate
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
  // Step 5 of 7: test with 'non const x'
  // ********************************************************************
  view_stride_adapter<_ViewTypeA, false> org_A("Org_A", _M, _N);
  Kokkos::deep_copy(org_A.d_base, A.d_base);
  Kokkos::deep_copy(org_A.h_view, A.h_view);

  if (test_x) {
    this->callKkSyr2AndCompareAgainstExpected(alpha, x.d_view, y.d_view, A, h_expected.d_view, "non const x");

    if ((_useAnalyticalResults == false) &&  // Just to save run time
        (_kkGerShouldThrowException == false)) {
      this->callKkGerAndCompareKkSyr2AgainstIt(alpha, x.d_view, y.d_view, org_A, A.h_view, "non const x");
    }
  }

  // ********************************************************************
  // Step 6 of 7: test with const x
  // ********************************************************************
  if (test_cx) {
    Kokkos::deep_copy(A.d_base, org_A.d_base);

    this->callKkSyr2AndCompareAgainstExpected(alpha, x.d_view_const, y.d_view_const, A, h_expected.d_view, "const x");
  }

  // ********************************************************************
  // Step 7 of 7: tests with invalid values on the first input parameter
  // ********************************************************************
  EXPECT_ANY_THROW(KokkosBlas::syr2(".", "U", alpha, x.d_view, y.d_view, A.d_view))
      << "Failed test: kk syr2 should have thrown an exception for mode '.'";
  EXPECT_ANY_THROW(KokkosBlas::syr2("", "U", alpha, x.d_view, y.d_view, A.d_view))
      << "Failed test: kk syr2 should have thrown an exception for mode ''";
  EXPECT_ANY_THROW(KokkosBlas::syr2("T", ".", alpha, x.d_view, y.d_view, A.d_view))
      << "Failed test: kk syr2 should have thrown an exception for uplo '.'";
  EXPECT_ANY_THROW(KokkosBlas::syr2("T", "", alpha, x.d_view, y.d_view, A.d_view))
      << "Failed test: kk syr2 should have thrown an exception for uplo ''";

#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Leaving Syr2Tester::test() - - - - - - - - - - - - - - - - - - "
               "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
               "- - - - - - - "
            << std::endl;
#endif
}

template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
void Syr2Tester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::populateVariables(
    ScalarA& alpha, view_stride_adapter<_ViewTypeX, false>& x, view_stride_adapter<_ViewTypeY, false>& y,
    view_stride_adapter<_ViewTypeA, false>& A, _ViewTypeExpected& h_expected, bool& expectedResultIsKnown) {
  expectedResultIsKnown = false;

  if (_useAnalyticalResults) {
    this->populateAnalyticalValues(alpha, x.h_view, y.h_view, A.h_view, h_expected);
    Kokkos::deep_copy(x.d_base, x.h_base);
    Kokkos::deep_copy(y.d_base, y.h_base);
    Kokkos::deep_copy(A.d_base, A.h_base);

    expectedResultIsKnown = true;
  } else if (_N == 1) {
    alpha = 3;

    x.h_view[0] = 2;

    y.h_view[0] = 4;

    A.h_view(0, 0) = 7;

    Kokkos::deep_copy(x.d_base, x.h_base);
    Kokkos::deep_copy(y.d_base, y.h_base);
    Kokkos::deep_copy(A.d_base, A.h_base);

    h_expected(0, 0)      = 55;
    expectedResultIsKnown = true;
  } else if (_N == 2) {
    alpha = 3;

    x.h_view[0] = -2;
    x.h_view[1] = 9;

    y.h_view[0] = 5;
    y.h_view[1] = -4;

    A.h_view(0, 0) = 17;
    A.h_view(0, 1) = -43;
    A.h_view(1, 0) = -43;
    A.h_view(1, 1) = 101;

    Kokkos::deep_copy(x.d_base, x.h_base);
    Kokkos::deep_copy(y.d_base, y.h_base);
    Kokkos::deep_copy(A.d_base, A.h_base);

    if (_useUpOption) {
      h_expected(0, 0) = -43;
      h_expected(0, 1) = 116;
      h_expected(1, 0) = -43;
      h_expected(1, 1) = -115;
    } else {
      h_expected(0, 0) = -43;
      h_expected(0, 1) = -43;
      h_expected(1, 0) = 116;
      h_expected(1, 1) = -115;
    }
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

    if (_useHermitianOption && _A_is_complex) {
      // ****************************************************************
      // Make h_A Hermitian
      // ****************************************************************
      for (int i(0); i < _N; ++i) {
        for (int j(i + 1); j < _N; ++j) {
          A.h_view(i, j) = _KAT_A::conj(A.h_view(j, i));
        }
      }

      for (int i(0); i < _N; ++i) {
        A.h_view(i, i) = 0.5 * (A.h_view(i, i) + _KAT_A::conj(A.h_view(i, i)));
      }
    } else {
      // ****************************************************************
      // Make h_A symmetric
      // ****************************************************************
      for (int i(0); i < _N; ++i) {
        for (int j(i + 1); j < _N; ++j) {
          A.h_view(i, j) = A.h_view(j, i);
        }
      }
    }
    Kokkos::deep_copy(A.d_base, A.h_base);
  }

#ifdef HAVE_KOKKOSKERNELS_DEBUG
  if (_N <= 2) {
    for (int i(0); i < _M; ++i) {
      for (int j(0); j < _N; ++j) {
        std::cout << "h_origA(" << i << "," << j << ") = " << A.h_view(i, j) << std::endl;
      }
    }
  }
#endif
}

// Code for complex values
template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
template <class T>
typename std::enable_if<
    std::is_same<T, Kokkos::complex<float>>::value || std::is_same<T, Kokkos::complex<double>>::value, void>::type
Syr2Tester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::populateAnalyticalValues(
    T& alpha, _HostViewTypeX& h_x, _HostViewTypeY& h_y, _HostViewTypeA& h_A, _ViewTypeExpected& h_expected) {
  alpha.real() = 1.4;
  alpha.imag() = -2.3;

  for (int i = 0; i < _M; ++i) {
    _AuxType auxI = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i));
    h_x[i].real() = sin(auxI);
    h_x[i].imag() = sin(auxI);
  }

  for (int i = 0; i < _M; ++i) {
    _AuxType auxI = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i));
    h_y[i].real() = cos(auxI);
    h_y[i].imag() = cos(auxI);
  }

  if (_useHermitianOption) {
    for (int i = 0; i < _M; ++i) {
      for (int j = 0; j < _N; ++j) {
        _AuxType auxIpJ = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i + j));
        _AuxType auxImJ = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i - j));
        if (((_useUpOption == true) && (i <= j)) || ((_useUpOption == false) && (i >= j))) {
          h_A(i, j).real() = sin(auxIpJ);
          h_A(i, j).imag() = -sin(auxImJ);
        } else {
          h_A(i, j).real() = sin(auxIpJ);
          h_A(i, j).imag() = sin(auxImJ);
        }
      }
    }
  } else {
    for (int i = 0; i < _M; ++i) {
      for (int j = 0; j < _N; ++j) {
        _AuxType auxIpJ  = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i + j));
        h_A(i, j).real() = sin(auxIpJ);
        h_A(i, j).imag() = sin(auxIpJ);
      }
    }
  }

  if (_useHermitianOption) {
    for (int i = 0; i < _M; ++i) {
      for (int j = 0; j < _N; ++j) {
        if (((_useUpOption == true) && (i <= j)) || ((_useUpOption == false) && (i >= j))) {
          _AuxType auxIpJ         = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i + j));
          _AuxType auxImJ         = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i - j));
          h_expected(i, j).real() = 3.8 * sin(auxIpJ);
          h_expected(i, j).imag() = -5.6 * sin(auxImJ);
        } else {
          h_expected(i, j).real() = h_A(i, j).real();
          h_expected(i, j).imag() = h_A(i, j).imag();
        }
      }
    }
  } else {
    for (int i = 0; i < _M; ++i) {
      for (int j = 0; j < _N; ++j) {
        if (((_useUpOption == true) && (i <= j)) || ((_useUpOption == false) && (i >= j))) {
          _AuxType auxIpJ         = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i + j));
          h_expected(i, j).real() = 5.6 * sin(auxIpJ);
          h_expected(i, j).imag() = 3.8 * sin(auxIpJ);
        } else {
          h_expected(i, j).real() = h_A(i, j).real();
          h_expected(i, j).imag() = h_A(i, j).imag();
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
Syr2Tester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::populateAnalyticalValues(
    T& alpha, _HostViewTypeX& h_x, _HostViewTypeY& h_y, _HostViewTypeA& h_A, _ViewTypeExpected& h_expected) {
  alpha = std::is_same<_AuxType, int>::value ? 1 : 1.1;

  for (int i = 0; i < _M; ++i) {
    _AuxType auxI = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i));
    h_x[i]        = sin(auxI);
  }

  for (int i = 0; i < _M; ++i) {
    _AuxType auxI = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i));
    h_y[i]        = cos(auxI);
  }

  for (int i = 0; i < _M; ++i) {
    for (int j = 0; j < _N; ++j) {
      _AuxType auxIpJ = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i + j));
      h_A(i, j)       = .1 * sin(auxIpJ);
    }
  }

  for (int i = 0; i < _M; ++i) {
    for (int j = 0; j < _N; ++j) {
      if (((_useUpOption == true) && (i <= j)) || ((_useUpOption == false) && (i >= j))) {
        _AuxType auxIpJ  = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i + j));
        h_expected(i, j) = 1.2 * sin(auxIpJ);
      } else {
        h_expected(i, j) = h_A(i, j);
      }
    }
  }
}

// Code for complex values
template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
template <class T>
typename std::enable_if<
    std::is_same<T, Kokkos::complex<float>>::value || std::is_same<T, Kokkos::complex<double>>::value, void>::type
Syr2Tester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::populateVanillaValues(
    const T& alpha, const _HostViewTypeX& h_x, const _HostViewTypeY& h_y, const _HostViewTypeA& h_A,
    _ViewTypeExpected& h_vanilla) {
  if (_vanillaUsesDifferentOrderOfOps) {
    if (_useHermitianOption) {
      for (int i = 0; i < _M; ++i) {
        for (int j = 0; j < _N; ++j) {
          if (((_useUpOption == true) && (i <= j)) || ((_useUpOption == false) && (i >= j))) {
            h_vanilla(i, j) =
                h_A(i, j) + alpha * _KAT_A::conj(h_y(j)) * h_x(i) + _KAT_A::conj(alpha) * _KAT_A::conj(h_x(j)) * h_y(i);
          } else {
            h_vanilla(i, j) = h_A(i, j);
          }
        }
      }
      for (int i = 0; i < _N; ++i) {
        h_vanilla(i, i).imag() = 0.;
      }
    } else {
      for (int i = 0; i < _M; ++i) {
        for (int j = 0; j < _N; ++j) {
          if (((_useUpOption == true) && (i <= j)) || ((_useUpOption == false) && (i >= j))) {
            h_vanilla(i, j) = h_A(i, j) + alpha * h_x(j) * h_y(i) + alpha * h_y(j) * h_x(i);
          } else {
            h_vanilla(i, j) = h_A(i, j);
          }
        }
      }
    }
  } else {
    if (_useHermitianOption) {
      for (int i = 0; i < _M; ++i) {
        for (int j = 0; j < _N; ++j) {
          if (((_useUpOption == true) && (i <= j)) || ((_useUpOption == false) && (i >= j))) {
            h_vanilla(i, j) =
                h_A(i, j) + alpha * h_x(i) * _KAT_A::conj(h_y(j)) + _KAT_A::conj(alpha) * h_y(i) * _KAT_A::conj(h_x(j));
          } else {
            h_vanilla(i, j) = h_A(i, j);
          }
        }
      }
      for (int i = 0; i < _N; ++i) {
        h_vanilla(i, i).imag() = 0.;
      }
    } else {
      for (int i = 0; i < _M; ++i) {
        for (int j = 0; j < _N; ++j) {
          if (((_useUpOption == true) && (i <= j)) || ((_useUpOption == false) && (i >= j))) {
            h_vanilla(i, j) = h_A(i, j) + alpha * h_x(i) * h_y(j) + alpha * h_y(i) * h_x(j);
          } else {
            h_vanilla(i, j) = h_A(i, j);
          }
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
Syr2Tester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::populateVanillaValues(
    const T& alpha, const _HostViewTypeX& h_x, const _HostViewTypeY& h_y, const _HostViewTypeA& h_A,
    _ViewTypeExpected& h_vanilla) {
  if (_useHermitianOption) {
    if (_vanillaUsesDifferentOrderOfOps) {
      for (int i = 0; i < _M; ++i) {
        for (int j = 0; j < _N; ++j) {
          if (((_useUpOption == true) && (i <= j)) || ((_useUpOption == false) && (i >= j))) {
            h_vanilla(i, j) =
                h_A(i, j) + alpha * h_x(j) * _KAT_A::conj(h_y(i)) + _KAT_A::conj(alpha) * h_y(j) * _KAT_A::conj(h_x(i));
          } else {
            h_vanilla(i, j) = h_A(i, j);
          }
        }
      }
    } else {
      for (int i = 0; i < _M; ++i) {
        for (int j = 0; j < _N; ++j) {
          if (((_useUpOption == true) && (i <= j)) || ((_useUpOption == false) && (i >= j))) {
            h_vanilla(i, j) =
                h_A(i, j) + alpha * h_x(i) * _KAT_A::conj(h_y(j)) + _KAT_A::conj(alpha) * h_y(i) * _KAT_A::conj(h_x(j));
          } else {
            h_vanilla(i, j) = h_A(i, j);
          }
        }
      }
    }
  } else {
    if (_vanillaUsesDifferentOrderOfOps) {
      for (int i = 0; i < _M; ++i) {
        for (int j = 0; j < _N; ++j) {
          if (((_useUpOption == true) && (i <= j)) || ((_useUpOption == false) && (i >= j))) {
            h_vanilla(i, j) = h_A(i, j) + alpha * h_x(j) * h_y(i) + alpha * h_y(j) * h_x(i);
          } else {
            h_vanilla(i, j) = h_A(i, j);
          }
        }
      }
    } else {
      for (int i = 0; i < _M; ++i) {
        for (int j = 0; j < _N; ++j) {
          if (((_useUpOption == true) && (i <= j)) || ((_useUpOption == false) && (i >= j))) {
            h_vanilla(i, j) = h_A(i, j) + alpha * h_x(i) * h_y(j) + alpha * h_y(i) * h_x(j);
          } else {
            h_vanilla(i, j) = h_A(i, j);
          }
        }
      }
    }
  }
}

template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
template <class T>
T Syr2Tester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::shrinkAngleToZeroTwoPiRange(
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
Syr2Tester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::compareVanillaAgainstExpected(
    const T& alpha, const _ViewTypeExpected& h_vanilla, const _ViewTypeExpected& h_expected) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  if (_N <= 2) {
    for (int i(0); i < _M; ++i) {
      for (int j(0); j < _N; ++j) {
        std::cout << "h_exp(" << i << "," << j << ") = " << h_expected(i, j) << ", h_van(" << i << "," << j
                  << ") = " << h_vanilla(i, j) << std::endl;
      }
    }
  }
#endif
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
          << ", _useUpOption = " << _useUpOption << ": vanilla differs too much from analytical on real components"
          << ", numErrorsRealAbs = " << numErrorsRealAbs << ", numErrorsRealRel = " << numErrorsRealRel
          << ", maxErrorRealRel = " << maxErrorRealRel << ", iForMaxErrorRealRel = " << iForMaxErrorRealRel
          << ", jForMaxErrorRealRel = " << jForMaxErrorRealRel << ", h_expected(i,j).real() = "
          << (((_M > 0) && (_N > 0)) ? h_expected(iForMaxErrorRealRel, jForMaxErrorRealRel).real() : 9.999e+99)
          << ", h_vanilla(i,j).real() = "
          << (((_M > 0) && (_N > 0)) ? h_vanilla(iForMaxErrorRealRel, jForMaxErrorRealRel).real() : 9.999e+99)
          << ", maxNumErrorsAllowed = " << maxNumErrorsAllowed;

      int numErrorsReal(numErrorsRealAbs + numErrorsRealRel);
#ifdef HAVE_KOKKOSKERNELS_DEBUG
      if (numErrorsReal > 0) {
        std::cout << "WARNING" << msg.str() << std::endl;
      }
#endif
      EXPECT_LE(numErrorsReal, maxNumErrorsAllowed) << "Failed test" << msg.str();
    }
    {
      std::ostringstream msg;
      msg << ", A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr << ", _A_is_ll = " << _A_is_ll
          << ", alpha type = " << typeid(alpha).name() << ", _useHermitianOption = " << _useHermitianOption
          << ", _useUpOption = " << _useUpOption << ": vanilla differs too much from analytical on imag components"
          << ", numErrorsImagAbs = " << numErrorsImagAbs << ", numErrorsImagRel = " << numErrorsImagRel
          << ", maxErrorImagRel = " << maxErrorImagRel << ", iForMaxErrorImagRel = " << iForMaxErrorImagRel
          << ", jForMaxErrorImagRel = " << jForMaxErrorImagRel << ", h_expected(i,j).imag() = "
          << (((_M > 0) && (_N > 0)) ? h_expected(iForMaxErrorImagRel, jForMaxErrorImagRel).imag() : 9.999e+99)
          << ", h_vanilla(i,j).imag() = "
          << (((_M > 0) && (_N > 0)) ? h_vanilla(iForMaxErrorImagRel, jForMaxErrorImagRel).imag() : 9.999e+99)
          << ", maxNumErrorsAllowed = " << maxNumErrorsAllowed;

      int numErrorsImag(numErrorsImagAbs + numErrorsImagRel);
#ifdef HAVE_KOKKOSKERNELS_DEBUG
      if (numErrorsImag > 0) {
        std::cout << "WARNING" << msg.str() << std::endl;
      }
#endif
      EXPECT_LE(numErrorsImag, maxNumErrorsAllowed) << "Failed test" << msg.str();
    }
  } else {
    int numErrorsReal(0);
    int numErrorsImag(0);

    for (int i(0); i < _M; ++i) {
      for (int j(0); j < _N; ++j) {
        if (h_expected(i, j).real() != h_vanilla(i, j).real()) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
          if (numErrorsReal == 0) {
            std::cout << "ERROR, i = " << i << ", j = " << j << ": h_expected(i,j).real() = " << h_expected(i, j).real()
                      << ", h_vanilla(i,j).real() = " << h_vanilla(i, j).real() << std::endl;
          }
#endif
          numErrorsReal++;
        }

        if (h_expected(i, j).imag() != h_vanilla(i, j).imag()) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
          if (numErrorsImag == 0) {
            std::cout << "ERROR, i = " << i << ", j = " << j << ": h_expected(i,j).imag() = " << h_expected(i, j).imag()
                      << ", h_vanilla(i,j).imag() = " << h_vanilla(i, j).imag() << std::endl;
          }
#endif
          numErrorsImag++;
        }
      }  // for j
    }    // for i
    EXPECT_EQ(numErrorsReal, 0) << "Failed test"
                                << ", A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr
                                << ", _A_is_ll = " << _A_is_ll << ", alpha type = " << typeid(alpha).name()
                                << ", _useHermitianOption = " << _useHermitianOption
                                << ", _useUpOption = " << _useUpOption
                                << ": vanilla result is incorrect on real components"
                                << ", numErrorsReal = " << numErrorsReal;
    EXPECT_EQ(numErrorsImag, 0) << "Failed test"
                                << ", A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr
                                << ", _A_is_ll = " << _A_is_ll << ", alpha type = " << typeid(alpha).name()
                                << ", _useHermitianOption = " << _useHermitianOption
                                << ", _useUpOption = " << _useUpOption
                                << ": vanilla result is incorrect on imag components"
                                << ", numErrorsImag = " << numErrorsImag;
  }
}

// Code for non-complex values
template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
template <class T>
typename std::enable_if<
    !std::is_same<T, Kokkos::complex<float>>::value && !std::is_same<T, Kokkos::complex<double>>::value, void>::type
Syr2Tester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::compareVanillaAgainstExpected(
    const T& alpha, const _ViewTypeExpected& h_vanilla, const _ViewTypeExpected& h_expected) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  if (_N <= 2) {
    for (int i(0); i < _M; ++i) {
      for (int j(0); j < _N; ++j) {
        std::cout << "h_exp(" << i << "," << j << ") = " << h_expected(i, j) << ", h_van(" << i << "," << j
                  << ") = " << h_vanilla(i, j) << std::endl;
      }
    }
  }
#endif
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
          << ", _useUpOption = " << _useUpOption << ": vanilla differs too much from expected"
          << ", numErrorsAbs = " << numErrorsAbs << ", numErrorsRel = " << numErrorsRel
          << ", maxErrorRel = " << maxErrorRel << ", iForMaxErrorRel = " << iForMaxErrorRel
          << ", jForMaxErrorRel = " << jForMaxErrorRel << ", h_expected(i,j) = "
          << (((_M > 0) && (_N > 0)) ? h_expected(iForMaxErrorRel, jForMaxErrorRel) : 9.999e+99)
          << ", h_vanilla(i,j) = " << (((_M > 0) && (_N > 0)) ? h_vanilla(iForMaxErrorRel, jForMaxErrorRel) : 9.999e+99)
          << ", maxNumErrorsAllowed = " << maxNumErrorsAllowed;

      int numErrors(numErrorsAbs + numErrorsRel);
#ifdef HAVE_KOKKOSKERNELS_DEBUG
      if (numErrors > 0) {
        std::cout << "WARNING" << msg.str() << std::endl;
      }
#endif
      EXPECT_LE(numErrors, maxNumErrorsAllowed) << "Failed test" << msg.str();
    }
  } else {
    int numErrors(0);

    for (int i(0); i < _M; ++i) {
      for (int j(0); j < _N; ++j) {
        if (h_expected(i, j) != h_vanilla(i, j)) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
          if (numErrors == 0) {
            std::cout << "ERROR, i = " << i << ", j = " << j << ": h_expected(i,j) = " << h_expected(i, j)
                      << ", h_vanilla(i,j) = " << h_vanilla(i, j) << std::endl;
          }
#endif
          numErrors++;
        }
      }  // for j
    }    // for i
    EXPECT_EQ(numErrors, 0) << "Failed test"
                            << ", A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr
                            << ", _A_is_ll = " << _A_is_ll << ", alpha type = " << typeid(alpha).name()
                            << ", _useHermitianOption = " << _useHermitianOption << ", _useUpOption = " << _useUpOption
                            << ": vanilla result is incorrect"
                            << ", numErrors = " << numErrors;
  }
}

// Code for complex values
template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
template <class T>
typename std::enable_if<
    std::is_same<T, Kokkos::complex<float>>::value || std::is_same<T, Kokkos::complex<double>>::value, void>::type
Syr2Tester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::compareKkSyr2AgainstReference(
    const T& alpha, const _HostViewTypeA& h_A, const _ViewTypeExpected& h_reference) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  if (_N <= 2) {
    for (int i(0); i < _M; ++i) {
      for (int j(0); j < _N; ++j) {
        std::cout << "h_exp(" << i << "," << j << ") = " << h_reference(i, j) << ", h_A(" << i << "," << j
                  << ") = " << h_A(i, j) << std::endl;
      }
    }
  }
#endif
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
      diff          = _KAT_A::abs(h_reference(i, j).real() - h_A(i, j).real());
      errorHappened = false;
      if (h_reference(i, j).real() == 0.) {
        diffThreshold = _KAT_A::abs(_absTol);
        if (diff > diffThreshold) {
          errorHappened = true;
          numErrorsRealAbs++;
        }
      } else {
        _AuxType aux = diff / _KAT_A::abs(h_reference(i, j).real());
        if (maxErrorRealRel < aux) {
          maxErrorRealRel     = aux;
          iForMaxErrorRealRel = i;
          jForMaxErrorRealRel = j;
        }

        diffThreshold = _KAT_A::abs(_relTol * h_reference(i, j).real());
        if (diff > diffThreshold) {
          errorHappened = true;
          numErrorsRealRel++;
        }
      }
      if (errorHappened && (numErrorsRealAbs + numErrorsRealRel == 1)) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
        std::cout << "ERROR, i = " << i << ", j = " << j << ": h_reference(i,j).real() = " << h_reference(i, j).real()
                  << ", h_A(i,j).real() = " << h_A(i, j).real()
                  << ", _KAT_A::abs(h_reference(i,j).real() - h_A(i,j).real()) = " << diff
                  << ", diffThreshold = " << diffThreshold << std::endl;
#endif
      }
      diff          = _KAT_A::abs(h_reference(i, j).imag() - h_A(i, j).imag());
      errorHappened = false;
      if (h_reference(i, j).imag() == 0.) {
        diffThreshold = _KAT_A::abs(_absTol);
        if (diff > diffThreshold) {
          errorHappened = true;
          numErrorsImagAbs++;
        }
      } else {
        _AuxType aux = diff / _KAT_A::abs(h_reference(i, j).imag());
        if (maxErrorImagRel < aux) {
          maxErrorImagRel     = aux;
          iForMaxErrorImagRel = i;
          jForMaxErrorImagRel = j;
        }

        diffThreshold = _KAT_A::abs(_relTol * h_reference(i, j).imag());
        if (diff > diffThreshold) {
          errorHappened = true;
          numErrorsImagRel++;
        }
      }
      if (errorHappened && (numErrorsImagAbs + numErrorsImagRel == 1)) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
        std::cout << "ERROR, i = " << i << ", j = " << j << ": h_reference(i,j).imag() = " << h_reference(i, j).imag()
                  << ", h_A(i,j).imag() = " << h_A(i, j).imag()
                  << ", _KAT_A::abs(h_reference(i,j).imag() - h_A(i,j).imag()) = " << diff
                  << ", diffThreshold = " << diffThreshold << std::endl;
#endif
      }
    }  // for j
  }    // for i

#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr << ", _A_is_ll = " << _A_is_ll
            << ", alpha type = " << typeid(alpha).name() << ", _useHermitianOption = " << _useHermitianOption
            << ", _useUpOption = " << _useUpOption << ", numErrorsRealAbs = " << numErrorsRealAbs
            << ", numErrorsRealRel = " << numErrorsRealRel << ", maxErrorRealRel = " << maxErrorRealRel
            << ", iForMaxErrorRealRel = " << iForMaxErrorRealRel << ", jForMaxErrorRealRel = " << jForMaxErrorRealRel
            << ", h_reference(i,j).real() = "
            << (((_M > 0) && (_N > 0)) ? h_reference(iForMaxErrorRealRel, jForMaxErrorRealRel).real() : 9.999e+99)
            << ", h_A(i,j).real() = "
            << (((_M > 0) && (_N > 0)) ? h_A(iForMaxErrorRealRel, jForMaxErrorRealRel).real() : 9.999e+99)
            << ", numErrorsImagAbs = " << numErrorsImagAbs << ", numErrorsImagRel = " << numErrorsImagRel
            << ", maxErrorImagRel = " << maxErrorImagRel << ", iForMaxErrorImagRel = " << iForMaxErrorImagRel
            << ", jForMaxErrorImagRel = " << jForMaxErrorImagRel << ", h_reference(i,j).imag() = "
            << (((_M > 0) && (_N > 0)) ? h_reference(iForMaxErrorImagRel, jForMaxErrorImagRel).imag() : 9.999e+99)
            << ", h_A(i,j).imag() = "
            << (((_M > 0) && (_N > 0)) ? h_A(iForMaxErrorImagRel, jForMaxErrorImagRel).imag() : 9.999e+99)
            << ", maxNumErrorsAllowed = " << maxNumErrorsAllowed << std::endl;
  if ((_M == 2131) && (_N == 2131)) {
    std::cout << "Information"
              << ": A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr << ", _A_is_ll = " << _A_is_ll
              << ", alpha type = " << typeid(alpha).name() << ", _useHermitianOption = " << _useHermitianOption
              << ", _useUpOption = " << _useUpOption << ", h_reference(11, 2119) = (" << h_reference(11, 2119).real()
              << ", " << h_reference(11, 2119).imag() << ")"
              << ", h_A(11, 2119) = (" << h_A(11, 2119).real() << ", " << h_A(11, 2119).imag() << ")" << std::endl;
    std::cout << "Information"
              << ": A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr << ", _A_is_ll = " << _A_is_ll
              << ", alpha type = " << typeid(alpha).name() << ", _useHermitianOption = " << _useHermitianOption
              << ", _useUpOption = " << _useUpOption << ", h_reference(710, 1065) = (" << h_reference(710, 1065).real()
              << ", " << h_reference(710, 1065).imag() << ")"
              << ", h_A(710, 1065) = (" << h_A(710, 1065).real() << ", " << h_A(710, 1065).imag() << ")" << std::endl;
  }
#endif
  {
    std::ostringstream msg;
    msg << ", A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr << ", _A_is_ll = " << _A_is_ll
        << ", alpha type = " << typeid(alpha).name() << ", _useHermitianOption = " << _useHermitianOption
        << ", _useUpOption = " << _useUpOption << ": syr2 result is incorrect on real components"
        << ", numErrorsRealAbs = " << numErrorsRealAbs << ", numErrorsRealRel = " << numErrorsRealRel
        << ", maxErrorRealRel = " << maxErrorRealRel << ", iForMaxErrorRealRel = " << iForMaxErrorRealRel
        << ", jForMaxErrorRealRel = " << jForMaxErrorRealRel << ", h_reference(i,j).real() = "
        << (((_M > 0) && (_N > 0)) ? h_reference(iForMaxErrorRealRel, jForMaxErrorRealRel).real() : 9.999e+99)
        << ", h_A(i,j).real() = "
        << (((_M > 0) && (_N > 0)) ? h_A(iForMaxErrorRealRel, jForMaxErrorRealRel).real() : 9.999e+99)
        << ", maxNumErrorsAllowed = " << maxNumErrorsAllowed;

    int numErrorsReal(numErrorsRealAbs + numErrorsRealRel);
#ifdef HAVE_KOKKOSKERNELS_DEBUG
    if (numErrorsReal > 0) {
      std::cout << "WARNING" << msg.str() << std::endl;
    }
#endif
    EXPECT_LE(numErrorsReal, maxNumErrorsAllowed) << "Failed test" << msg.str();
  }
  {
    std::ostringstream msg;
    msg << ", A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr << ", _A_is_ll = " << _A_is_ll
        << ", alpha type = " << typeid(alpha).name() << ", _useHermitianOption = " << _useHermitianOption
        << ", _useUpOption = " << _useUpOption << ": syr2 result is incorrect on imag components"
        << ", numErrorsImagAbs = " << numErrorsImagAbs << ", numErrorsImagRel = " << numErrorsImagRel
        << ", maxErrorImagRel = " << maxErrorImagRel << ", iForMaxErrorImagRel = " << iForMaxErrorImagRel
        << ", jForMaxErrorImagRel = " << jForMaxErrorImagRel << ", h_reference(i,j).imag() = "
        << (((_M > 0) && (_N > 0)) ? h_reference(iForMaxErrorImagRel, jForMaxErrorImagRel).imag() : 9.999e+99)
        << ", h_A(i,j).imag() = "
        << (((_M > 0) && (_N > 0)) ? h_A(iForMaxErrorImagRel, jForMaxErrorImagRel).imag() : 9.999e+99)
        << ", maxNumErrorsAllowed = " << maxNumErrorsAllowed;

    int numErrorsImag(numErrorsImagAbs + numErrorsImagRel);
#ifdef HAVE_KOKKOSKERNELS_DEBUG
    if (numErrorsImag > 0) {
      std::cout << "WARNING" << msg.str() << std::endl;
    }
#endif
    EXPECT_LE(numErrorsImag, maxNumErrorsAllowed) << "Failed test" << msg.str();
  }
}

// Code for non-complex values
template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
template <class T>
typename std::enable_if<
    !std::is_same<T, Kokkos::complex<float>>::value && !std::is_same<T, Kokkos::complex<double>>::value, void>::type
Syr2Tester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::compareKkSyr2AgainstReference(
    const T& alpha, const _HostViewTypeA& h_A, const _ViewTypeExpected& h_reference) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  if (_N <= 2) {
    for (int i(0); i < _M; ++i) {
      for (int j(0); j < _N; ++j) {
        std::cout << "h_exp(" << i << "," << j << ") = " << h_reference(i, j) << ", h_A(" << i << "," << j
                  << ") = " << h_A(i, j) << std::endl;
      }
    }
  }
#endif
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
      diff          = _KAT_A::abs(h_reference(i, j) - h_A(i, j));
      errorHappened = false;
      if (h_reference(i, j) == 0.) {
        diffThreshold = _KAT_A::abs(_absTol);
        if (diff > diffThreshold) {
          errorHappened = true;
          numErrorsAbs++;
        }
      } else {
        _AuxType aux = diff / _KAT_A::abs(h_reference(i, j));
        if (maxErrorRel < aux) {
          maxErrorRel     = aux;
          iForMaxErrorRel = i;
          jForMaxErrorRel = j;
        }

        diffThreshold = _KAT_A::abs(_relTol * h_reference(i, j));
        if (diff > diffThreshold) {
          errorHappened = true;
          numErrorsRel++;
        }
      }
      if (errorHappened && (numErrorsAbs + numErrorsRel == 1)) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
        std::cout << "ERROR, i = " << i << ", j = " << j << ": h_reference(i,j) = " << h_reference(i, j)
                  << ", h_A(i,j) = " << h_A(i, j) << ", _KAT_A::abs(h_reference(i,j) - h_A(i,j)) = " << diff
                  << ", diffThreshold = " << diffThreshold << std::endl;
#endif
      }
    }  // for j
  }    // for i
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr << ", _A_is_ll = " << _A_is_ll
            << ", alpha type = " << typeid(alpha).name() << ", _useHermitianOption = " << _useHermitianOption
            << ", _useUpOption = " << _useUpOption << ", numErrorsAbs = " << numErrorsAbs
            << ", numErrorsRel = " << numErrorsRel << ", maxErrorRel = " << maxErrorRel
            << ", iForMaxErrorRel = " << iForMaxErrorRel << ", jForMaxErrorRel = " << jForMaxErrorRel
            << ", h_reference(i,j) = "
            << (((_M > 0) && (_N > 0)) ? h_reference(iForMaxErrorRel, jForMaxErrorRel) : 9.999e+99)
            << ", h_A(i,j) = " << (((_M > 0) && (_N > 0)) ? h_A(iForMaxErrorRel, jForMaxErrorRel) : 9.999e+99)
            << ", maxNumErrorsAllowed = " << maxNumErrorsAllowed << std::endl;
#endif
  {
    std::ostringstream msg;
    msg << ", A is " << _M << " by " << _N << ", _A_is_lr = " << _A_is_lr << ", _A_is_ll = " << _A_is_ll
        << ", alpha type = " << typeid(alpha).name() << ", _useHermitianOption = " << _useHermitianOption
        << ", _useUpOption = " << _useUpOption << ": syr2 result is incorrect"
        << ", numErrorsAbs = " << numErrorsAbs << ", numErrorsRel = " << numErrorsRel
        << ", maxErrorRel = " << maxErrorRel << ", iForMaxErrorRel = " << iForMaxErrorRel
        << ", jForMaxErrorRel = " << jForMaxErrorRel << ", h_reference(i,j) = "
        << (((_M > 0) && (_N > 0)) ? h_reference(iForMaxErrorRel, jForMaxErrorRel) : 9.999e+99)
        << ", h_A(i,j) = " << (((_M > 0) && (_N > 0)) ? h_A(iForMaxErrorRel, jForMaxErrorRel) : 9.999e+99)
        << ", maxNumErrorsAllowed = " << maxNumErrorsAllowed;

    int numErrors(numErrorsAbs + numErrorsRel);
#ifdef HAVE_KOKKOSKERNELS_DEBUG
    if (numErrors > 0) {
      std::cout << "WARNING" << msg.str() << std::endl;
    }
#endif
    EXPECT_LE(numErrors, maxNumErrorsAllowed) << "Failed test" << msg.str();
  }
}

template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
template <class TX, class TY>
void Syr2Tester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::callKkSyr2AndCompareAgainstExpected(
    const ScalarA& alpha, TX& x, TY& y, view_stride_adapter<_ViewTypeA, false>& A, const _ViewTypeExpected& h_expected,
    const std::string& situation) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "In Test_Blas2_syr2, '" << situation << "', alpha = " << alpha << std::endl;
  std::cout << "In Test_Blas2_syr2.hpp, right before calling KokkosBlas::syr2()"
            << ": ViewTypeA = " << typeid(_ViewTypeA).name()
            << ", _kkSyr2ShouldThrowException = " << _kkSyr2ShouldThrowException << std::endl;
#endif
  std::string mode = _useHermitianOption ? "H" : "T";
  std::string uplo = _useUpOption ? "U" : "L";
  bool gotStdException(false);
  bool gotUnknownException(false);
  try {
    KokkosBlas::syr2(mode.c_str(), uplo.c_str(), alpha, x, y, A.d_view);
    Kokkos::fence();
  } catch (const std::exception& e) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
    std::cout << "In Test_Blas2_syr2, '" << situation << "': caught exception, e.what() = " << e.what() << std::endl;
#endif
    gotStdException = true;
  } catch (...) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
    std::cout << "In Test_Blas2_syr2, '" << situation << "': caught unknown exception" << std::endl;
#endif
    gotUnknownException = true;
  }

  EXPECT_EQ(gotUnknownException, false) << "Failed test, '" << situation
                                        << "': unknown exception should not have happened";

  EXPECT_EQ(gotStdException, _kkSyr2ShouldThrowException)
      << "Failed test, '" << situation << "': kk syr2() should" << (_kkSyr2ShouldThrowException ? " " : " not ")
      << "have thrown a std::exception";

  if ((gotStdException == false) && (gotUnknownException == false)) {
    Kokkos::deep_copy(A.h_base, A.d_base);
    this->compareKkSyr2AgainstReference(alpha, A.h_view, h_expected);
  }
}

template <class ScalarX, class tLayoutX, class ScalarY, class tLayoutY, class ScalarA, class tLayoutA, class Device>
template <class TX, class TY>
void Syr2Tester<ScalarX, tLayoutX, ScalarY, tLayoutY, ScalarA, tLayoutA, Device>::callKkGerAndCompareKkSyr2AgainstIt(
    const ScalarA& alpha, TX& x, TY& y, view_stride_adapter<_ViewTypeA, false>& org_A, const _HostViewTypeA& h_A_syr2,
    const std::string& situation) {
  view_stride_adapter<_ViewTypeA, false> A_ger("A_ger", _M, _N);
  Kokkos::deep_copy(A_ger.d_base, org_A.d_base);

  // ********************************************************************
  // Call ger()
  // ********************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "In Test_Blas2_syr2, '" << situation << "', alpha = " << alpha << std::endl;
  std::cout << "In Test_Blas2_syr2.hpp, right before calling KokkosBlas::ger()"
            << ": ViewTypeA = " << typeid(_ViewTypeA).name()
            << ", _kkGerShouldThrowException = " << _kkGerShouldThrowException << std::endl;
#endif
  std::string mode = _useHermitianOption ? "H" : "T";
  bool gotStdException(false);
  bool gotUnknownException(false);
  try {
    KokkosBlas::ger(mode.c_str(), alpha, x, y, A_ger.d_view);
    Kokkos::fence();
  } catch (const std::exception& e) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
    std::cout << "In Test_Blas2_syr2, '" << situation << "', ger() call 1: caught exception, e.what() = " << e.what()
              << std::endl;
#endif
    gotStdException = true;
  } catch (...) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
    std::cout << "In Test_Blas2_syr2, '" << situation << "', ger() call 1: caught unknown exception" << std::endl;
#endif
    gotUnknownException = true;
  }

  EXPECT_EQ(gotUnknownException, false) << "Failed test, '" << situation
                                        << "': unknown exception should not have happened for ger() call 1";

  EXPECT_EQ(gotStdException, false) << "Failed test, '" << situation
                                    << "': kk ger() 1 should not have thrown a std::exception";

  // ********************************************************************
  // Call ger() again
  // ********************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "In Test_Blas2_syr2.hpp, right before calling KokkosBlas::ger() again";
#endif
  try {
    if (_useHermitianOption) {
      KokkosBlas::ger(mode.c_str(), _KAT_A::conj(alpha), y, x, A_ger.d_view);
    } else {
      KokkosBlas::ger(mode.c_str(), alpha, y, x, A_ger.d_view);
    }
    Kokkos::fence();
  } catch (const std::exception& e) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
    std::cout << "In Test_Blas2_syr2, '" << situation << "', ger() call 2: caught exception, e.what() = " << e.what()
              << std::endl;
#endif
    gotStdException = true;
  } catch (...) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
    std::cout << "In Test_Blas2_syr2, '" << situation << "', ger() call 2: caught unknown exception" << std::endl;
#endif
    gotUnknownException = true;
  }

  EXPECT_EQ(gotUnknownException, false) << "Failed test, '" << situation
                                        << "': unknown exception should not have happened for ger() call 2";

  EXPECT_EQ(gotStdException, false) << "Failed test, '" << situation
                                    << "': kk ger() 2 should not have thrown a std::exception";

  // ********************************************************************
  // Prepare h_ger_reference to be compared against h_A_syr2
  // ********************************************************************
  view_stride_adapter<_ViewTypeExpected, true> h_ger_reference("h_ger_reference", _M, _N);
  Kokkos::deep_copy(h_ger_reference.d_base, A_ger.d_base);
  Kokkos::deep_copy(h_ger_reference.h_base, h_ger_reference.d_base);

  std::string uplo = _useUpOption ? "U" : "L";
  for (int i = 0; i < _M; ++i) {
    for (int j = 0; j < _N; ++j) {
      if (((_useUpOption == true) && (i <= j)) || ((_useUpOption == false) && (i >= j))) {
        // Keep h_ger_reference as already computed
      } else {
        h_ger_reference.h_view(i, j) = org_A.h_view(i, j);
      }
    }
  }
  if (_useHermitianOption && _A_is_complex) {
    for (int i(0); i < _N; ++i) {
      h_ger_reference.h_view(i, i) = 0.5 * (h_ger_reference.h_view(i, i) + _KAT_A::conj(h_ger_reference.h_view(i, i)));
    }
  }

  // ********************************************************************
  // Compare
  // ********************************************************************
  this->compareKkSyr2AgainstReference(alpha, h_A_syr2, h_ger_reference.h_view);
}

}  // namespace Test

template <class ScalarX, class ScalarY, class ScalarA, class Device>
#ifdef HAVE_KOKKOSKERNELS_DEBUG
int test_syr2(const std::string& caseName) {
  std::cout << "+=============================================================="
               "============"
            << std::endl;
  std::cout << "Starting " << caseName << "..." << std::endl;
#else
int test_syr2(const std::string& /*caseName*/) {
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
  std::cout << "+--------------------------------------------------------------"
               "------------"
            << std::endl;
  std::cout << "Starting " << caseName << " for LAYOUTLEFT ..." << std::endl;
#endif
  if (true) {
    Test::Syr2Tester<ScalarX, Kokkos::LayoutLeft, ScalarY, Kokkos::LayoutLeft, ScalarA, Kokkos::LayoutLeft, Device>
        tester;
    tester.test(0, 0);
    tester.test(1, 0);
    tester.test(2, 0);
    tester.test(13, 0);
    tester.test(1024, 0);

    if (useAnalyticalResults) {
      tester.test(1024, 0, true, false, false);
      tester.test(1024, 0, true, false, true);
      tester.test(1024, 0, true, true, false);
      tester.test(1024, 0, true, true, true);
    }

    tester.test(2, 0, false, false, true);
    tester.test(50, 0, false, false, true);
    tester.test(2, 0, false, true, false);
    tester.test(50, 0, false, true, false);
    tester.test(2, 0, false, true, true);
    tester.test(50, 0, false, true, true);

    tester.test(50, 4);
    tester.test(2131, 0);
  }

#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Finished " << caseName << " for LAYOUTLEFT" << std::endl;
  std::cout << "+--------------------------------------------------------------"
               "------------"
            << std::endl;
#endif
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "+--------------------------------------------------------------"
               "------------"
            << std::endl;
  std::cout << "Starting " << caseName << " for LAYOUTRIGHT ..." << std::endl;
#endif
  if (true) {
    Test::Syr2Tester<ScalarX, Kokkos::LayoutRight, ScalarY, Kokkos::LayoutRight, ScalarA, Kokkos::LayoutRight, Device>
        tester;
    tester.test(0, 0);
    tester.test(1, 0);
    tester.test(2, 0);
    tester.test(13, 0);
    tester.test(1024, 0);

    if (useAnalyticalResults) {
      tester.test(1024, 0, true, false, false);
      tester.test(1024, 0, true, false, true);
      tester.test(1024, 0, true, true, false);
      tester.test(1024, 0, true, true, true);
    }

    tester.test(2, 0, false, false, true);
    tester.test(50, 0, false, false, true);
    tester.test(2, 0, false, true, false);
    tester.test(50, 0, false, true, false);
    tester.test(2, 0, false, true, true);
    tester.test(50, 0, false, true, true);

    tester.test(50, 4);
    tester.test(2131, 0);
  }

#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Finished " << caseName << " for LAYOUTRIGHT" << std::endl;
  std::cout << "+--------------------------------------------------------------"
               "------------"
            << std::endl;
#endif
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "+--------------------------------------------------------------"
               "------------"
            << std::endl;
  std::cout << "Starting " << caseName << " for LAYOUTSTRIDE ..." << std::endl;
#endif
  if (true) {
    Test::Syr2Tester<ScalarX, Kokkos::LayoutStride, ScalarY, Kokkos::LayoutStride, ScalarA, Kokkos::LayoutStride,
                     Device>
        tester;
    tester.test(0, 0);
    tester.test(1, 0);
    tester.test(2, 0);
    tester.test(13, 0);
    tester.test(1024, 0);

    if (useAnalyticalResults) {
      tester.test(1024, 0, true, false, false);
      tester.test(1024, 0, true, false, true);
      tester.test(1024, 0, true, true, false);
      tester.test(1024, 0, true, true, true);
    }

    tester.test(2, 0, false, false, true);
    tester.test(50, 0, false, false, true);
    tester.test(2, 0, false, true, false);
    tester.test(50, 0, false, true, false);
    tester.test(2, 0, false, true, true);
    tester.test(50, 0, false, true, true);

    tester.test(50, 4);
    tester.test(2131, 0);
  }

#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Finished " << caseName << " for LAYOUTSTRIDE" << std::endl;
  std::cout << "+--------------------------------------------------------------"
               "------------"
            << std::endl;
#endif
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "+--------------------------------------------------------------"
               "------------"
            << std::endl;
  std::cout << "Starting " << caseName << " for MIXED LAYOUTS ..." << std::endl;
#endif
  if (true) {
    Test::Syr2Tester<ScalarX, Kokkos::LayoutStride, ScalarY, Kokkos::LayoutLeft, ScalarA, Kokkos::LayoutRight, Device>
        tester;
    tester.test(1, 0);
    tester.test(2, 0);
    tester.test(1024, 0);

    if (useAnalyticalResults) {
      tester.test(1024, 0, true, false, true);
      tester.test(1024, 0, true, true, true);
    }

    tester.test(2, 0, false, false, true);
    tester.test(50, 0, false, false, true);
    tester.test(2, 0, false, true, true);
    tester.test(50, 0, false, true, true);
  }

  if (true) {
    Test::Syr2Tester<ScalarX, Kokkos::LayoutLeft, ScalarX, Kokkos::LayoutStride, ScalarA, Kokkos::LayoutRight, Device>
        tester;
    tester.test(1024, 0);
  }

#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Finished " << caseName << " for MIXED LAYOUTS" << std::endl;
  std::cout << "+--------------------------------------------------------------"
               "------------"
            << std::endl;
#endif
#endif

#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Finished " << caseName << std::endl;
  std::cout << "+=============================================================="
               "============"
            << std::endl;
#endif
  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, syr2_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::syr2_float");
  test_syr2<float, float, float, TestDevice>("test case syr2_float");
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, syr2_complex_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::syr2_complex_float");
  test_syr2<Kokkos::complex<float>, Kokkos::complex<float>, Kokkos::complex<float>, TestDevice>(
      "test case syr2_complex_float");
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, syr2_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::syr2_double");
  test_syr2<double, double, double, TestDevice>("test case syr2_double");
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, syr2_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::syr2_complex_double");
  test_syr2<Kokkos::complex<double>, Kokkos::complex<double>, Kokkos::complex<double>, TestDevice>(
      "test case syr2_complex_double");
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, syr2_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::syr2_int");
  test_syr2<int, int, int, TestDevice>("test case syr2_int");
  Kokkos::Profiling::popRegion();
}
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
TEST_F(TestCategory, syr2_int_float_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::syr2_int_float_double");
  test_syr2<int, float, double, TestDevice>("test case syr2_mixed_types");
  Kokkos::Profiling::popRegion();
}
#endif
