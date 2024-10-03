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
// the operation A += alpha * x * x^{T,H}:
// 01) Type of 'x' components: float, double, complex, ...
// 02) Type of 'A' components: float, double, complex, ...
// 03) Execution space: serial, threads, OpenMP, Cuda, ...
// 04) Layout of 'x'
// 05) Layout of 'A'
// 06) Dimension of 'A'
// 07) Options 'const' or 'non const' for x view, when calling syr()
// 08) Usage of analytical results in the tests
// 09) Options 'T' or 'H' when calling syr()
// 10) Options 'U' or 'L' when calling syr()
//
// Choices (01)-(03) are selected in the routines TEST_F() at the
// very bottom of the file, when calling test_syr<...>().
//
// Choices (04)-(10) are selected in routine test_syr<...>(),
// when calling the method test() of class Test::SyrTester<...>.
//
// The class Test::SyrTester<...> represents the "core" of the test
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
#include <KokkosBlas2_syr.hpp>
#include <Kokkos_MathematicalConstants.hpp>

namespace Test {

template <class ScalarX, class tLayoutX, class ScalarA, class tLayoutA, class Device>
class SyrTester {
 public:
  SyrTester();

  ~SyrTester();

  void test(const int N, const int nonConstConstCombinations, const bool useAnalyticalResults = false,
            const bool useHermitianOption = false, const bool useUpOption = false);

 private:
  using _ViewTypeX = Kokkos::View<ScalarX*, tLayoutX, Device>;
  using _ViewTypeA = Kokkos::View<ScalarA**, tLayoutA, Device>;

  using _HostViewTypeX    = typename _ViewTypeX::HostMirror;
  using _HostViewTypeA    = typename _ViewTypeA::HostMirror;
  using _ViewTypeExpected = Kokkos::View<ScalarA**, tLayoutA, Kokkos::HostSpace>;

  using _KAT_A   = Kokkos::ArithTraits<ScalarA>;
  using _AuxType = typename _KAT_A::mag_type;

  void populateVariables(ScalarA& alpha, view_stride_adapter<_ViewTypeX, false>& x,
                         view_stride_adapter<_ViewTypeA, false>& A, _ViewTypeExpected& h_expected,
                         bool& expectedResultIsKnown);

  template <class T>
  typename std::enable_if<
      std::is_same<T, Kokkos::complex<float>>::value || std::is_same<T, Kokkos::complex<double>>::value, void>::type
  populateAnalyticalValues(T& alpha, _HostViewTypeX& h_x, _HostViewTypeA& h_A, _ViewTypeExpected& h_expected);

  template <class T>
  typename std::enable_if<
      !std::is_same<T, Kokkos::complex<float>>::value && !std::is_same<T, Kokkos::complex<double>>::value, void>::type
  populateAnalyticalValues(T& alpha, _HostViewTypeX& h_x, _HostViewTypeA& h_A, _ViewTypeExpected& h_expected);

  template <class T>
  typename std::enable_if<
      std::is_same<T, Kokkos::complex<float>>::value || std::is_same<T, Kokkos::complex<double>>::value, void>::type
  populateVanillaValues(const T& alpha, const _HostViewTypeX& h_x, const _HostViewTypeA& h_A,
                        _ViewTypeExpected& h_vanilla);

  template <class T>
  typename std::enable_if<
      !std::is_same<T, Kokkos::complex<float>>::value && !std::is_same<T, Kokkos::complex<double>>::value, void>::type
  populateVanillaValues(const T& alpha, const _HostViewTypeX& h_x, const _HostViewTypeA& h_A,
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
  compareKkSyrAgainstReference(const T& alpha, const _HostViewTypeA& h_A, const _ViewTypeExpected& h_reference);

  template <class T>
  typename std::enable_if<
      !std::is_same<T, Kokkos::complex<float>>::value && !std::is_same<T, Kokkos::complex<double>>::value, void>::type
  compareKkSyrAgainstReference(const T& alpha, const _HostViewTypeA& h_A, const _ViewTypeExpected& h_reference);

  template <class T>
  T shrinkAngleToZeroTwoPiRange(const T input);

  template <class TX>
  void callKkSyrAndCompareAgainstExpected(const ScalarA& alpha, TX& x, view_stride_adapter<_ViewTypeA, false>& A,
                                          const _ViewTypeExpected& h_expected, const std::string& situation);

  template <class TX>
  void callKkGerAndCompareKkSyrAgainstIt(const ScalarA& alpha, TX& x, view_stride_adapter<_ViewTypeA, false>& org_A,
                                         const _HostViewTypeA& h_A_syr, const std::string& situation);

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
  bool _kkSyrShouldThrowException;
  bool _kkGerShouldThrowException;
};

template <class ScalarX, class tLayoutX, class ScalarA, class tLayoutA, class Device>
SyrTester<ScalarX, tLayoutX, ScalarA, tLayoutA, Device>::SyrTester()
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
      _kkSyrShouldThrowException(false),
      _kkGerShouldThrowException(false) {
}

template <class ScalarX, class tLayoutX, class ScalarA, class tLayoutA, class Device>
SyrTester<ScalarX, tLayoutX, ScalarA, tLayoutA, Device>::~SyrTester() {
  // Nothing to do
}

template <class ScalarX, class tLayoutX, class ScalarA, class tLayoutA, class Device>
void SyrTester<ScalarX, tLayoutX, ScalarA, tLayoutA, Device>::test(const int N, const int nonConstConstCombinations,
                                                                   const bool useAnalyticalResults,
                                                                   const bool useHermitianOption,
                                                                   const bool useUpOption) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Entering SyrTester::test()... - - - - - - - - - - - - - - - - "
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
  _kkSyrShouldThrowException = false;

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
  view_stride_adapter<_ViewTypeA, false> A("A", _M, _N);

  view_stride_adapter<_ViewTypeExpected, true> h_expected("expected A += alpha * x * x^{t,h}", _M, _N);
  bool expectedResultIsKnown = false;

  ScalarA alpha(_KAT_A::zero());

  // ********************************************************************
  // Step 2 of 7: populate alpha, h_x, h_A, h_expected, x, A
  // ********************************************************************
  this->populateVariables(alpha, x, A, h_expected.d_view, expectedResultIsKnown);

  // ********************************************************************
  // Step 3 of 7: populate h_vanilla
  // ********************************************************************
  view_stride_adapter<_ViewTypeExpected, true> h_vanilla("vanilla = A + alpha * x * x^{t,h}", _M, _N);
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  Kokkos::printf("In Test_Blas2_syr.hpp, computing vanilla A with alpha type = %s\n", typeid(alpha).name());
#endif
  this->populateVanillaValues(alpha, x.h_view, A.h_view, h_vanilla.d_view);

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
    this->callKkSyrAndCompareAgainstExpected(alpha, x.d_view, A, h_expected.d_view, "non const x");

    if ((_useAnalyticalResults == false) &&  // Just to save run time
        (_kkGerShouldThrowException == false)) {
      this->callKkGerAndCompareKkSyrAgainstIt(alpha, x.d_view, org_A, A.h_view, "non const x");
    }
  }

  // ********************************************************************
  // Step 6 of 7: test with const x
  // ********************************************************************
  if (test_cx) {
    Kokkos::deep_copy(A.d_base, org_A.d_base);

    this->callKkSyrAndCompareAgainstExpected(alpha, x.d_view_const, A, h_expected.d_view, "const x");
  }

  // ********************************************************************
  // Step 7 of 7: tests with invalid values on the first input parameter
  // ********************************************************************
  EXPECT_ANY_THROW(KokkosBlas::syr(".", "U", alpha, x.d_view, A.d_view))
      << "Failed test: kk syr should have thrown an exception for mode '.'";
  EXPECT_ANY_THROW(KokkosBlas::syr("", "U", alpha, x.d_view, A.d_view))
      << "Failed test: kk syr should have thrown an exception for mode ''";
  EXPECT_ANY_THROW(KokkosBlas::syr("T", ".", alpha, x.d_view, A.d_view))
      << "Failed test: kk syr should have thrown an exception for uplo '.'";
  EXPECT_ANY_THROW(KokkosBlas::syr("T", "", alpha, x.d_view, A.d_view))
      << "Failed test: kk syr should have thrown an exception for uplo ''";

#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Leaving SyrTester::test() - - - - - - - - - - - - - - - - - - "
               "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
               "- - - - - - - "
            << std::endl;
#endif
}

template <class ScalarX, class tLayoutX, class ScalarA, class tLayoutA, class Device>
void SyrTester<ScalarX, tLayoutX, ScalarA, tLayoutA, Device>::populateVariables(
    ScalarA& alpha, view_stride_adapter<_ViewTypeX, false>& x, view_stride_adapter<_ViewTypeA, false>& A,
    _ViewTypeExpected& h_expected, bool& expectedResultIsKnown) {
  expectedResultIsKnown = false;

  if (_useAnalyticalResults) {
    this->populateAnalyticalValues(alpha, x.h_view, A.h_view, h_expected);
    Kokkos::deep_copy(x.d_base, x.h_base);
    Kokkos::deep_copy(A.d_base, A.h_base);

    expectedResultIsKnown = true;
  } else if (_N == 1) {
    alpha = 3;

    x.h_view[0] = 2;

    A.h_view(0, 0) = 7;

    Kokkos::deep_copy(x.d_base, x.h_base);
    Kokkos::deep_copy(A.d_base, A.h_base);

    h_expected(0, 0)      = 19;
    expectedResultIsKnown = true;
  } else if (_N == 2) {
    alpha = 3;

    x.h_view[0] = -2;
    x.h_view[1] = 9;

    A.h_view(0, 0) = 17;
    A.h_view(0, 1) = -43;
    A.h_view(1, 0) = -43;
    A.h_view(1, 1) = 101;

    Kokkos::deep_copy(x.d_base, x.h_base);
    Kokkos::deep_copy(A.d_base, A.h_base);

    if (_useUpOption) {
      h_expected(0, 0) = 29;
      h_expected(0, 1) = -97;
      h_expected(1, 0) = -43;
      h_expected(1, 1) = 344;
    } else {
      h_expected(0, 0) = 29;
      h_expected(0, 1) = -43;
      h_expected(1, 0) = -97;
      h_expected(1, 1) = 344;
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
      ScalarA randStart, randEnd;
      Test::getRandomBounds(1.0, randStart, randEnd);
      Kokkos::fill_random(A.d_view, rand_pool, randStart, randEnd);
    }

    Kokkos::deep_copy(x.h_base, x.d_base);
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
        std::cout << "h_origA(" << i << "," << j << ")=" << A.h_view(i, j) << std::endl;
      }
    }
  }
#endif
}

// Code for complex values
template <class ScalarX, class tLayoutX, class ScalarA, class tLayoutA, class Device>
template <class T>
typename std::enable_if<
    std::is_same<T, Kokkos::complex<float>>::value || std::is_same<T, Kokkos::complex<double>>::value, void>::type
SyrTester<ScalarX, tLayoutX, ScalarA, tLayoutA, Device>::populateAnalyticalValues(T& alpha, _HostViewTypeX& h_x,
                                                                                  _HostViewTypeA& h_A,
                                                                                  _ViewTypeExpected& h_expected) {
  if (_useHermitianOption) {
    alpha.real() = 1.;
    alpha.imag() = 0.;
  } else {
    alpha.real() = 1.;
    alpha.imag() = -1.;
  }

  for (int i = 0; i < _M; ++i) {
    _AuxType auxI = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i));
    h_x[i].real() = sin(auxI);
    h_x[i].imag() = cos(auxI);
  }

  if (_useHermitianOption) {
    for (int i = 0; i < _M; ++i) {
      for (int j = 0; j < _N; ++j) {
        _AuxType auxImJ = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i - j));
        if (((_useUpOption == true) && (i <= j)) || ((_useUpOption == false) && (i >= j))) {
          h_A(i, j).real() = cos(auxImJ);
          h_A(i, j).imag() = -sin(auxImJ);
        } else {
          h_A(i, j).real() = cos(auxImJ);
          h_A(i, j).imag() = sin(auxImJ);
        }
      }
    }
  } else {
    for (int i = 0; i < _M; ++i) {
      for (int j = 0; j < _N; ++j) {
        _AuxType auxIpJ  = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i + j));
        h_A(i, j).real() = sin(auxIpJ) + cos(auxIpJ);
        h_A(i, j).imag() = sin(auxIpJ) - cos(auxIpJ);
      }
    }
  }

  if (_useHermitianOption) {
    for (int i = 0; i < _M; ++i) {
      for (int j = 0; j < _N; ++j) {
        if (((_useUpOption == true) && (i <= j)) || ((_useUpOption == false) && (i >= j))) {
          _AuxType auxImJ         = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i - j));
          h_expected(i, j).real() = 2. * cos(auxImJ);
          h_expected(i, j).imag() = -2. * sin(auxImJ);
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
          h_expected(i, j).real() = 2. * sin(auxIpJ);
          h_expected(i, j).imag() = 2. * sin(auxIpJ);
        } else {
          h_expected(i, j).real() = h_A(i, j).real();
          h_expected(i, j).imag() = h_A(i, j).imag();
        }
      }
    }
  }
}

// Code for non-complex values
template <class ScalarX, class tLayoutX, class ScalarA, class tLayoutA, class Device>
template <class T>
typename std::enable_if<
    !std::is_same<T, Kokkos::complex<float>>::value && !std::is_same<T, Kokkos::complex<double>>::value, void>::type
SyrTester<ScalarX, tLayoutX, ScalarA, tLayoutA, Device>::populateAnalyticalValues(T& alpha, _HostViewTypeX& h_x,
                                                                                  _HostViewTypeA& h_A,
                                                                                  _ViewTypeExpected& h_expected) {
  alpha = 2;

  for (int i = 0; i < _M; ++i) {
    _AuxType auxI = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i));
    h_x[i]        = sin(auxI);
  }

  for (int i = 0; i < _M; ++i) {
    _AuxType auxI = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i));
    for (int j = 0; j < _N; ++j) {
      _AuxType auxJ = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(j));
      h_A(i, j)     = 2 * cos(auxI) * cos(auxJ);
    }
  }

  for (int i = 0; i < _M; ++i) {
    for (int j = 0; j < _N; ++j) {
      if (((_useUpOption == true) && (i <= j)) || ((_useUpOption == false) && (i >= j))) {
        _AuxType auxImJ  = this->shrinkAngleToZeroTwoPiRange(static_cast<_AuxType>(i - j));
        h_expected(i, j) = 2 * cos(auxImJ);
      } else {
        h_expected(i, j) = h_A(i, j);
      }
    }
  }
}

// Code for complex values
template <class ScalarX, class tLayoutX, class ScalarA, class tLayoutA, class Device>
template <class T>
typename std::enable_if<
    std::is_same<T, Kokkos::complex<float>>::value || std::is_same<T, Kokkos::complex<double>>::value, void>::type
SyrTester<ScalarX, tLayoutX, ScalarA, tLayoutA, Device>::populateVanillaValues(const T& alpha,
                                                                               const _HostViewTypeX& h_x,
                                                                               const _HostViewTypeA& h_A,
                                                                               _ViewTypeExpected& h_vanilla) {
  if (_vanillaUsesDifferentOrderOfOps) {
    if (_useHermitianOption) {
      for (int i = 0; i < _M; ++i) {
        for (int j = 0; j < _N; ++j) {
          if (((_useUpOption == true) && (i <= j)) || ((_useUpOption == false) && (i >= j))) {
            h_vanilla(i, j) = h_A(i, j) + alpha * _KAT_A::conj(h_x(j)) * h_x(i);
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
            h_vanilla(i, j) = h_A(i, j) + alpha * h_x(j) * h_x(i);
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
            h_vanilla(i, j) = h_A(i, j) + alpha * h_x(i) * _KAT_A::conj(h_x(j));
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
            h_vanilla(i, j) = h_A(i, j) + alpha * h_x(i) * h_x(j);
          } else {
            h_vanilla(i, j) = h_A(i, j);
          }
        }
      }
    }
  }
}

// Code for non-complex values
template <class ScalarX, class tLayoutX, class ScalarA, class tLayoutA, class Device>
template <class T>
typename std::enable_if<
    !std::is_same<T, Kokkos::complex<float>>::value && !std::is_same<T, Kokkos::complex<double>>::value, void>::type
SyrTester<ScalarX, tLayoutX, ScalarA, tLayoutA, Device>::populateVanillaValues(const T& alpha,
                                                                               const _HostViewTypeX& h_x,
                                                                               const _HostViewTypeA& h_A,
                                                                               _ViewTypeExpected& h_vanilla) {
  if (_vanillaUsesDifferentOrderOfOps) {
    for (int i = 0; i < _M; ++i) {
      for (int j = 0; j < _N; ++j) {
        if (((_useUpOption == true) && (i <= j)) || ((_useUpOption == false) && (i >= j))) {
          h_vanilla(i, j) = h_A(i, j) + alpha * h_x(j) * h_x(i);
        } else {
          h_vanilla(i, j) = h_A(i, j);
        }
      }
    }
  } else {
    for (int i = 0; i < _M; ++i) {
      for (int j = 0; j < _N; ++j) {
        if (((_useUpOption == true) && (i <= j)) || ((_useUpOption == false) && (i >= j))) {
          h_vanilla(i, j) = h_A(i, j) + alpha * h_x(i) * h_x(j);
        } else {
          h_vanilla(i, j) = h_A(i, j);
        }
      }
    }
  }
}

template <class ScalarX, class tLayoutX, class ScalarA, class tLayoutA, class Device>
template <class T>
T SyrTester<ScalarX, tLayoutX, ScalarA, tLayoutA, Device>::shrinkAngleToZeroTwoPiRange(const T input) {
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
template <class ScalarX, class tLayoutX, class ScalarA, class tLayoutA, class Device>
template <class T>
typename std::enable_if<
    std::is_same<T, Kokkos::complex<float>>::value || std::is_same<T, Kokkos::complex<double>>::value, void>::type
SyrTester<ScalarX, tLayoutX, ScalarA, tLayoutA, Device>::compareVanillaAgainstExpected(
    const T& alpha, const _ViewTypeExpected& h_vanilla, const _ViewTypeExpected& h_expected) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  if (_N <= 2) {
    for (int i(0); i < _M; ++i) {
      for (int j(0); j < _N; ++j) {
        std::cout << "h_exp(" << i << "," << j << ")=" << h_expected(i, j) << ", h_van(" << i << "," << j
                  << ")=" << h_vanilla(i, j) << std::endl;
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
          << ", _useUpOption = " << _useUpOption << ": vanilla differs too much from analytical on imag components"
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
template <class ScalarX, class tLayoutX, class ScalarA, class tLayoutA, class Device>
template <class T>
typename std::enable_if<
    !std::is_same<T, Kokkos::complex<float>>::value && !std::is_same<T, Kokkos::complex<double>>::value, void>::type
SyrTester<ScalarX, tLayoutX, ScalarA, tLayoutA, Device>::compareVanillaAgainstExpected(
    const T& alpha, const _ViewTypeExpected& h_vanilla, const _ViewTypeExpected& h_expected) {
  if (_N <= 2) {
    for (int i(0); i < _M; ++i) {
      for (int j(0); j < _N; ++j) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
        std::cout << "h_exp(" << i << "," << j << ")=" << h_expected(i, j) << ", h_van(" << i << "," << j
                  << ")=" << h_vanilla(i, j) << std::endl;
#endif
      }
    }
  }

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
                            << ", _useHermitianOption = " << _useHermitianOption << ", _useUpOption = " << _useUpOption
                            << ": vanilla result is incorrect"
                            << ", numErrors = " << numErrors;
  }
}

// Code for complex values
template <class ScalarX, class tLayoutX, class ScalarA, class tLayoutA, class Device>
template <class T>
typename std::enable_if<
    std::is_same<T, Kokkos::complex<float>>::value || std::is_same<T, Kokkos::complex<double>>::value, void>::type
SyrTester<ScalarX, tLayoutX, ScalarA, tLayoutA, Device>::compareKkSyrAgainstReference(
    const T& alpha, const _HostViewTypeA& h_A, const _ViewTypeExpected& h_reference) {
  if (_N <= 2) {
    for (int i(0); i < _M; ++i) {
      for (int j(0); j < _N; ++j) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
        std::cout << "h_exp(" << i << "," << j << ")=" << h_reference(i, j) << ", h_A(" << i << "," << j
                  << ")=" << h_A(i, j) << std::endl;
#endif
      }
    }
  }

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
        << ", _useUpOption = " << _useUpOption << ": syr result is incorrect on real components"
        << ", numErrorsRealAbs = " << numErrorsRealAbs << ", numErrorsRealRel = " << numErrorsRealRel
        << ", maxErrorRealRel = " << maxErrorRealRel << ", iForMaxErrorRealRel = " << iForMaxErrorRealRel
        << ", jForMaxErrorRealRel = " << jForMaxErrorRealRel << ", h_reference(i,j).real() = "
        << (((_M > 0) && (_N > 0)) ? h_reference(iForMaxErrorRealRel, jForMaxErrorRealRel).real() : 9.999e+99)
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
        << ", _useUpOption = " << _useUpOption << ": syr result is incorrect on imag components"
        << ", numErrorsImagAbs = " << numErrorsImagAbs << ", numErrorsImagRel = " << numErrorsImagRel
        << ", maxErrorImagRel = " << maxErrorImagRel << ", iForMaxErrorImagRel = " << iForMaxErrorImagRel
        << ", jForMaxErrorImagRel = " << jForMaxErrorImagRel << ", h_reference(i,j).imag() = "
        << (((_M > 0) && (_N > 0)) ? h_reference(iForMaxErrorImagRel, jForMaxErrorImagRel).imag() : 9.999e+99)
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
template <class ScalarX, class tLayoutX, class ScalarA, class tLayoutA, class Device>
template <class T>
typename std::enable_if<
    !std::is_same<T, Kokkos::complex<float>>::value && !std::is_same<T, Kokkos::complex<double>>::value, void>::type
SyrTester<ScalarX, tLayoutX, ScalarA, tLayoutA, Device>::compareKkSyrAgainstReference(
    const T& alpha, const _HostViewTypeA& h_A, const _ViewTypeExpected& h_reference) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  if (_N <= 2) {
    for (int i(0); i < _M; ++i) {
      for (int j(0); j < _N; ++j) {
        std::cout << "h_exp(" << i << "," << j << ")=" << h_reference(i, j) << ", h_A(" << i << "," << j
                  << ")=" << h_A(i, j) << std::endl;
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
        << ", _useUpOption = " << _useUpOption << ": syr result is incorrect"
        << ", numErrorsAbs = " << numErrorsAbs << ", numErrorsRel = " << numErrorsRel
        << ", maxErrorRel = " << maxErrorRel << ", iForMaxErrorRel = " << iForMaxErrorRel
        << ", jForMaxErrorRel = " << jForMaxErrorRel << ", h_reference(i,j) = "
        << (((_M > 0) && (_N > 0)) ? h_reference(iForMaxErrorRel, jForMaxErrorRel) : 9.999e+99)
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

template <class ScalarX, class tLayoutX, class ScalarA, class tLayoutA, class Device>
template <class TX>
void SyrTester<ScalarX, tLayoutX, ScalarA, tLayoutA, Device>::callKkSyrAndCompareAgainstExpected(
    const ScalarA& alpha, TX& x, view_stride_adapter<_ViewTypeA, false>& A, const _ViewTypeExpected& h_expected,
    const std::string& situation) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "In Test_Blas2_syr, '" << situation << "', alpha = " << alpha << std::endl;
  Kokkos::printf(
      "In Test_Blas2_syr.hpp, right before calling KokkosBlas::syr(): "
      "ViewTypeA = %s, _kkSyrShouldThrowException=%d\n",
      typeid(_ViewTypeA).name(), _kkSyrShouldThrowException);
#endif
  std::string mode = _useHermitianOption ? "H" : "T";
  std::string uplo = _useUpOption ? "U" : "L";
  bool gotStdException(false);
  bool gotUnknownException(false);
  try {
    KokkosBlas::syr(mode.c_str(), uplo.c_str(), alpha, x, A.d_view);
  } catch (const std::exception& e) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
    std::cout << "In Test_Blas2_syr, '" << situation << "': caught exception, e.what() = " << e.what() << std::endl;
#endif
    gotStdException = true;
  } catch (...) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
    std::cout << "In Test_Blas2_syr, '" << situation << "': caught unknown exception" << std::endl;
#endif
    gotUnknownException = true;
  }

  EXPECT_EQ(gotUnknownException, false) << "Failed test, '" << situation
                                        << "': unknown exception should not have happened";

  EXPECT_EQ(gotStdException, _kkSyrShouldThrowException)
      << "Failed test, '" << situation << "': kk syr() should" << (_kkSyrShouldThrowException ? " " : " not ")
      << "have thrown a std::exception";

  if ((gotStdException == false) && (gotUnknownException == false)) {
    Kokkos::deep_copy(A.h_base, A.d_base);
    this->compareKkSyrAgainstReference(alpha, A.h_view, h_expected);
  }
}

template <class ScalarX, class tLayoutX, class ScalarA, class tLayoutA, class Device>
template <class TX>
void SyrTester<ScalarX, tLayoutX, ScalarA, tLayoutA, Device>::callKkGerAndCompareKkSyrAgainstIt(
    const ScalarA& alpha, TX& x, view_stride_adapter<_ViewTypeA, false>& org_A, const _HostViewTypeA& h_A_syr,
    const std::string& situation) {
  view_stride_adapter<_ViewTypeA, false> A_ger("A_ger", _M, _N);
  Kokkos::deep_copy(A_ger.d_base, org_A.d_base);

  // ********************************************************************
  // Call ger()
  // ********************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "In Test_Blas2_syr, '" << situation << "', alpha = " << alpha << std::endl;
  Kokkos::printf(
      "In Test_Blas2_syr.hpp, right before calling KokkosBlas::ger(): "
      "ViewTypeA = %s, _kkGerShouldThrowException=%d\n",
      typeid(_ViewTypeA).name(), _kkGerShouldThrowException);
#endif
  std::string mode = _useHermitianOption ? "H" : "T";
  bool gotStdException(false);
  bool gotUnknownException(false);
  try {
    KokkosBlas::ger(mode.c_str(), alpha, x, x, A_ger.d_view);
  } catch (const std::exception& e) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
    std::cout << "In Test_Blas2_syr, '" << situation << "', ger() call: caught exception, e.what() = " << e.what()
              << std::endl;
#endif
    gotStdException = true;
  } catch (...) {
#ifdef HAVE_KOKKOSKERNELS_DEBUG
    std::cout << "In Test_Blas2_syr, '" << situation << "', ger() call: caught unknown exception" << std::endl;
#endif
    gotUnknownException = true;
  }

  EXPECT_EQ(gotUnknownException, false) << "Failed test, '" << situation
                                        << "': unknown exception should not have happened for ger() call";

  EXPECT_EQ(gotStdException, false) << "Failed test, '" << situation
                                    << "': kk ger() should not have thrown a std::exception";

  // ********************************************************************
  // Prepare h_ger_reference to be compared against h_A_syr
  // ********************************************************************
  view_stride_adapter<_ViewTypeExpected, true> h_ger_reference("h_ger_reference", _M, _N);
  Kokkos::deep_copy(h_ger_reference.d_base, A_ger.d_base);

  std::string uplo = _useUpOption ? "U" : "L";
  for (int i = 0; i < _M; ++i) {
    for (int j = 0; j < _N; ++j) {
      if (((_useUpOption == true) && (i <= j)) || ((_useUpOption == false) && (i >= j))) {
        // Keep h_ger_reference as already computed
      } else {
        h_ger_reference.d_view(i, j) = org_A.h_view(i, j);
      }
    }
  }
  if (_useHermitianOption && _A_is_complex) {
    for (int i(0); i < _N; ++i) {
      h_ger_reference.d_view(i, i) = 0.5 * (h_ger_reference.d_view(i, i) + _KAT_A::conj(h_ger_reference.d_view(i, i)));
    }
  }

  // ********************************************************************
  // Compare
  // ********************************************************************
  this->compareKkSyrAgainstReference(alpha, h_A_syr, h_ger_reference.d_view);
}

}  // namespace Test

template <class ScalarX, class ScalarA, class Device>
#ifdef HAVE_KOKKOSKERNELS_DEBUG
int test_syr(const std::string& caseName) {
  Kokkos::printf(
      "+======================================================================="
      "===\n");
  Kokkos::printf("Starting %s ...\n", caseName.c_str());
#else
int test_syr(const std::string& /*caseName*/) {
#endif
  bool xBool = std::is_same<ScalarX, float>::value || std::is_same<ScalarX, double>::value ||
               std::is_same<ScalarX, Kokkos::complex<float>>::value ||
               std::is_same<ScalarX, Kokkos::complex<double>>::value;
  bool aBool = std::is_same<ScalarA, float>::value || std::is_same<ScalarA, double>::value ||
               std::is_same<ScalarA, Kokkos::complex<float>>::value ||
               std::is_same<ScalarA, Kokkos::complex<double>>::value;
  bool useAnalyticalResults = xBool && aBool;

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  Kokkos::printf(
      "+-----------------------------------------------------------------------"
      "---\n");
  Kokkos::printf("Starting %s for LAYOUTLEFT ...\n", caseName.c_str());
#endif
  if (true) {
    Test::SyrTester<ScalarX, Kokkos::LayoutLeft, ScalarA, Kokkos::LayoutLeft, Device> tester;
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
    Test::SyrTester<ScalarX, Kokkos::LayoutRight, ScalarA, Kokkos::LayoutRight, Device> tester;
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
  Kokkos::printf("Finished %s for LAYOUTRIGHT\n", caseName.c_str());
  Kokkos::printf(
      "+-----------------------------------------------------------------------"
      "---\n");
#endif
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  Kokkos::printf(
      "+-----------------------------------------------------------------------"
      "---\n");
  Kokkos::printf("Starting %s for LAYOUTSTRIDE ...\n", caseName.c_str());
#endif
  if (true) {
    Test::SyrTester<ScalarX, Kokkos::LayoutStride, ScalarA, Kokkos::LayoutStride, Device> tester;
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
    Test::SyrTester<ScalarX, Kokkos::LayoutStride, ScalarA, Kokkos::LayoutRight, Device> tester;
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
    Test::SyrTester<ScalarX, Kokkos::LayoutLeft, ScalarA, Kokkos::LayoutRight, Device> tester;
    tester.test(1024, 0);
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
TEST_F(TestCategory, syr_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::syr_float");
  test_syr<float, float, TestDevice>("test case syr_float");
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, syr_complex_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::syr_complex_float");
  test_syr<Kokkos::complex<float>, Kokkos::complex<float>, TestDevice>("test case syr_complex_float");
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, syr_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::syr_double");
  test_syr<double, double, TestDevice>("test case syr_double");
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, syr_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::syr_complex_double");
  test_syr<Kokkos::complex<double>, Kokkos::complex<double>, TestDevice>("test case syr_complex_double");
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, syr_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::syr_int");
  test_syr<int, int, TestDevice>("test case syr_int");
  Kokkos::Profiling::popRegion();
}
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
TEST_F(TestCategory, syr_int_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::syr_int_float");
  test_syr<int, float, TestDevice>("test case syr_int_float");
  Kokkos::Profiling::popRegion();
}
#endif
