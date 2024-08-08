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

/// \file ArithTraitsTest.hpp
/// \brief Templated test for Kokkos::ArithTraits
///
/// This header file is an implementation detail of the tests for
/// Kokkos::ArithTraits.  Users must not rely on it existing,
/// or on its contents.  This header file should <i>not</i> be
/// installed with Kokkos' other header files.
///
/// On the other hand, this header file does give examples of how to
/// use Kokkos::ArithTraits, so it may be useful for users to
/// read it.

#ifndef KOKKOS_ARITHTRAITSTEST_HPP
#define KOKKOS_ARITHTRAITSTEST_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_ArithTraits.hpp"
#include <limits>    // std::numeric_limits
#include <typeinfo>  // typeid (T)
#include <cstdio>

#define FAILURE()                                                        \
  {                                                                      \
    Kokkos::printf("%s:%s:%d: Failure\n", __FILE__, __func__, __LINE__); \
    success = 0;                                                         \
  }

#if 0
#define TRACE() Kokkos::printf("%s:%s:%d: Trace\n", __FILE__, __func__, __LINE__);
#else
#define TRACE()
#endif

namespace {
// Whether Kokkos::ArithTraits<ScalarType> implements
// transcendental functions.  These include sqrt, pow, log, and
// log10.
template <class ScalarType>
struct HasTranscendentals {
  static const bool value = false;
};

template <>
struct HasTranscendentals<float> {
  static const bool value = true;
};

template <>
struct HasTranscendentals<double> {
  static const bool value = true;
};

template <>
struct HasTranscendentals<long double> {
  static const bool value = true;
};

// template<>
// struct HasTranscendentals< ::Kokkos::complex<float> > {
//   static const bool value = true;
// };

// template<>
// struct HasTranscendentals< ::Kokkos::complex<double> > {
//   static const bool value = true;
// };

// template<>
// struct HasTranscendentals< ::Kokkos::complex<long double> > {
//   static const bool value = true;
// };

}  // namespace

/// \class ArithTraitsTesterBase
/// \brief Base class providing tests for Kokkos::ArithTraits
/// \tparam ScalarType Any type for which Kokkos::ArithTraits
///   has a specialization, and which can be executed on the parallel
///   device.
/// \tparam DeviceType A Kokkos parallel device type.
///
/// This class is really an implementation detail of ArithTraitsTester
/// (see below).  ArithTraitsTester works by inheriting "hooks" from
/// ArithTraitsTesterBase and the chain of subclasses in between.
/// ArithTraitsTesterBase provides basic tests that work for all
/// <tt>ScalarType</tt>, and the subclasses provide additional tests
/// relevant to things like complex-valued types or floating-point
/// types.
///
/// This class provides a Kokkos reduction operator for testing
/// Kokkos::ArithTraits.  This test works for any type
/// <tt>ScalarType</tt> for which Kokkos::ArithTraits has a
/// specialization, and which can be executed on the parallel device.
///
/// The tests include those suitable for execution on the parallel
/// device (operator()) and those suitable for execution on the host
/// (testHost()).  The device-based test is a reduction over redundant
/// executions of the test.  All redundant executions must return
/// '1' (passed).
template <class ScalarType, class DeviceType>
class ArithTraitsTesterBase {
 public:
  typedef typename DeviceType::execution_space execution_space;
  typedef typename execution_space::size_type size_type;
  //! Type of the result of the reduction.
  typedef int value_type;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterBase() {}

  /// \brief Set the initial value (\c 1) of the reduction.
  ///
  /// Subclasses need not and must not override this method.
  KOKKOS_INLINE_FUNCTION void init(value_type& dst) const { dst = 1; }

  /// \brief Combine two intermediate reduction results into \c dst.
  ///
  /// Subclasses need not and must not override this method.
  KOKKOS_INLINE_FUNCTION void join(value_type& dst, const value_type& src) const {
    dst = dst && src;
    // dst = 1;
  }

  /// \brief The "parallel for" part of the reduction.
  ///
  /// This is the method that actually runs the tests on the device.
  /// It runs through a sequence of tests, and produces a \c 1
  /// result if all the tests pass.
  ///
  /// Subclasses must override this to implement their own tests.
  /// They must always call their parent class' version.  Refer to the
  /// implementations of operator() in ArithTraitsTesterComplexBase
  /// for examples.  Subclass' implementations must ignore \c work,
  /// and set the \c dst argument to the logical AND of all the tests'
  /// results.
  ///
  /// \param iwork [in] Ignored.
  /// \param dst [in/out] On input: The result of any tests run thus
  ///   far.  On output: The result of the tests run in this method.
  ///   The result of more than one test is the logical AND of each
  ///   test's result.
  KOKKOS_INLINE_FUNCTION void operator()(size_type iwork, value_type& dst) const {
    TRACE();
    typedef Kokkos::ArithTraits<ScalarType> AT;
    (void)iwork;  // not using this argument
    int success = 1;

    // Make sure that the typedef exists.
    typedef typename AT::mag_type mag_type;

    // mfh 14 Feb 2014: In order to avoid a warning for an unused by
    // declared typedef, we declare an instance of mag_type, and mark
    // it with "(void)" to prevent a warning for the unused variable.
    {
      mag_type thing;
      (void)thing;
    }

    // ArithTraits should not even compile if it's not specialized for
    // T, but we check for this int constant for compatibility with
    // std::numeric_limits.
    if (!AT::is_specialized) {
      Kokkos::printf("! AT::is_specialized\n");
      FAILURE();
    }

    // It's OK to refer to std::numeric_limits constants in a device
    // function, just not to its class methods (which are not marked
    // as device functions).
    if (AT::is_integer != std::numeric_limits<ScalarType>::is_integer) {
      Kokkos::printf("AT::is_integer not same as numeric_limits\n");
      FAILURE();
    }
    if (AT::is_exact != std::numeric_limits<ScalarType>::is_exact) {
      Kokkos::printf("AT::is_exact not same as numeric_limits\n");
      FAILURE();
    }

    const ScalarType zero = AT::zero();
    const ScalarType one  = AT::one();

    // Test properties of the arithmetic and multiplicative identities.
    if (zero + zero != zero) {
      Kokkos::printf("0 + 0 != 0\n");
      FAILURE();
    }
    if (zero + one != one) {
      Kokkos::printf("0 + 1 != 1\n");
      FAILURE();
    }
    if (one - one != zero) {
      Kokkos::printf("1 - 1 != 0\n");
      FAILURE();
    }
    // This is technically 1 even of Z_2, since in that field, one
    // is its own inverse (so -one == one).
    if ((one + one) - one != one) {
      Kokkos::printf("(1 + 1) - 1 != 1\n");
      FAILURE();
    }

    if (AT::abs(zero) != zero) {
      Kokkos::printf("AT::abs(0) != 0\n");
      FAILURE();
    }
    if (AT::abs(one) != one) {
      Kokkos::printf("AT::abs(1) != 1\n");
      FAILURE();
    }
    if (AT::is_signed && AT::abs(-one) != one) {
      Kokkos::printf("AT::is_signed and AT::abs(-1) != 1\n");
      FAILURE();
    }
    // Need enable_if to test whether T can be compared using <=.
    // However, mag_type should always be comparable using <=.
    //
    // These are very mild ordering properties.
    // They should work even for a set only containing zero.
    if (AT::abs(zero) > AT::abs(AT::max())) {
      Kokkos::printf("AT::abs(0) > AT::abs (AT::max ())\n");
      FAILURE();
    }

    dst = dst && success;
  }

 protected:
  /// \brief Hook for subclasses to add their own host-based tests.
  ///
  /// We use this to add complex-arithmetic tests, if appropriate for
  /// \c ScalarType.  You may use it for other tests that are specific
  /// to \c ScalarType.
  ///
  /// The default implementation does nothing.  (That's what makes
  /// this a "hook.")
  ///
  /// \return \c 1 if all tests succeeded, else \c 0.
  int testHostImpl(std::ostream& /*out*/) const {
    return 1;  // there are no tests, so trivially, all the tests pass
  }

 public:
  /// \brief Run the tests on the host.
  ///
  /// This method only works on the host.  It's helpful for debugging,
  /// because it prints a message (to the given output stream) for
  /// each test that fails.
  ///
  /// \param out [out] Output stream to which to print error messages.
  ///
  /// \return \c 1 if all the tests pass, else \c 0.
  int testHost(std::ostream& out) const {
    typedef Kokkos::ArithTraits<ScalarType> AT;
    using std::endl;
    int success = 1;

    // Make sure that the typedef exists.
    typedef typename AT::mag_type mag_type;

    // mfh 14 Feb 2014: In order to avoid a warning for an unused by
    // declared typedef, we declare an instance of mag_type, and mark
    // it with "(void)" to prevent a warning for the unused variable.
    {
      mag_type thing;
      (void)thing;
    }

    // ArithTraits should not even compile if it's not specialized for
    // T, but we check for this int constant for compatibility with
    // std::numeric_limits.
    if (!AT::is_specialized) {
      out << "ArithTraits is not specialized for T" << endl;
      FAILURE();
    }

    if (AT::is_integer != std::numeric_limits<ScalarType>::is_integer) {
      out << "AT::is_integer != std::numeric_limits<ScalarType>::is_integer" << endl;
      FAILURE();
    }

    if (AT::is_exact != std::numeric_limits<ScalarType>::is_exact) {
      out << "AT::is_exact != std::numeric_limits<ScalarType>::is_exact" << endl;
      FAILURE();
    }

    const ScalarType zero = AT::zero();
    const ScalarType one  = AT::one();
    // Test properties of the arithmetic and multiplicative identities.

    if (zero + zero != zero) {
      out << "zero + zero != zero" << endl;
      FAILURE();
    }
    if (zero + one != one) {
      out << "zero + one != one" << endl;
      FAILURE();
    }
    if (one - one != zero) {
      out << "one - one != zero" << endl;
      FAILURE();
    }
    // This is technically 1 even of Z_2, since in that field, one
    // is its own inverse (so -one == one).
    if ((one + one) - one != one) {
      out << "(one + one) - one != one" << endl;
      FAILURE();
    }

    if (AT::abs(zero) != zero) {
      out << "AT::abs (zero) != zero" << endl;
      FAILURE();
    }
    if (AT::abs(one) != one) {
      out << "AT::abs (one) != one" << endl;
      FAILURE();
    }
    if (AT::is_signed) {
      if (AT::abs(-one) != one) {
        out << "AT::abs (-one) != one" << endl;
        FAILURE();
      }
    }
    // Need enable_if to test whether T can be compared using <=.
    // However, mag_type should always be comparable using <=.
    //
    // // These are very mild ordering properties.
    // // They should work even for a set only containing zero.
    if (AT::abs(zero) > AT::abs(AT::max())) {
      out << "AT::abs (zero) > AT::abs (AT::max ())" << endl;
      FAILURE();
    }

    if (AT::has_infinity) {
// Compiler intrinsic casts from inf of type half_t / bhalf_t to inf
// of type float in CUDA, SYCL and HIP do not work yet.
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ENABLE_HIP)
      namespace KE = Kokkos::Experimental;
      if constexpr (!std::is_same<ScalarType, KE::half_t>::value && !std::is_same<ScalarType, KE::bhalf_t>::value) {
#else
      {
#endif  // KOKKOS_ENABLE_CUDA || KOKKOS_ENABLE_SYCL || KOKKOS_ENABLE_HIP
        if (!AT::isInf(AT::infinity())) {
          out << "AT::isInf (inf) != true" << endl;
          FAILURE();
        }
      }
    }
    if (!std::is_same<ScalarType, decltype(AT::infinity())>::value) {
      std::cout << "AT::infinity() return value has wrong type" << endl;
      FAILURE();
    }

    // Run the parent class' remaining tests, if any.
    const int parentSuccess = testHostImpl(out);
    success                 = success && parentSuccess;

    return success;
  }
};

/// \class ArithTraitsTesterTranscendentalBase
/// \brief Base class of ArithTraitsTester that exercises
///   transcendental functions, if and only if ArithTraits<ScalarType>
///   implements them.
/// \tparam ScalarType Any type for which Kokkos::ArithTraits
///   implements transcendental functions, along with the requirements
///   imposed by ArithTraitsTesterBase.
/// \tparam DeviceType A Kokkos parallel device type.
/// \tparam hasTranscendentals Whether ArithTraits<ScalarType>
///   implements transcendental functions.
///
/// Some tests will be executed whether or not ArithTraits<ScalarType>
/// implements transcendental functions, but the specific tests that
/// are run will depend on \c ScalarType.
template <class ScalarType, class DeviceType,
          const int has_transcendentals = (HasTranscendentals<ScalarType>::value ? 1 : 0)>
class ArithTraitsTesterTranscendentalBase : public ArithTraitsTesterBase<ScalarType, DeviceType> {
 private:
  //! The base class of this class.
  typedef ArithTraitsTesterBase<ScalarType, DeviceType> base_type;

 public:
  typedef DeviceType execution_space;
  typedef typename execution_space::size_type size_type;
  //! Type of the result of the reduction.
  typedef int value_type;

  /// \brief The "parallel for" part of the reduction.
  ///
  /// See comments of ArithTraitsTesterBase's operator().
  KOKKOS_INLINE_FUNCTION void operator()(size_type iwork, value_type& dst) const;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterTranscendentalBase();

 protected:
  // The host hook gets implemented in the "transcendental functions
  // are implemented" specialization of this class.
  virtual int testHostImpl(std::ostream& out) const;
};

//
// Specialization of ArithTraitsTesterTranscendentalBase when
// ArithTraits<ScalarType> does NOT implement transcendentals.
//
template <class ScalarType, class DeviceType>
class ArithTraitsTesterTranscendentalBase<ScalarType, DeviceType, 0>
    : public ArithTraitsTesterBase<ScalarType, DeviceType> {
 private:
  //! The base class of this class.
  typedef ArithTraitsTesterBase<ScalarType, DeviceType> base_type;

 public:
  typedef typename DeviceType::execution_space execution_space;
  typedef typename execution_space::size_type size_type;
  //! Type of the result of the reduction.
  typedef int value_type;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterTranscendentalBase() {}

  KOKKOS_INLINE_FUNCTION void operator()(size_type iwork, value_type& dst) const {
    TRACE();
    // typedef Kokkos::ArithTraits<ScalarType> AT;
    (void)iwork;  // forestall compiler warning for unused variable
    int success = 1;

    if (HasTranscendentals<ScalarType>::value) {
      FAILURE();
    }

    // Call the base class' implementation.  Every subclass'
    // implementation of operator() must do this, in order to include
    // the parent class' tests.
    int baseResult = 1;
    base_type::operator()(iwork, baseResult);
    success = success && baseResult;

    dst = dst && success;
  }

 protected:
  virtual int testHostImpl(std::ostream& out) const {
    using std::endl;
    // typedef Kokkos::ArithTraits<ScalarType> AT;
    int success = 1;

    if (HasTranscendentals<ScalarType>::value) {
      out << "HasTranscendentals<T>::value is true" << endl;
      FAILURE();
    }

    // Call the base class' implementation.  Every subclass'
    // implementation of testHostImpl() should (must) do this, in
    // order to include the parent class' tests.  In the case of this
    // particular class, the base class' implementation doesn't do
    // anything, but that's OK.
    const int parentSuccess = base_type::testHostImpl(out);
    success                 = success && parentSuccess;

    return success;
  }
};

//
// Specialization of ArithTraitsTesterTranscendentalBase when
// ArithTraits<ScalarType> DOES implement transcendentals.
//
template <class ScalarType, class DeviceType>
class ArithTraitsTesterTranscendentalBase<ScalarType, DeviceType, 1>
    : public ArithTraitsTesterBase<ScalarType, DeviceType> {
 private:
  //! The base class of this class.
  typedef ArithTraitsTesterBase<ScalarType, DeviceType> base_type;

  KOKKOS_INLINE_FUNCTION
  bool equal(const ScalarType& a, const ScalarType& b) const {
    if (b != Kokkos::ArithTraits<ScalarType>::zero()) {
      if (a > b)
        return (a - b) / b < 2 * Kokkos::ArithTraits<ScalarType>::epsilon();
      else
        return (b - a) / b < 2 * Kokkos::ArithTraits<ScalarType>::epsilon();
    } else {
      if (a > b)
        return (a - b) < 2 * Kokkos::ArithTraits<ScalarType>::epsilon();
      else
        return (b - a) < 2 * Kokkos::ArithTraits<ScalarType>::epsilon();
    }
  }

 public:
  typedef typename DeviceType::execution_space execution_space;
  typedef typename execution_space::size_type size_type;
  //! Type of the result of the reduction.
  typedef int value_type;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterTranscendentalBase() {}

  KOKKOS_INLINE_FUNCTION void operator()(size_type iwork, value_type& dst) const {
    TRACE();
    typedef Kokkos::ArithTraits<ScalarType> AT;
    (void)iwork;  // forestall compiler warning for unused variable
    int success = 1;

    if (!HasTranscendentals<ScalarType>::value) {
      FAILURE();
    }

    const ScalarType zero        = AT::zero();
    const ScalarType one         = AT::one();
    const ScalarType two         = one + one;
    const ScalarType three       = one + one + one;
    const ScalarType four        = two * two;
    const ScalarType five        = four + one;
    const ScalarType six         = three * two;
    const ScalarType seven       = four + three;
    const ScalarType eight       = four * two;
    const ScalarType nine        = eight + one;
    const ScalarType eleven      = five + six;
    const ScalarType twentySeven = nine * three;
    const ScalarType thirtySix   = six * six;
    const ScalarType fortyTwo    = six * seven;
    const ScalarType sixtyThree  = eight * eight - one;
    const ScalarType sixtyFour   = eight * eight;
    // max char value, for 8-bit char
    const ScalarType oneTwentySeven = sixtyFour + sixtyThree;

    ScalarType result;

    // This fails inexplicably for complex numbers on gcc 4.2.1 on Mac.
    if (!AT::is_complex) {
      result = AT::pow(two, three);
      if (!equal(result, eight)) {
        Kokkos::printf("AT::pow(2,3) != 8\n");
        FAILURE();
      }
    }
    if (!equal(AT::pow(three, zero), one)) {
      Kokkos::printf("AT::pow(3,0) != 1\n");
      FAILURE();
    }
    if (!equal(AT::pow(three, one), three)) {
      Kokkos::printf("AT::pow(3,1) != 3\n");
      FAILURE();
    }
    if (!equal(AT::pow(three, two), nine)) {
      Kokkos::printf("AT::pow(3,2) != 9\n");
      FAILURE();
    }

    // This fails inexplicably for complex numbers on gcc 4.2.1 on Mac.
    if (!AT::is_complex) {
      result = AT::pow(three, three);
      if (!equal(result, twentySeven)) {
        Kokkos::printf("AT::pow(3,3) != 27\n");
        FAILURE();
      }
    }

    // These fail inexplicably for complex numbers on gcc 4.2.1 on Mac.
    if (AT::is_signed && !AT::is_complex) {
      result = AT::pow(-three, one);
      if (!equal(result, -three)) {
        Kokkos::printf("AT::pow(-3,1) != -3\n");
        FAILURE();
      }
      result = AT::pow(-three, two);
      if (!equal(result, nine)) {
        Kokkos::printf("AT::pow(-3,2) != 9\n");
        FAILURE();
      }
      result = AT::pow(-three, three);
      if (!equal(result, -twentySeven)) {
        Kokkos::printf("AT::pow(-3,3) != 27\n");
        FAILURE();
      }
    }

    if (!equal(AT::sqrt(zero), zero)) {
      Kokkos::printf("AT::sqrt(0) != 0\n");
      FAILURE();
    }
    if (!equal(AT::sqrt(one), one)) {
      Kokkos::printf("AT::sqrt(1) != 1\n");
      FAILURE();
    }
    if (!equal(AT::sqrt(thirtySix), six)) {
      Kokkos::printf("AT::sqrt(36) != 6\n");
      FAILURE();
    }
    if (!equal(AT::sqrt(sixtyFour), eight)) {
      Kokkos::printf("AT::sqrt(64) != 8\n");
      FAILURE();
    }
    if (AT::is_integer) {
      if (!equal(AT::sqrt(fortyTwo), six)) {
        Kokkos::printf("AT:sqrt(42) != 6\n");
        FAILURE();
      }
      if (!equal(AT::sqrt(oneTwentySeven), eleven)) {
        Kokkos::printf("AT::sqrt(127) != 11\n");
        FAILURE();
      }
    }

    if (!equal(AT::cbrt(zero), zero)) {
      Kokkos::printf("AT::cbrt(0) != 0\n");
      FAILURE();
    }
    if (!equal(AT::cbrt(one), one)) {
      Kokkos::printf("AT::cbrt(1) != 1\n");
      FAILURE();
    }
    if (!equal(AT::cbrt(twentySeven), three)) {
      Kokkos::printf("AT::cbrt(27) != 3\n");
      FAILURE();
    }
    if (!equal(AT::cbrt(sixtyFour), four)) {
      Kokkos::printf("AT::cbrt(64) != 4\n");
      FAILURE();
    }
    if (AT::is_integer) {
      if (!equal(AT::cbrt(fortyTwo), three)) {
        Kokkos::printf("AT:cbrt(42) != 3\n");
        FAILURE();
      }
      if (!equal(AT::cbrt(oneTwentySeven), five)) {
        Kokkos::printf("AT::cbrt(127) != 5\n");
        FAILURE();
      }
    }

    if (!equal(AT::exp(zero), one)) {
      Kokkos::printf("AT::cbrt(0) != 1\n");
      FAILURE();
    }
    if (AT::is_complex) {
      const ScalarType val = two;  //(two.real(), two.real());
      if (!equal(AT::conj(AT::exp(val)), AT::exp(AT::conj(val)))) {
        Kokkos::printf("AT::conj(exp(complex(2,2))) != AT::exp(conj(complex(2,2)))\n");
        FAILURE();
      }
    }
    if (!equal(AT::log(one), zero)) {
      Kokkos::printf("AT::log(1) != 0\n");
      FAILURE();
    }
    if (!equal(AT::log10(one), zero)) {
      Kokkos::printf("AT::log10(1) != 0\n");
      FAILURE();
    }

    if (AT::is_complex) {
      ScalarType val     = two;  //(two, two);
      const auto val_sin = AT::sin(val);
      const auto val_cos = AT::cos(val);
      if (!equal(val_sin * val_sin + val_cos * val_cos, one)) {
        Kokkos::printf("AT(complex):: sin(val)*sin(val) + cos(val)*cos(val) != 1\n");
        FAILURE();
      }
      if (!equal(val_sin / val_cos, AT::tan(val))) {
        Kokkos::printf("AT(complex):: sin(val)/cos(val) != AT(real)::tan(val)\n");
        FAILURE();
      }
    } else {
      ScalarType val     = three;
      const auto val_sin = AT::sin(val);
      const auto val_cos = AT::cos(val);
      if (!equal(val_sin * val_sin + val_cos * val_cos, one)) {
        Kokkos::printf("AT(real):: sin(val)*sin(val) + cos(a)*cos(a) != 1\n");
        FAILURE();
      }
      if (!equal(val_sin / val_cos, AT::tan(val))) {
        Kokkos::printf("AT(real):: sin(val)/cos(val) != AT(real)::tan(val)\n");
        FAILURE();
      }
    }

    if (!equal(AT::asin(AT::sin(one)), one)) {
      Kokkos::printf("AT::asin(sin(1)) != 1\n");
      FAILURE();
    }
    if (!equal(AT::acos(AT::cos(one)), one)) {
      Kokkos::printf("AT::acos(cos(1)) != 1\n");
      FAILURE();
    }
    if (!equal(AT::atan(AT::tan(one)), one)) {
      Kokkos::printf("AT::atan(tan(1)) != 1\n");
      FAILURE();
    }

    // Call the base class' implementation.  Every subclass'
    // implementation of operator() must do this, in order to include
    // the parent class' tests.
    int baseResult = 1;
    base_type::operator()(iwork, baseResult);
    success = success && baseResult;

    dst = dst && success;
  }

 protected:
  virtual int testHostImpl(std::ostream& out) const {
    using std::endl;
    typedef Kokkos::ArithTraits<ScalarType> AT;
    int success = 1;

    if (!HasTranscendentals<ScalarType>::value) {
      out << "HasTranscendentals<T>::value is false" << endl;
      FAILURE();
    }

    const ScalarType zero        = AT::zero();
    const ScalarType one         = AT::one();
    const ScalarType two         = one + one;
    const ScalarType three       = one + one + one;
    const ScalarType four        = two * two;
    const ScalarType five        = four + one;
    const ScalarType six         = three * two;
    const ScalarType seven       = four + three;
    const ScalarType eight       = four * two;
    const ScalarType nine        = eight + one;
    const ScalarType eleven      = five + six;
    const ScalarType twentySeven = nine * three;
    const ScalarType thirtySix   = six * six;
    const ScalarType fortyTwo    = six * seven;
    const ScalarType sixtyThree  = eight * eight - one;
    const ScalarType sixtyFour   = eight * eight;
    // max char value, for 8-bit char
    const ScalarType oneTwentySeven = sixtyFour + sixtyThree;

    ScalarType result;

    // This fails inexplicably for complex numbers on gcc 4.2.1 on Mac.
    if (!AT::is_complex) {
      result = AT::pow(two, three);
      if (result != eight) {
        out << "AT::pow (two, three) != eight" << endl;
        FAILURE();
      }
    }
    if (AT::pow(three, zero) != one) {
      out << "AT::pow (three, zero) != one" << endl;
      FAILURE();
    }
    if (AT::pow(three, one) != three) {
      out << "AT::pow (three, one) != three" << endl;
      FAILURE();
    }
    if (AT::pow(three, two) != nine) {
      out << "AT::pow (three, two) != nine" << endl;
      FAILURE();
    }

    // This fails inexplicably for complex numbers on gcc 4.2.1 on Mac.
    if (!AT::is_complex) {
      result = AT::pow(three, three);
      if (result != twentySeven) {
        out << "AT::pow (three, three) = " << result << " != twentySeven = " << twentySeven << endl;
        FAILURE();
      }
    }

    // These fail inexplicably for complex numbers on gcc 4.2.1 on Mac.
    if (AT::is_signed && !AT::is_complex) {
      result = AT::pow(-three, one);
      if (result != -three) {
        out << "AT::pow (-three, one) = " << result << " != -three = " << -three << endl;
        FAILURE();
      }
      result = AT::pow(-three, two);
      if (result != nine) {
        out << "AT::pow (-three, two) = " << result << " != nine = " << nine << endl;
        FAILURE();
      }
      result = AT::pow(-three, three);
      if (result != -twentySeven) {
        out << "AT::pow (-three, three) = " << result << " != -twentySeven = " << twentySeven << endl;
        FAILURE();
      }
    }

    if (AT::sqrt(zero) != zero) {
      out << "AT::sqrt (zero) != zero" << endl;
      FAILURE();
    }
    if (AT::sqrt(one) != one) {
      out << "AT::sqrt (one) != one" << endl;
      FAILURE();
    }
    if (AT::sqrt(thirtySix) != six) {
      out << "AT::sqrt (thirtySix) != six" << endl;
      FAILURE();
    }
    if (AT::sqrt(sixtyFour) != eight) {
      out << "AT::sqrt (sixtyFour) != eight" << endl;
      FAILURE();
    }
    if (AT::is_integer) {
      if (AT::sqrt(fortyTwo) != six) {
        out << "AT::sqrt (fortyTwo) != six" << endl;
        FAILURE();
      }
      if (AT::sqrt(oneTwentySeven) != eleven) {
        out << "AT::sqrt (oneTwentySeven) != eleven" << endl;
        FAILURE();
      }
    }

    if (!equal(AT::cbrt(zero), zero)) {
      Kokkos::printf("AT::cbrt(0) != 0\n");
      FAILURE();
    }
    if (!equal(AT::cbrt(one), one)) {
      Kokkos::printf("AT::cbrt(1) != 1\n");
      FAILURE();
    }
    if (!equal(AT::cbrt(twentySeven), three)) {
      Kokkos::printf("AT::cbrt(27) != 3\n");
      FAILURE();
    }
    if (!equal(AT::cbrt(sixtyFour), four)) {
      Kokkos::printf("AT::cbrt(64) != 4\n");
      FAILURE();
    }
    if (AT::is_integer) {
      if (!equal(AT::cbrt(fortyTwo), three)) {
        Kokkos::printf("AT:cbrt(42) != 3\n");
        FAILURE();
      }
      if (!equal(AT::cbrt(oneTwentySeven), five)) {
        Kokkos::printf("AT::cbrt(127) != 5\n");
        FAILURE();
      }
    }

    if (!equal(AT::exp(zero), one)) {
      Kokkos::printf("AT::cbrt(0) != 1\n");
      FAILURE();
    }
    if (AT::is_complex) {
      const ScalarType val = two;  //(two.real(), two.real());
      if (!equal(AT::conj(AT::exp(val)), AT::exp(AT::conj(val)))) {
        Kokkos::printf("AT::conj(exp(complex(2,0))) != AT::exp(conj(complex(2,0)))\n");
        FAILURE();
      }
    }
    if (AT::log(one) != zero) {
      out << "AT::log (one) != zero" << endl;
      FAILURE();
    }
    if (AT::log10(one) != zero) {
      out << "AT::log10 (one) != zero" << endl;
      FAILURE();
    }

    if (AT::is_complex) {
      const ScalarType val = two;  // (two.real(), two.real());
      const auto val_sin   = AT::sin(val);
      const auto val_cos   = AT::cos(val);
      if (!equal(val_sin * val_sin + val_cos * val_cos, one)) {
        Kokkos::printf("AT(complex):: sin(val)*sin(val) + cos(val)*cos(val) != 1\n");
        FAILURE();
      }
      if (!equal(val_sin / val_cos, AT::tan(val))) {
        Kokkos::printf("AT(complex):: sin(val)/cos(val) != AT(real)::tan(val)\n");
        FAILURE();
      }
    } else {
      const ScalarType val = three;
      const auto val_sin   = AT::sin(val);
      const auto val_cos   = AT::cos(val);
      if (!equal(val_sin * val_sin + val_cos * val_cos, one)) {
        Kokkos::printf("AT(real):: sin(val)*sin(val) + cos(a)*cos(a) != 1\n");
        FAILURE();
      }
      if (!equal(val_sin / val_cos, AT::tan(val))) {
        Kokkos::printf("AT(real):: sin(val)/cos(val) != AT(real)::tan(val)\n");
        FAILURE();
      }
    }

    if (!equal(AT::asin(AT::sin(three)), three)) {
      Kokkos::printf("AT::asin(sin(3)) != 3\n");
      FAILURE();
    }
    if (!equal(AT::acos(AT::cos(three)), three)) {
      Kokkos::printf("AT::acos(cos(3)) != 3\n");
      FAILURE();
    }
    if (!equal(AT::atan(AT::tan(three)), three)) {
      Kokkos::printf("AT::atan(tan(3)) != 3\n");
      FAILURE();
    }

    // Call the base class' implementation.  Every subclass'
    // implementation of testHostImpl() should (must) do this, in
    // order to include the parent class' tests.  In the case of this
    // particular class, the base class' implementation doesn't do
    // anything, but that's OK.
    const int parentSuccess = base_type::testHostImpl(out);
    success                 = success && parentSuccess;

    return success;
  }
};

/// \class ArithTraitsTesterComplexBase
/// \brief Execute Kokkos::ArithTraits tests relevant to
///   complex numbers (whether or not \c ScalarType is itself a
///   complex-valued type).
///
/// \tparam ScalarType The template parameter of ArithTraits
/// \tparam DeviceType The Kokkos device type over which to execute tests
/// \tparam is_complex Whether \c ScalarType is a complex-valued type
///
/// Some tests will be executed whether or not <tt>ScalarType</tt> is
/// complex, but the specific tests that are run will depend on
/// <tt>ScalarType</tt>.
template <class ScalarType, class DeviceType, const int is_complex = Kokkos::ArithTraits<ScalarType>::is_complex>
class ArithTraitsTesterComplexBase : public ArithTraitsTesterTranscendentalBase<ScalarType, DeviceType> {
 private:
  //! The base class of this class.
  typedef ArithTraitsTesterTranscendentalBase<ScalarType, DeviceType> base_type;

 public:
  typedef DeviceType execution_space;
  typedef typename execution_space::size_type size_type;
  //! Type of the result of the reduction.
  typedef int value_type;

  /// \brief The "parallel for" part of the reduction.
  ///
  /// See comments of ArithTraitsTesterBase's operator().
  KOKKOS_INLINE_FUNCTION void operator()(size_type iwork, value_type& dst) const;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterComplexBase();

 protected:
  // The host hook gets implemented in the complex-arithmetic
  // specialization of this class.
  virtual int testHostImpl(std::ostream& out) const;
};

//
// Specialization of ArithTraitsTesterComplexBase for real T.
//
template <class ScalarType, class DeviceType>
class ArithTraitsTesterComplexBase<ScalarType, DeviceType, 0>
    : public ArithTraitsTesterTranscendentalBase<ScalarType, DeviceType> {
 private:
  //! The base class of this class.
  typedef ArithTraitsTesterTranscendentalBase<ScalarType, DeviceType> base_type;

 public:
  typedef typename DeviceType::execution_space execution_space;
  typedef typename execution_space::size_type size_type;
  //! Type of the result of the reduction.
  typedef int value_type;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterComplexBase() {}

  KOKKOS_INLINE_FUNCTION void operator()(size_type iwork, value_type& dst) const {
    TRACE();
    typedef Kokkos::ArithTraits<ScalarType> AT;
    (void)iwork;  // forestall compiler warning for unused variable
    int success = 1;

    // Apparently, std::numeric_limits<ScalarType>::is_signed is 1
    // only for real numbers.
#if defined(KOKKOS_HALF_T_IS_FLOAT)
    if (std::is_same<ScalarType, Kokkos::Experimental::half_t>::value) {
      if (AT::is_signed != 0x1) FAILURE();
    } else
#else
    {
      if (AT::is_signed != std::numeric_limits<ScalarType>::is_signed) {
        Kokkos::printf(
            "AT::is_signed = 0x%x, std::numeric_limits<ScalarType>::is_signed "
            "= 0x%x\n",
            AT::is_signed, std::numeric_limits<ScalarType>::is_signed);
        FAILURE();
      }
    }
#endif  // KOKKOS_HALF_T_IS_FLOAT

      if (AT::is_complex) {
        FAILURE();
      }

    // Call the base class' implementation.  Every subclass'
    // implementation of operator() must do this, in order to include
    // the parent class' tests.
    int baseResult = 1;
    base_type::operator()(iwork, baseResult);
    success = success && baseResult;

    dst = dst && success;
  }

 protected:
  virtual int testHostImpl(std::ostream& out) const {
    using std::endl;
    typedef Kokkos::ArithTraits<ScalarType> AT;

    int success = 1;
    // Apparently, std::numeric_limits<ScalarType>::is_signed is 1 only for real
    // numbers.
    if (AT::is_signed != std::numeric_limits<ScalarType>::is_signed) {
      out << "ArithTraits<T>::is_signed != "
             "std::numeric_limits<ScalarType>::is_signed"
          << endl;
      FAILURE();
    }
    if (AT::is_complex) {
      out << "ArithTraits<T>::is_complex is wrong" << endl;
      FAILURE();
    }
    // Call the base class' implementation.  Every subclass'
    // implementation of testHostImpl() should (must) do this, in
    // order to include the parent class' tests.  In the case of this
    // particular class, the base class' implementation doesn't do
    // anything, but that's OK.
    const int parentSuccess = base_type::testHostImpl(out);
    success                 = success && parentSuccess;

    return success;
  }
};

// Specialization for complex T.
template <class ScalarType, class DeviceType>
class ArithTraitsTesterComplexBase<ScalarType, DeviceType, 1>
    : public ArithTraitsTesterTranscendentalBase<ScalarType, DeviceType> {
 private:
  //! The base class of this class.
  typedef ArithTraitsTesterTranscendentalBase<ScalarType, DeviceType> base_type;

 public:
  typedef typename DeviceType::execution_space execution_space;
  typedef typename execution_space::size_type size_type;
  //! Type of the result of the reduction.
  typedef int value_type;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterComplexBase() {}

  KOKKOS_INLINE_FUNCTION void operator()(size_type iwork, value_type& dst) const {
    TRACE();
    typedef Kokkos::ArithTraits<ScalarType> AT;
    (void)iwork;  // forestall compiler warning for unused variable
    int success = 1;

    if (!AT::is_complex) {
      FAILURE();
    }
    typedef typename AT::mag_type mag_type;
    const mag_type one = Kokkos::ArithTraits<mag_type>::one();

    // This presumes that ScalarType, being a complex number, has a
    // constructor which takes two mag_type arguments.
    const ScalarType oneMinusOne(one, -one);
    const ScalarType onePlusOne(one, one);

    // Test conjugation.
    if (AT::conj(oneMinusOne) != onePlusOne || AT::conj(onePlusOne) != oneMinusOne) {
      FAILURE();
    }

    // Call the base class' implementation.  Every subclass'
    // implementation of operator() must do this, in order to include
    // the parent class' tests.
    int baseResult = 1;
    base_type::operator()(iwork, baseResult);
    success = success && baseResult;

    dst = dst && success;
  }

 protected:
  virtual int testHostImpl(std::ostream& out) const {
    using std::endl;
    typedef Kokkos::ArithTraits<ScalarType> AT;
    int success = 1;

    if (!AT::is_complex) {
      out << "ArithTraits<T>::is_complex is wrong" << endl;
      FAILURE();
    }
    typedef typename AT::mag_type mag_type;
    const mag_type one = Kokkos::ArithTraits<mag_type>::one();

    // This presumes that ScalarType, being a complex number, has a
    // constructor which takes two mag_type arguments.
    const ScalarType oneMinusOne(one, -one);
    const ScalarType onePlusOne(one, one);

    // Test conjugation.
    if (AT::conj(oneMinusOne) != onePlusOne) {
      out << "AT::conj ((1, -1)) != (1, 1)" << endl;
      FAILURE();
    }
    if (AT::conj(onePlusOne) != oneMinusOne) {
      out << "AT::conj ((1, 1)) != (1, -1)" << endl;
      FAILURE();
    }
    // Call the base class' implementation.  Every subclass'
    // implementation of testHostImpl() should (must) do this, in
    // order to include the parent class' tests.  In the case of this
    // particular class, the base class' implementation doesn't do
    // anything, but that's OK.
    const int parentSuccess = base_type::testHostImpl(out);
    success                 = success && parentSuccess;

    return success;
  }
};

/// \class ArithTraitsTesterFloatingPointBase
/// \brief Kokkos reduction functor for testing those attributes of
///   ArithTraits suitable for floating-point types.
/// \tparam ScalarType A type suitable for execution on the parallel
///   device \c DeviceType.
/// \tparam DeviceType A Kokkos parallel device type.
///
/// Kokkos reduction operator for testing those attributes of
/// Kokkos::ArithTraits relevant to floating-point types.
///
/// The tests include those suitable for execution on the parallel
/// device (operator()) and those suitable for execution on the host
/// (testHost()).  The device-based test is a reduction over redundant
/// executions of the test.  All redundant executions must return
/// '1' (passed).
template <class ScalarType, class DeviceType, const int is_exact = Kokkos::ArithTraits<ScalarType>::is_exact>
class ArithTraitsTesterFloatingPointBase
    : public ArithTraitsTesterComplexBase<ScalarType, DeviceType, Kokkos::ArithTraits<ScalarType>::is_complex> {
 private:
  //! The base class of this class.
  typedef ArithTraitsTesterComplexBase<ScalarType, DeviceType, Kokkos::ArithTraits<ScalarType>::is_complex> base_type;

 public:
  typedef DeviceType execution_space;
  typedef typename execution_space::size_type size_type;
  //! Type of the result of the reduction.
  typedef int value_type;

  /// \brief The "parallel for" part of the reduction.
  ///
  /// See comments of ArithTraitsTesterBase's operator().
  KOKKOS_INLINE_FUNCTION void operator()(size_type iwork, value_type& dst) const;

 protected:
  virtual int testHostImpl(std::ostream& out) const;
};

//
// Specialization for is_exact = 0 (i.e., ScalarType is a
// floating-point type).
//
template <class ScalarType, class DeviceType>
class ArithTraitsTesterFloatingPointBase<ScalarType, DeviceType, 0>
    : public ArithTraitsTesterComplexBase<ScalarType, DeviceType, Kokkos::ArithTraits<ScalarType>::is_complex> {
 private:
  //! The base class of this class.
  typedef ArithTraitsTesterComplexBase<ScalarType, DeviceType, Kokkos::ArithTraits<ScalarType>::is_complex> base_type;

 public:
  typedef typename DeviceType::execution_space execution_space;
  typedef typename execution_space::size_type size_type;
  //! Type of the result of the reduction.
  typedef int value_type;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterFloatingPointBase() {}

  KOKKOS_INLINE_FUNCTION void operator()(size_type iwork, value_type& dst) const {
    TRACE();
    typedef Kokkos::ArithTraits<ScalarType> AT;
    (void)iwork;  // forestall compiler warning for unused variable
    int success = 1;

    if (AT::is_exact) {
      Kokkos::printf("AT::is_exact is 1\n");
      FAILURE();
    }

// Compiler intrinsic casts from nan of type half_t / bhalf_t to nan
// of type float in CUDA, SYCL and HIP do not work yet.
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ENABLE_HIP)
    namespace KE = Kokkos::Experimental;
    if constexpr (!std::is_same<ScalarType, KE::half_t>::value && !std::is_same<ScalarType, KE::bhalf_t>::value) {
#else
    {
#endif  // KOKKOS_ENABLE_CUDA || KOKKOS_ENABLE_SYCL || KOKKOS_ENABLE_HIP
      if (!AT::isNan(AT::nan())) {
        Kokkos::printf("NaN is not NaN\n");
        FAILURE();
      }
    }

    const ScalarType zero = AT::zero();
    const ScalarType one  = AT::one();

    if (AT::isInf(zero)) {
      Kokkos::printf("0 is Inf\n");
      FAILURE();
    }
    if (AT::isInf(one)) {
      Kokkos::printf("1 is Inf\n");
      FAILURE();
    }
#if defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ENABLE_HIP)  // FIXME_SYCL, FIXME_HIP
    if constexpr (!std::is_same_v<ScalarType, Kokkos::Experimental::half_t>) {
      if (AT::isNan(zero)) {
        Kokkos::printf("0 is NaN\n");
        FAILURE();
      }
      if (AT::isNan(one)) {
        Kokkos::printf("1 is NaN\n");
        FAILURE();
      }
    }
#else
    if (AT::isNan(zero)) {
      Kokkos::printf("0 is NaN\n");
      FAILURE();
    }
    if (AT::isNan(one)) {
      Kokkos::printf("1 is NaN\n");
      FAILURE();
    }
#endif

    // Call the base class' implementation.  Every subclass'
    // implementation of operator() must do this, in order to include
    // the parent class' tests.
    int baseResult = 1;
    base_type::operator()(iwork, baseResult);
    success = success && baseResult;

    dst = dst && success;
  }

 protected:
  virtual int testHostImpl(std::ostream& out) const {
    typedef Kokkos::ArithTraits<ScalarType> AT;
    using std::endl;
    int success = 1;

    if (AT::is_exact) {
      out << "AT::is_exact is wrong" << endl;
      FAILURE();
    }

    // if (std::numeric_limits<ScalarType>::is_iec559) {
    // success = success && AT::isInf (AT::inf ());
#if defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ENABLE_HIP)
    if constexpr (!std::is_same_v<ScalarType, Kokkos::Experimental::half_t>) {
      if (!AT::isNan(AT::nan())) {
        out << "isNan or nan failed" << endl;
        FAILURE();
      }
    }
#else
    if (!AT::isNan(AT::nan())) {
      out << "isNan or nan failed" << endl;
      FAILURE();
    }
#endif
    //}

    const ScalarType zero = AT::zero();
    const ScalarType one  = AT::one();

    if (AT::isInf(zero)) {
      out << "isInf(zero) is 1" << endl;
      FAILURE();
    }
    if (AT::isInf(one)) {
      out << "isInf(one) is 1" << endl;
      FAILURE();
    }
#if defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ENABLE_HIP)
    if constexpr (!std::is_same_v<ScalarType, Kokkos::Experimental::half_t>) {
      if (AT::isNan(zero)) {
        out << "isNan(zero) is 1" << endl;
        FAILURE();
      }
      if (AT::isNan(one)) {
        out << "isNan(one) is 1" << endl;
        FAILURE();
      }
    }
#else
    if (AT::isNan(zero)) {
      out << "isNan(zero) is 1" << endl;
      FAILURE();
    }
    if (AT::isNan(one)) {
      out << "isNan(one) is 1" << endl;
      FAILURE();
    }
#endif

    // Call the base class' implementation.  Every subclass'
    // implementation of testHostImpl() should (must) do this, in
    // order to include the parent class' tests.
    const int parentSuccess = base_type::testHostImpl(out);
    success                 = success && parentSuccess;

    return success;
  }
};

//
// Specialization for is_exact = 1 (i.e., ScalarType is <i>not</i>
// a floating-point type).
//
template <class ScalarType, class DeviceType>
class ArithTraitsTesterFloatingPointBase<ScalarType, DeviceType, 1>
    : public ArithTraitsTesterComplexBase<ScalarType, DeviceType, Kokkos::ArithTraits<ScalarType>::is_complex> {
 private:
  //! The base class of this class.
  typedef ArithTraitsTesterComplexBase<ScalarType, DeviceType, Kokkos::ArithTraits<ScalarType>::is_complex> base_type;

 public:
  typedef typename DeviceType::execution_space execution_space;
  typedef typename execution_space::size_type size_type;
  //! Type of the result of the reduction.
  typedef int value_type;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterFloatingPointBase() {}

  KOKKOS_INLINE_FUNCTION void operator()(size_type iwork, value_type& dst) const {
    TRACE();
    typedef Kokkos::ArithTraits<ScalarType> AT;
    (void)iwork;  // forestall compiler warning for unused variable
    int success = 1;

    if (!AT::is_exact) {
      Kokkos::printf("! AT:is_exact\n");
      FAILURE();
    }

    // Call the base class' implementation.  Every subclass'
    // implementation of operator() must do this, in order to include
    // the parent class' tests.
    int baseResult = 1;
    base_type::operator()(iwork, baseResult);
    success = success && baseResult;

    dst = dst && success;
  }

 protected:
  virtual int testHostImpl(std::ostream& out) const {
    typedef Kokkos::ArithTraits<ScalarType> AT;
    using std::endl;
    int success = 1;

    if (!AT::is_exact) {
      out << "AT::is_exact is wrong" << endl;
      FAILURE();
    }
    // Call the base class' implementation.  Every subclass'
    // implementation of testHostImpl() should (must) do this, in
    // order to include the parent class' tests.
    const int parentSuccess = base_type::testHostImpl(out);
    success                 = success && parentSuccess;

    return success;
  }
};

/// \class ArithTraitsTester
/// \brief Tests for Kokkos::ArithTraits
/// \tparam ScalarType Any type for which Kokkos::ArithTraits
///   has a specialization, and which can be executed on the parallel
///   device.
/// \tparam DeviceType A Kokkos parallel device type.
///
/// This class works by inheriting "hooks" from ArithTraitsTesterBase
/// and the chain of subclasses in between.  ArithTraitsTesterBase
/// provides basic tests that work for all <tt>ScalarType</tt>, and
/// the subclasses provide additional tests relevant to things like
/// complex-valued types or floating-point types.  The hooks for
/// device functions do <i>not</i> use run-time polymorphism, since
/// this does not always work with CUDA device functions.  The hooks
/// for host functions <i>do</i> use run-time polymorphism.
///
/// This class (through its base class) provides a Kokkos reduction
/// operator for testing Kokkos::ArithTraits.  This test
/// works for any type <tt>ScalarType</tt> for which
/// Kokkos::ArithTraits has a specialization, and which can
/// be executed on the parallel device.
///
/// The tests include those suitable for execution on the parallel
/// device (operator()) and those suitable for execution on the host
/// (testHost()).  The device-based test is a reduction over redundant
/// executions of the test.  All redundant executions must return
/// '1' (passed).
template <class ScalarType, class DeviceType>
class ArithTraitsTester : public ArithTraitsTesterFloatingPointBase<ScalarType, DeviceType> {
 public:
  typedef typename DeviceType::execution_space execution_space;
  typedef typename execution_space::size_type size_type;
  //! Type of the result of the reduction.
  typedef int value_type;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTester() {}
};

/// \brief Run the Kokkos::ArithTraits tests on the parallel device.
/// \tparam ScalarType Any type for which Kokkos::ArithTraits
///   has a specialization, and which can be executed on the parallel
///   device.
/// \tparam DeviceType A Kokkos parallel device type.
///
/// This function must be called on the host, but it executes on the
/// device with (redundant) parallelism.
///
/// \return \c 1 if all redundant executions pass, else \c 0.
template <class ScalarType, class DeviceType>
int testArithTraitsOnDevice(std::ostream& out, const int verbose) {
  using std::endl;
  typedef ArithTraitsTester<ScalarType, DeviceType> functor_type;
  int success = 1;  // output argument of parallel_reduce
  Kokkos::parallel_reduce("KokkosKernels::Common::Test::ArithTraitsOnDevice", 1, functor_type(), success);
  if (success) {
    if (verbose) out << Kokkos::ArithTraits<ScalarType>::name() << " passed" << endl;
  } else {
    out << Kokkos::ArithTraits<ScalarType>::name() << " FAILED" << endl;
  }
  return success;
}

/// \brief Run the Kokkos::ArithTraits tests on the host.
/// \tparam ScalarType Any type for which Kokkos::ArithTraits
///   has a specialization.
/// \tparam DeviceType A Kokkos parallel device type.
///
/// This function must be called on the host, and executes on the host.
///
/// \return \c 1 if all tests pass, else \c 0.
template <class ScalarType, class DeviceType>
int testArithTraitsOnHost(std::ostream& out, const int verbose) {
  using std::endl;
  ArithTraitsTester<ScalarType, DeviceType> f;
  const int localSuccess = f.testHost(out);

  if (localSuccess) {
    if (verbose) out << Kokkos::ArithTraits<ScalarType>::name() << " passed" << endl;
  } else {
    out << Kokkos::ArithTraits<ScalarType>::name() << " FAILED" << endl;
  }
  return localSuccess;
}

/// \brief Run the Kokkos::ArithTraits tests for all (valid)
///   scalar types, on the given parallel device.
/// \tparam DeviceType A Kokkos parallel device type.
///
/// This is the "outward-facing" function meant to be called by
/// main().  This function must be called on the host, but it executes
/// on the device with (redundant) parallelism.
///
/// \return \c 1 if all tests pass, else \c 0.
template <class DeviceType>
int runAllArithTraitsDeviceTests(std::ostream& out, const int verbose) {
  int success    = 1;
  int curSuccess = 1;
  //
  // Built-in char(acter) types
  //

  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<char, DeviceType>(out, verbose);
  // Interestingly enough, char and int8_t are different types, but
  // signed char and int8_t are the same (on my system).
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<signed char, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<unsigned char, DeviceType>(out, verbose);

  //
  // Built-in integer types
  //

  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<short, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<unsigned short, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<int8_t, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<uint8_t, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<int16_t, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<uint16_t, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<int32_t, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<uint32_t, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<int, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<unsigned int, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<int64_t, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<uint64_t, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<long, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<unsigned long, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<long long, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<unsigned long long, DeviceType>(out, verbose);

  //
  // Built-in real floating-point types
  //

#if defined(KOKKOS_HALF_T_IS_FLOAT)
  TRACE();
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<Kokkos::Experimental::half_t, DeviceType>(out, verbose);
#endif  // KOKKOS_HALF_T_IS_FLOAT
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<float, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<double, DeviceType>(out, verbose);

  //
  // Kokkos' complex floating-point types
  //

  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<Kokkos::complex<float>, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnDevice<Kokkos::complex<double>, DeviceType>(out, verbose);

  return success && curSuccess;
}

/// \brief Run the Kokkos::ArithTraits tests for all scalar
///   types, on the host.
/// \tparam DeviceType A Kokkos parallel device type.
///
/// This is the "outward-facing" function meant to be called by
/// main().  This function must be called on the host, and executes on
/// the host.
///
/// \return \c 1 if all tests pass, else \c 0.
template <class DeviceType>
int runAllArithTraitsHostTests(std::ostream& out, const int verbose) {
  int success    = 1;
  int curSuccess = 1;

  //
  // Built-in char(acter) types
  //

  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<char, DeviceType>(out, verbose);
  // Interestingly enough, char and int8_t are different types, but
  // signed char and int8_t are the same (on my system).
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<signed char, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<unsigned char, DeviceType>(out, verbose);

  //
  // Built-in integer types
  //

  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<short, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<unsigned short, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<int8_t, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<uint8_t, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<int16_t, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<uint16_t, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<int32_t, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<uint32_t, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<int, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<unsigned int, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<int64_t, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<uint64_t, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<long, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<unsigned long, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<long long, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<unsigned long long, DeviceType>(out, verbose);

  //
  // Built-in real and complex floating-point types
  //

  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<float, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<double, DeviceType>(out, verbose);
#if !defined(KOKKOS_ENABLE_CUDA) && !defined(KOKKOS_ENABLE_HIP) && !defined(KOKKOS_ENABLE_SYCL)
  // This would spill tons of warnings about host device stuff otherwise
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<long double, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<std::complex<float>, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<std::complex<double>, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<std::complex<long double>, DeviceType>(out, verbose);
#endif
  //
  // Kokkos' complex floating-point types
  //

#if defined(KOKKOS_HALF_T_IS_FLOAT)
  success = success && curSuccess;
  TRACE();
  curSuccess = testArithTraitsOnHost<Kokkos::Experimental::half_t, DeviceType>(out, verbose);
#endif  // KOKKOS_HALF_T_IS_FLOAT
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<Kokkos::complex<float>, DeviceType>(out, verbose);
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<Kokkos::complex<double>, DeviceType>(out, verbose);
  // success = success && curSuccess; curSuccess =
  // testArithTraitsOnHost<Kokkos::complex<long double>, DeviceType> (out,
  // verbose);

#if defined(KOKKOS_ENABLE_LIBQUADMATH)
  success    = success && curSuccess;
  curSuccess = testArithTraitsOnHost<__float128, DeviceType>(out, verbose);
#endif
  return success && curSuccess;
}

template <typename device>
void test_ArithTraits() {
  using std::endl;

  class NullBuffer : public std::streambuf {
   public:
    int overflow(int c) { return c; }
  };
  NullBuffer null_buffer;
  std::ostream& out = std::cerr;
  // std::ostream out(&null_buffer);

  bool success = true;
  success      = runAllArithTraitsDeviceTests<device>(out, 0);
  EXPECT_TRUE(success);
  success = runAllArithTraitsHostTests<device>(out, 0);
  EXPECT_TRUE(success);
}
TEST_F(TestCategory, common_ArithTraits) { test_ArithTraits<TestDevice>(); }

#endif  // KOKKOS_ARITHTRAITSTEST_HPP
