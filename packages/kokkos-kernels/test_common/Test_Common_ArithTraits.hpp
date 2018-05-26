/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
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

/// \file ArithTraitsTest.hpp
/// \brief Templated test for Kokkos::Details::ArithTraits
///
/// This header file is an implementation detail of the tests for
/// Kokkos::Details::ArithTraits.  Users must not rely on it existing,
/// or on its contents.  This header file should <i>not</i> be
/// installed with Kokkos' other header files.
///
/// On the other hand, this header file does give examples of how to
/// use Kokkos::Details::ArithTraits, so it may be useful for users to
/// read it.

#ifndef KOKKOS_ARITHTRAITSTEST_HPP
#define KOKKOS_ARITHTRAITSTEST_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_ArithTraits.hpp"
#include <limits> // std::numeric_limits
#include <typeinfo> // typeid (T)
#include <cstdio>


namespace {
  // Whether Kokkos::Details::ArithTraits<ScalarType> implements
  // transcendental functions.  These include sqrt, pow, log, and
  // log10.
  template<class ScalarType>
  struct HasTranscendentals {
    static const bool value = false;
  };

  template<>
  struct HasTranscendentals<float> {
    static const bool value = true;
  };

  template<>
  struct HasTranscendentals<double> {
    static const bool value = true;
  };

  template<>
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

} // namespace (anonymous)




/// \class ArithTraitsTesterBase
/// \brief Base class providing tests for Kokkos::Details::ArithTraits
/// \tparam ScalarType Any type for which Kokkos::Details::ArithTraits
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
/// Kokkos::Details::ArithTraits.  This test works for any type
/// <tt>ScalarType</tt> for which Kokkos::Details::ArithTraits has a
/// specialization, and which can be executed on the parallel device.
///
/// The tests include those suitable for execution on the parallel
/// device (operator()) and those suitable for execution on the host
/// (testHost()).  The device-based test is a reduction over redundant
/// executions of the test.  All redundant executions must return
/// '1' (passed).
template<class ScalarType, class DeviceType>
class ArithTraitsTesterBase {
public:
  typedef DeviceType execution_space;
  typedef typename execution_space::size_type size_type;
  //! Type of the result of the reduction.
  typedef int value_type;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterBase () {}

  /// \brief Set the initial value (\c 1) of the reduction.
  ///
  /// Subclasses need not and must not override this method.
  KOKKOS_INLINE_FUNCTION void init ( value_type& dst) const {
    dst = 1;
  }

  /// \brief Combine two intermediate reduction results into \c dst.
  ///
  /// Subclasses need not and must not override this method.
  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& dst,
        const volatile value_type& src) const
  {
    dst = dst && src;
    //dst = 1;
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
  KOKKOS_INLINE_FUNCTION void
  operator () (size_type iwork, value_type& dst) const
  {
    typedef Kokkos::Details::ArithTraits<ScalarType> AT;
    (void) iwork; // not using this argument
    int success = 1;

    // Make sure that the typedef exists.
    typedef typename AT::mag_type mag_type;

    // mfh 14 Feb 2014: In order to avoid a warning for an unused by
    // declared typedef, we declare an instance of mag_type, and mark
    // it with "(void)" to prevent a warning for the unused variable.
    {
      mag_type thing;
      (void) thing;
    }

    // ArithTraits should not even compile if it's not specialized for
    // T, but we check for this int constant for compatibility with
    // std::numeric_limits.
    if (! AT::is_specialized) {
      printf ("! AT::is_specialized\n");
      success = 0;
    }

    // It's OK to refer to std::numeric_limits constants in a device
    // function, just not to its class methods (which are not marked
    // as device functions).
    if (AT::is_integer != std::numeric_limits<ScalarType>::is_integer) {
      printf ("AT::is_integer not same as numeric_limits\n");
      success = 0;
    }
    if (AT::is_exact != std::numeric_limits<ScalarType>::is_exact) {
      printf ("AT::is_exact not same as numeric_limits\n");
      success = 0;
    }

    const ScalarType zero = AT::zero ();
    const ScalarType one = AT::one ();

    // Test properties of the arithmetic and multiplicative identities.
    if (zero + zero != zero) {
      printf ("0 + 0 != 0\n");
      success = 0;
    }
    if (zero + one != one) {
      printf ("0 + 1 != 1\n");
      success = 0;
    }
    if (one - one != zero) {
      printf ("1 - 1 != 0\n");
      success = 0;
    }
    // This is technically 1 even of Z_2, since in that field, one
    // is its own inverse (so -one == one).
    if ((one + one) - one != one) {
      printf ("(1 + 1) - 1 != 1\n");
      success = 0;
    }

    if (AT::abs (zero) != zero) {
      printf ("AT::abs(0) != 0\n");
      success = 0;
    }
    if (AT::abs (one) != one) {
      printf ("AT::abs(1) != 1\n");
      success = 0;
    }
    if (AT::is_signed && AT::abs (-one) != one) {
      printf ("AT::is_signed and AT::abs(-1) != 1\n");
      success = 0;
    }
    // Need enable_if to test whether T can be compared using <=.
    // However, mag_type should always be comparable using <=.
    //
    // These are very mild ordering properties.
    // They should work even for a set only containing zero.
    if (AT::abs (zero) > AT::abs (AT::max ())) {
      printf ("AT::abs(0) > AT::abs (AT::max ())\n");
      success = 0;
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
  int testHostImpl (std::ostream& out) const {
    return 1; // there are no tests, so trivially, all the tests pass
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
  int testHost (std::ostream& out) const {
    typedef Kokkos::Details::ArithTraits<ScalarType> AT;
    using std::endl;
    int success = 1;

    // Make sure that the typedef exists.
    typedef typename AT::mag_type mag_type;

    // mfh 14 Feb 2014: In order to avoid a warning for an unused by
    // declared typedef, we declare an instance of mag_type, and mark
    // it with "(void)" to prevent a warning for the unused variable.
    {
      mag_type thing;
      (void) thing;
    }

    // ArithTraits should not even compile if it's not specialized for
    // T, but we check for this int constant for compatibility with
    // std::numeric_limits.
    if (! AT::is_specialized) {
      out << "ArithTraits is not specialized for T" << endl;
      success = 0;
    }

    if (AT::is_integer != std::numeric_limits<ScalarType>::is_integer) {
      out << "AT::is_integer != std::numeric_limits<ScalarType>::is_integer" << endl;
      success = 0;
    }

    if (AT::is_exact != std::numeric_limits<ScalarType>::is_exact) {
      out << "AT::is_exact != std::numeric_limits<ScalarType>::is_exact" << endl;
      success = 0;
    }

    const ScalarType zero = AT::zero ();
    const ScalarType one = AT::one ();
    // Test properties of the arithmetic and multiplicative identities.

    if (zero + zero != zero) {
      out << "zero + zero != zero" << endl;
      success = 0;
    }
    if (zero + one != one) {
      out << "zero + one != one" << endl;
      success = 0;
    }
    if (one - one != zero) {
      out << "one - one != zero" << endl;
      success = 0;
    }
    // This is technically 1 even of Z_2, since in that field, one
    // is its own inverse (so -one == one).
    if ((one + one) - one != one) {
      out << "(one + one) - one != one" << endl;
      success = 0;
    }

    if (AT::abs (zero) != zero) {
      out << "AT::abs (zero) != zero" << endl;
      success = 0;
    }
    if (AT::abs (one) != one) {
      out << "AT::abs (one) != one" << endl;
      success = 0;
    }
    if (AT::is_signed) {
      if (AT::abs (-one) != one) {
        out << "AT::abs (-one) != one" << endl;
        success = 0;
      }
    }
    // Need enable_if to test whether T can be compared using <=.
    // However, mag_type should always be comparable using <=.
    //
    // // These are very mild ordering properties.
    // // They should work even for a set only containing zero.
    if (AT::abs (zero) > AT::abs (AT::max ())) {
      out << "AT::abs (zero) > AT::abs (AT::max ())" << endl;
      success = 0;
    }

    // Run the parent class' remaining tests, if any.
    const int parentSuccess = testHostImpl (out);
    success = success && parentSuccess;

    return success;
  }
};


/// \class ArithTraitsTesterTranscendentalBase
/// \brief Base class of ArithTraitsTester that exercises
///   transcendental functions, if and only if ArithTraits<ScalarType>
///   implements them.
/// \tparam ScalarType Any type for which Kokkos::Details::ArithTraits
///   implements transcendental functions, along with the requirements
///   imposed by ArithTraitsTesterBase.
/// \tparam DeviceType A Kokkos parallel device type.
/// \tparam hasTranscendentals Whether ArithTraits<ScalarType>
///   implements transcendental functions.
///
/// Some tests will be executed whether or not ArithTraits<ScalarType>
/// implements transcendental functions, but the specific tests that
/// are run will depend on \c ScalarType.
template<class ScalarType,
         class DeviceType,
         const int has_transcendentals =
         (HasTranscendentals<ScalarType>::value ? 1 : 0) >
class ArithTraitsTesterTranscendentalBase :
  public ArithTraitsTesterBase<ScalarType, DeviceType> {
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
  KOKKOS_INLINE_FUNCTION void
  operator () (size_type iwork, value_type& dst) const;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterTranscendentalBase ();

protected:
  // The host hook gets implemented in the "transcendental functions
  // are implemented" specialization of this class.
  virtual int testHostImpl (std::ostream& out) const;
};


//
// Specialization of ArithTraitsTesterTranscendentalBase when
// ArithTraits<ScalarType> does NOT implement transcendentals.
//
template<class ScalarType,
         class DeviceType>
class ArithTraitsTesterTranscendentalBase<ScalarType, DeviceType, 0> :
  public ArithTraitsTesterBase<ScalarType, DeviceType> {
private:
  //! The base class of this class.
  typedef ArithTraitsTesterBase<ScalarType, DeviceType> base_type;

public:
  typedef DeviceType execution_space;
  typedef typename execution_space::size_type size_type;
  //! Type of the result of the reduction.
  typedef int value_type;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterTranscendentalBase () {}

  KOKKOS_INLINE_FUNCTION void
  operator () (size_type iwork, value_type& dst) const {
    //typedef Kokkos::Details::ArithTraits<ScalarType> AT;
    (void) iwork; // forestall compiler warning for unused variable
    int success = 1;

    if (HasTranscendentals<ScalarType>::value) {
      success = 0;
    }

    // Call the base class' implementation.  Every subclass'
    // implementation of operator() must do this, in order to include
    // the parent class' tests.
    int baseResult = 1;
    base_type::operator () (iwork, baseResult);
    success = success && baseResult;

    dst = dst && success;
  }

protected:
  virtual int testHostImpl (std::ostream& out) const {
    using std::endl;
    //typedef Kokkos::Details::ArithTraits<ScalarType> AT;
    int success = 1;

    if (HasTranscendentals<ScalarType>::value) {
      out << "HasTranscendentals<T>::value is true" << endl;
      success = 0;
    }

    // Call the base class' implementation.  Every subclass'
    // implementation of testHostImpl() should (must) do this, in
    // order to include the parent class' tests.  In the case of this
    // particular class, the base class' implementation doesn't do
    // anything, but that's OK.
    const int parentSuccess = base_type::testHostImpl (out);
    success = success && parentSuccess;

    return success;
  }
};


//
// Specialization of ArithTraitsTesterTranscendentalBase when
// ArithTraits<ScalarType> DOES implement transcendentals.
//
template<class ScalarType,
         class DeviceType>
class ArithTraitsTesterTranscendentalBase<ScalarType, DeviceType, 1> :
  public ArithTraitsTesterBase<ScalarType, DeviceType> {
private:
  //! The base class of this class.
  typedef ArithTraitsTesterBase<ScalarType, DeviceType> base_type;

  KOKKOS_INLINE_FUNCTION
  bool equal(const ScalarType& a, const ScalarType& b) const {
    if(b!=Kokkos::Details::ArithTraits<ScalarType>::zero()) {
      if(a>b)
        return (a-b)/b < 2 * Kokkos::Details::ArithTraits<ScalarType>::epsilon();
      else
        return (b-a)/b < 2 * Kokkos::Details::ArithTraits<ScalarType>::epsilon();
    } else {
      if(a>b)
        return (a-b) < 2 * Kokkos::Details::ArithTraits<ScalarType>::epsilon();
      else
        return (b-a) < 2 * Kokkos::Details::ArithTraits<ScalarType>::epsilon();
    }
  }

public:
  typedef DeviceType execution_space;
  typedef typename execution_space::size_type size_type;
  //! Type of the result of the reduction.
  typedef int value_type;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterTranscendentalBase () {}

  KOKKOS_INLINE_FUNCTION void
  operator () (size_type iwork, value_type& dst) const {
    typedef Kokkos::Details::ArithTraits<ScalarType> AT;
    (void) iwork; // forestall compiler warning for unused variable
    int success = 1;

    if (! HasTranscendentals<ScalarType>::value) {
      success = 0;
    }

    const ScalarType zero = AT::zero ();
    const ScalarType one = AT::one ();
    const ScalarType two = one + one;
    const ScalarType three = one + one + one;
    const ScalarType four = two * two;
    const ScalarType five = four + one;
    const ScalarType six = three * two;
    const ScalarType seven = four + three;
    const ScalarType eight = four * two;
    const ScalarType nine = eight + one;
    const ScalarType eleven = five + six;
    const ScalarType twentySeven = nine * three;
    const ScalarType thirtySix = six * six;
    const ScalarType fortyTwo = six * seven;
    const ScalarType sixtyThree = eight * eight - one;
    const ScalarType sixtyFour = eight * eight;
    // max char value, for 8-bit char
    const ScalarType oneTwentySeven = sixtyFour + sixtyThree;

    ScalarType result;

    // This fails inexplicably for complex numbers on gcc 4.2.1 on Mac.
    if (! AT::is_complex) {
      result = AT::pow (two, three);
      if (!equal(result,eight)) {
        printf ("AT::pow(2,3) != 8\n");
        success = 0;
      }
    }
    if (!equal(AT::pow (three, zero) , one)) {
      printf ("AT::pow(3,0) != 1\n");
      success = 0;
    }
    if (!equal(AT::pow (three, one) , three)) {
      printf ("AT::pow(3,1) != 3\n");
      success = 0;
    }
    if (!equal(AT::pow (three, two) , nine)) {
      printf ("AT::pow(3,2) != 9\n");
      success = 0;
    }

    // This fails inexplicably for complex numbers on gcc 4.2.1 on Mac.
    if (! AT::is_complex) {
      result = AT::pow (three, three);
      if (!equal(result , twentySeven)) {
        printf ("AT::pow(3,3) != 27\n");
        success = 0;
      }
    }

    // These fail inexplicably for complex numbers on gcc 4.2.1 on Mac.
    if (AT::is_signed && ! AT::is_complex) {
      result = AT::pow (-three, one);
      if (!equal(result , -three)) {
        printf ("AT::pow(-3,1) != -3\n");
        success = 0;
      }
      result = AT::pow (-three, two);
      if (!equal(result , nine)) {
        printf ("AT::pow(-3,2) != 9\n");
        success = 0;
      }
      result = AT::pow (-three, three);
      if (!equal(result , -twentySeven)) {
        printf ("AT::pow(-3,3) != 27\n");
        success = 0;
      }
    }

    if (!equal(AT::sqrt (zero) , zero)) {
      printf ("AT::sqrt(0) != 0\n");
      success = 0;
    }
    if (!equal(AT::sqrt (one) , one)) {
      printf ("AT::sqrt(1) != 1\n");
      success = 0;
    }
    if (!equal(AT::sqrt (thirtySix) , six)) {
      printf ("AT::sqrt(36) != 6\n");
      success = 0;
    }
    if (!equal(AT::sqrt (sixtyFour) , eight)) {
      printf ("AT::sqrt(64) != 8\n");
      success = 0;
    }
    if (AT::is_integer) {
      if (!equal(AT::sqrt (fortyTwo) , six)) {
        printf ("AT:sqrt(42) != 6\n");
        success = 0;
      }
      if (!equal(AT::sqrt (oneTwentySeven) , eleven)) {
        printf ("AT::sqrt(127) != 11\n");
        success = 0;
      }
    }

    if (!equal(AT::cbrt (zero) , zero)) {
      printf ("AT::cbrt(0) != 0\n");
      success = 0;
    }
    if (!equal(AT::cbrt (one) , one)) {
      printf ("AT::cbrt(1) != 1\n");
      success = 0;
    }
    if (!equal(AT::cbrt (twentySeven) , three)) {
      printf ("AT::cbrt(27) != 3\n");
      success = 0;
    }
    if (!equal(AT::cbrt (sixtyFour) , four)) {
      printf ("AT::cbrt(64) != 4\n");
      success = 0;
    }
    if (AT::is_integer) {
      if (!equal(AT::cbrt (fortyTwo) , three)) {
        printf ("AT:cbrt(42) != 3\n");
        success = 0;
      }
      if (!equal(AT::cbrt (oneTwentySeven) , five)) {
        printf ("AT::cbrt(127) != 5\n");
        success = 0;
      }
    }

    if (!equal(AT::exp (zero) , one)) {
      printf ("AT::cbrt(0) != 1\n");
      success = 0;
    }
    if (AT::is_complex) {
      const ScalarType val = two; //(two.real(), two.real());
      if (!equal(AT::conj (AT::exp  (val)) , 
                 AT::exp  (AT::conj (val)))) {
        printf ("AT::conj(exp(complex(2,2))) != AT::exp(conj(complex(2,2)))\n");
        success = 0;
      }
    }
    if (!equal(AT::log (one) , zero)) {
      printf ("AT::log(1) != 0\n");
      success = 0;
    }
    if (!equal(AT::log10 (one) , zero)) {
      printf ("AT::log10(1) != 0\n");
      success = 0;
    }

    if (AT::is_complex) {
      ScalarType val = two; //(two, two);
      const auto val_sin = AT::sin (val);
      const auto val_cos = AT::cos (val);
      if (!equal(val_sin*val_sin + val_cos*val_cos , one)) {
        printf ("AT(complex):: sin(val)*sin(val) + cos(val)*cos(val) != 1\n");
        success = 0;
      } 
      if (!equal(val_sin/val_cos , AT::tan(val))) {
        printf ("AT(complex):: sin(val)/cos(val) != AT(real)::tan(val)\n");
        success = 0;
      } 
    } else {
      ScalarType val = three; 
      const auto val_sin = AT::sin (val);
      const auto val_cos = AT::cos (val);
      if (!equal(val_sin*val_sin + val_cos*val_cos , one)) {
        printf ("AT(real):: sin(val)*sin(val) + cos(a)*cos(a) != 1\n");
        success = 0;
      } 
      if (!equal(val_sin/val_cos , AT::tan(val))) {
        printf ("AT(real):: sin(val)/cos(val) != AT(real)::tan(val)\n");
        success = 0;
      } 
    }

    if (!equal(AT::asin (AT::sin (one)), one)) {
      printf ("AT::asin(sin(1)) != 1\n");
      success = 0;
    } 
    if (!equal(AT::acos (AT::cos (one)), one)) {
      printf ("AT::acos(cos(1)) != 1\n");
      success = 0;
    } 
    if (!equal(AT::atan (AT::tan (one)), one)) {
      printf ("AT::atan(tan(1)) != 1\n");
      success = 0;
    } 

    // Call the base class' implementation.  Every subclass'
    // implementation of operator() must do this, in order to include
    // the parent class' tests.
    int baseResult = 1;
    base_type::operator () (iwork, baseResult);
    success = success && baseResult;

    dst = dst && success;
  }

protected:
  virtual int testHostImpl (std::ostream& out) const {
    using std::endl;
    typedef Kokkos::Details::ArithTraits<ScalarType> AT;
    int success = 1;

    if (! HasTranscendentals<ScalarType>::value) {
      out << "HasTranscendentals<T>::value is false" << endl;
      success = 0;
    }

    const ScalarType zero = AT::zero ();
    const ScalarType one = AT::one ();
    const ScalarType two = one + one;
    const ScalarType three = one + one + one;
    const ScalarType four = two * two;
    const ScalarType five = four + one;
    const ScalarType six = three * two;
    const ScalarType seven = four + three;
    const ScalarType eight = four * two;
    const ScalarType nine = eight + one;
    const ScalarType eleven = five + six;
    const ScalarType twentySeven = nine * three;
    const ScalarType thirtySix = six * six;
    const ScalarType fortyTwo = six * seven;
    const ScalarType sixtyThree = eight * eight - one;
    const ScalarType sixtyFour = eight * eight;
    // max char value, for 8-bit char
    const ScalarType oneTwentySeven = sixtyFour + sixtyThree;

    ScalarType result;

    // This fails inexplicably for complex numbers on gcc 4.2.1 on Mac.
    if (! AT::is_complex) {
      result = AT::pow (two, three);
      if (result != eight) {
        out << "AT::pow (two, three) != eight" << endl;
        success = 0;
      }
    }
    if (AT::pow (three, zero) != one) {
      out << "AT::pow (three, zero) != one" << endl;
      success = 0;
    }
    if (AT::pow (three, one) != three) {
      out << "AT::pow (three, one) != three" << endl;
      success = 0;
    }
    if (AT::pow (three, two) != nine) {
      out << "AT::pow (three, two) != nine" << endl;
      success = 0;
    }

    // This fails inexplicably for complex numbers on gcc 4.2.1 on Mac.
    if (! AT::is_complex) {
      result = AT::pow (three, three);
      if (result != twentySeven) {
        out << "AT::pow (three, three) = " << result
            << " != twentySeven = " << twentySeven << endl;
        success = 0;
      }
    }

    // These fail inexplicably for complex numbers on gcc 4.2.1 on Mac.
    if (AT::is_signed && ! AT::is_complex) {
      result = AT::pow (-three, one);
      if (result != -three) {
        out << "AT::pow (-three, one) = " << result
            << " != -three = " << -three << endl;
        success = 0;
      }
      result = AT::pow (-three, two);
      if (result != nine) {
        out << "AT::pow (-three, two) = " << result
            << " != nine = " << nine << endl;
        success = 0;
      }
      result = AT::pow (-three, three);
      if (result != -twentySeven) {
        out << "AT::pow (-three, three) = " << result
            << " != -twentySeven = " << twentySeven << endl;
        success = 0;
      }
    }

    if (AT::sqrt (zero) != zero) {
      out << "AT::sqrt (zero) != zero" << endl;
      success = 0;
    }
    if (AT::sqrt (one) != one) {
      out << "AT::sqrt (one) != one" << endl;
      success = 0;
    }
    if (AT::sqrt (thirtySix) != six) {
      out << "AT::sqrt (thirtySix) != six" << endl;
      success = 0;
    }
    if (AT::sqrt (sixtyFour) != eight) {
      out << "AT::sqrt (sixtyFour) != eight" << endl;
      success = 0;
    }
    if (AT::is_integer) {
      if (AT::sqrt (fortyTwo) != six) {
        out << "AT::sqrt (fortyTwo) != six" << endl;
        success = 0;
      }
      if (AT::sqrt (oneTwentySeven) != eleven) {
        out << "AT::sqrt (oneTwentySeven) != eleven" << endl;
        success = 0;
      }
    }

    if (!equal(AT::cbrt (zero) , zero)) {
      printf ("AT::cbrt(0) != 0\n");
      success = 0;
    }
    if (!equal(AT::cbrt (one) , one)) {
      printf ("AT::cbrt(1) != 1\n");
      success = 0;
    }
    if (!equal(AT::cbrt (twentySeven) , three)) {
      printf ("AT::cbrt(27) != 3\n");
      success = 0;
    }
    if (!equal(AT::cbrt (sixtyFour) , four)) {
      printf ("AT::cbrt(64) != 4\n");
      success = 0;
    }
    if (AT::is_integer) {
      if (!equal(AT::cbrt (fortyTwo) , three)) {
        printf ("AT:cbrt(42) != 3\n");
        success = 0;
      }
      if (!equal(AT::cbrt (oneTwentySeven) , five)) {
        printf ("AT::cbrt(127) != 5\n");
        success = 0;
      }
    }

    if (!equal(AT::exp (zero) , one)) {
      printf ("AT::cbrt(0) != 1\n");
      success = 0;
    }
    if (AT::is_complex) {
      const ScalarType val = two; //(two.real(), two.real());
      if (!equal(AT::conj (AT::exp  (val)) , 
                 AT::exp  (AT::conj (val)))) {
        printf ("AT::conj(exp(complex(2,0))) != AT::exp(conj(complex(2,0)))\n");
        success = 0;
      }
    }
    if (AT::log (one) != zero) {
      out << "AT::log (one) != zero" << endl;
      success = 0;
    }
    if (AT::log10 (one) != zero) {
      out << "AT::log10 (one) != zero" << endl;
      success = 0;
    }

    if (AT::is_complex) {
      const ScalarType val = two; // (two.real(), two.real());
      const auto val_sin = AT::sin (val);
      const auto val_cos = AT::cos (val);
      if (!equal(val_sin*val_sin + val_cos*val_cos , one)) {
        printf ("AT(complex):: sin(val)*sin(val) + cos(val)*cos(val) != 1\n");
        success = 0;
      } 
      if (!equal(val_sin/val_cos , AT::tan(val))) {
        printf ("AT(complex):: sin(val)/cos(val) != AT(real)::tan(val)\n");
        success = 0;
      } 
    } else {
      const ScalarType val = three; 
      const auto val_sin = AT::sin (val);
      const auto val_cos = AT::cos (val);
      if (!equal(val_sin*val_sin + val_cos*val_cos , one)) {
        printf ("AT(real):: sin(val)*sin(val) + cos(a)*cos(a) != 1\n");
        success = 0;
      } 
      if (!equal(val_sin/val_cos , AT::tan(val))) {
        printf ("AT(real):: sin(val)/cos(val) != AT(real)::tan(val)\n");
        success = 0;
      } 
    }

    if (!equal(AT::asin (AT::sin (three)), three)) {
      printf ("AT::asin(sin(3)) != 3\n");
      success = 0;
    } 
    if (!equal(AT::acos (AT::cos (three)), three)) {
      printf ("AT::acos(cos(3)) != 3\n");
      success = 0;
    } 
    if (!equal(AT::atan (AT::tan (three)), three)) {
      printf ("AT::atan(tan(3)) != 3\n");
      success = 0;
    } 

    // Call the base class' implementation.  Every subclass'
    // implementation of testHostImpl() should (must) do this, in
    // order to include the parent class' tests.  In the case of this
    // particular class, the base class' implementation doesn't do
    // anything, but that's OK.
    const int parentSuccess = base_type::testHostImpl (out);
    success = success && parentSuccess;

    return success;
  }
};


/// \class ArithTraitsTesterComplexBase
/// \brief Execute Kokkos::Details::ArithTraits tests relevant to
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
template<class ScalarType,
         class DeviceType,
         const int is_complex = Kokkos::Details::ArithTraits<ScalarType>::is_complex>
class ArithTraitsTesterComplexBase :
  public ArithTraitsTesterTranscendentalBase<ScalarType, DeviceType> {
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
  KOKKOS_INLINE_FUNCTION void
  operator () (size_type iwork, value_type& dst) const;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterComplexBase ();

protected:
  // The host hook gets implemented in the complex-arithmetic
  // specialization of this class.
  virtual int testHostImpl (std::ostream& out) const;
};

//
// Specialization of ArithTraitsTesterComplexBase for real T.
//
template<class ScalarType,
         class DeviceType>
class ArithTraitsTesterComplexBase<ScalarType, DeviceType, 0> :
  public ArithTraitsTesterTranscendentalBase<ScalarType, DeviceType> {
private:
  //! The base class of this class.
  typedef ArithTraitsTesterTranscendentalBase<ScalarType, DeviceType> base_type;

public:
  typedef DeviceType execution_space;
  typedef typename execution_space::size_type size_type;
  //! Type of the result of the reduction.
  typedef int value_type;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterComplexBase () {}

  KOKKOS_INLINE_FUNCTION void
  operator () (size_type iwork, value_type& dst) const {
    typedef Kokkos::Details::ArithTraits<ScalarType> AT;
    (void) iwork; // forestall compiler warning for unused variable
    int success = 1;

    // Apparently, std::numeric_limits<ScalarType>::is_signed is 1
    // only for real numbers.
    if (AT::is_signed != std::numeric_limits<ScalarType>::is_signed) {
      success = 0;
    }
    if (AT::is_complex) {
      success = 0;
    }

    // Call the base class' implementation.  Every subclass'
    // implementation of operator() must do this, in order to include
    // the parent class' tests.
    int baseResult = 1;
    base_type::operator () (iwork, baseResult);
    success = success && baseResult;

    dst = dst && success;
  }

protected:
  virtual int testHostImpl (std::ostream& out) const {
    using std::endl;
    typedef Kokkos::Details::ArithTraits<ScalarType> AT;

    int success = 1;
    // Apparently, std::numeric_limits<ScalarType>::is_signed is 1 only for real numbers.
    if (AT::is_signed != std::numeric_limits<ScalarType>::is_signed) {
      out << "ArithTraits<T>::is_signed != std::numeric_limits<ScalarType>::is_signed" << endl;
      success = 0;
    }
    if (AT::is_complex) {
      out << "ArithTraits<T>::is_complex is wrong" << endl;
      success = 0;
    }
    // Call the base class' implementation.  Every subclass'
    // implementation of testHostImpl() should (must) do this, in
    // order to include the parent class' tests.  In the case of this
    // particular class, the base class' implementation doesn't do
    // anything, but that's OK.
    const int parentSuccess = base_type::testHostImpl (out);
    success = success && parentSuccess;

    return success;
  }
};

// Specialization for complex T.
template<class ScalarType,
         class DeviceType>
class ArithTraitsTesterComplexBase<ScalarType, DeviceType, 1> :
  public ArithTraitsTesterTranscendentalBase<ScalarType, DeviceType> {
private:
  //! The base class of this class.
  typedef ArithTraitsTesterTranscendentalBase<ScalarType, DeviceType> base_type;

public:
  typedef DeviceType execution_space;
  typedef typename execution_space::size_type size_type;
  //! Type of the result of the reduction.
  typedef int value_type;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterComplexBase () {}

  KOKKOS_INLINE_FUNCTION void
  operator () (size_type iwork, value_type& dst) const {
    typedef Kokkos::Details::ArithTraits<ScalarType> AT;
    (void) iwork; // forestall compiler warning for unused variable
    int success = 1;

    if (! AT::is_complex) {
      success = 0;
    }
    typedef typename AT::mag_type mag_type;
    const mag_type one = Kokkos::Details::ArithTraits<mag_type>::one ();

    // This presumes that ScalarType, being a complex number, has a
    // constructor which takes two mag_type arguments.
    const ScalarType oneMinusOne (one, -one);
    const ScalarType onePlusOne (one, one);

    // Test conjugation.
    if (AT::conj (oneMinusOne) != onePlusOne ||
        AT::conj (onePlusOne) != oneMinusOne) {
      success = 0;
    }

    // Call the base class' implementation.  Every subclass'
    // implementation of operator() must do this, in order to include
    // the parent class' tests.
    int baseResult = 1;
    base_type::operator () (iwork, baseResult);
    success = success && baseResult;

    dst = dst && success;
  }

protected:
  virtual int testHostImpl (std::ostream& out) const {
    using std::endl;
    typedef Kokkos::Details::ArithTraits<ScalarType> AT;
    int success = 1;

    if (! AT::is_complex) {
      out << "ArithTraits<T>::is_complex is wrong" << endl;
      success = 0;
    }
    typedef typename AT::mag_type mag_type;
    const mag_type one = Kokkos::Details::ArithTraits<mag_type>::one ();

    // This presumes that ScalarType, being a complex number, has a
    // constructor which takes two mag_type arguments.
    const ScalarType oneMinusOne (one, -one);
    const ScalarType onePlusOne (one, one);

    // Test conjugation.
    if (AT::conj (oneMinusOne) != onePlusOne) {
      out << "AT::conj ((1, -1)) != (1, 1)" << endl;
      success = 0;
    }
    if (AT::conj (onePlusOne) != oneMinusOne) {
      out << "AT::conj ((1, 1)) != (1, -1)" << endl;
      success = 0;
    }
    // Call the base class' implementation.  Every subclass'
    // implementation of testHostImpl() should (must) do this, in
    // order to include the parent class' tests.  In the case of this
    // particular class, the base class' implementation doesn't do
    // anything, but that's OK.
    const int parentSuccess = base_type::testHostImpl (out);
    success = success && parentSuccess;

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
/// Kokkos::Details::ArithTraits relevant to floating-point types.
///
/// The tests include those suitable for execution on the parallel
/// device (operator()) and those suitable for execution on the host
/// (testHost()).  The device-based test is a reduction over redundant
/// executions of the test.  All redundant executions must return
/// '1' (passed).
template<class ScalarType,
         class DeviceType,
         const int is_exact = Kokkos::Details::ArithTraits<ScalarType>::is_exact>
class ArithTraitsTesterFloatingPointBase :
  public ArithTraitsTesterComplexBase<ScalarType,
                                      DeviceType,
                                      Kokkos::Details::ArithTraits<ScalarType>::is_complex>
{
private:
  //! The base class of this class.
  typedef ArithTraitsTesterComplexBase<ScalarType,
                                       DeviceType,
                                       Kokkos::Details::ArithTraits<ScalarType>::is_complex> base_type;
public:
  typedef DeviceType execution_space;
  typedef typename execution_space::size_type size_type;
  //! Type of the result of the reduction.
  typedef int value_type;

  /// \brief The "parallel for" part of the reduction.
  ///
  /// See comments of ArithTraitsTesterBase's operator().
  KOKKOS_INLINE_FUNCTION void
  operator () (size_type iwork, value_type& dst) const;

protected:
  virtual int testHostImpl (std::ostream& out) const;
};


//
// Specialization for is_exact = 0 (i.e., ScalarType is a
// floating-point type).
//
template<class ScalarType, class DeviceType>
class ArithTraitsTesterFloatingPointBase<ScalarType, DeviceType, 0> :
  public ArithTraitsTesterComplexBase<ScalarType,
                                      DeviceType,
                                      Kokkos::Details::ArithTraits<ScalarType>::is_complex>
{
private:
  //! The base class of this class.
  typedef ArithTraitsTesterComplexBase<ScalarType,
                                       DeviceType,
                                       Kokkos::Details::ArithTraits<ScalarType>::is_complex> base_type;
public:
  typedef DeviceType execution_space;
  typedef typename execution_space::size_type size_type;
  //! Type of the result of the reduction.
  typedef int value_type;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterFloatingPointBase () {}

  KOKKOS_INLINE_FUNCTION void
  operator () (size_type iwork, value_type& dst) const {
    typedef Kokkos::Details::ArithTraits<ScalarType> AT;
    (void) iwork; // forestall compiler warning for unused variable
    int success = 1;

    if (AT::is_exact) {
      printf ("AT::is_exact is 1\n");
      success = 0;
    }
    if (! AT::isNan (AT::nan ())) {
      printf ("NaN is not NaN\n");
      success = 0;
    }

    const ScalarType zero = AT::zero ();
    const ScalarType one = AT::one ();

    if (AT::isInf (zero)) {
      printf ("0 is Inf\n");
      success = 0;
    }
    if (AT::isInf (one)) {
      printf ("1 is Inf\n");
      success = 0;
    }
    if (AT::isNan (zero)) {
      printf ("0 is NaN\n");
      success = 0;
    }
    if (AT::isNan (one)) {
      printf ("1 is NaN\n");
      success = 0;
    }

    // Call the base class' implementation.  Every subclass'
    // implementation of operator() must do this, in order to include
    // the parent class' tests.
    int baseResult = 1;
    base_type::operator () (iwork, baseResult);
    success = success && baseResult;

    dst = dst && success;
  }

protected:
  virtual int testHostImpl (std::ostream& out) const {
    typedef Kokkos::Details::ArithTraits<ScalarType> AT;
    using std::endl;
    int success = 1;

    if (AT::is_exact) {
      out << "AT::is_exact is wrong" << endl;
      success = 0;
    }

    //if (std::numeric_limits<ScalarType>::is_iec559) {
    //success = success && AT::isInf (AT::inf ());
    if (! AT::isNan (AT::nan ())) {
      out << "isNan or nan failed" << endl;
      success = 0;
    }
    //}

    const ScalarType zero = AT::zero ();
    const ScalarType one = AT::one ();

    if (AT::isInf (zero)) {
      out << "isInf(zero) is 1" << endl;
      success = 0;
    }
    if (AT::isInf (one)) {
      out << "isInf(one) is 1" << endl;
      success = 0;
    }
    if (AT::isNan (zero)) {
      out << "isNan(zero) is 1" << endl;
      success = 0;
    }
    if (AT::isNan (one)) {
      out << "isNan(one) is 1" << endl;
      success = 0;
    }

    // Call the base class' implementation.  Every subclass'
    // implementation of testHostImpl() should (must) do this, in
    // order to include the parent class' tests.
    const int parentSuccess = base_type::testHostImpl (out);
    success = success && parentSuccess;

    return success;
  }
};


//
// Specialization for is_exact = 1 (i.e., ScalarType is <i>not</i>
// a floating-point type).
//
template<class ScalarType, class DeviceType>
class ArithTraitsTesterFloatingPointBase<ScalarType, DeviceType, 1> :
  public ArithTraitsTesterComplexBase<ScalarType,
                                      DeviceType,
                                      Kokkos::Details::ArithTraits<ScalarType>::is_complex>
{
private:
  //! The base class of this class.
  typedef ArithTraitsTesterComplexBase<ScalarType,
                                       DeviceType,
                                       Kokkos::Details::ArithTraits<ScalarType>::is_complex> base_type;
public:
  typedef DeviceType execution_space;
  typedef typename execution_space::size_type size_type;
  //! Type of the result of the reduction.
  typedef int value_type;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterFloatingPointBase () {}

  KOKKOS_INLINE_FUNCTION void
  operator () (size_type iwork, value_type& dst) const {
    typedef Kokkos::Details::ArithTraits<ScalarType> AT;
    (void) iwork; // forestall compiler warning for unused variable
    int success = 1;

    if (! AT::is_exact) {
      printf ("! AT:is_exact\n");
      success = 0;
    }

    // Call the base class' implementation.  Every subclass'
    // implementation of operator() must do this, in order to include
    // the parent class' tests.
    int baseResult = 1;
    base_type::operator () (iwork, baseResult);
    success = success && baseResult;

    dst = dst && success;
  }

protected:
  virtual int testHostImpl (std::ostream& out) const {
    typedef Kokkos::Details::ArithTraits<ScalarType> AT;
    using std::endl;
    int success = 1;

    if (! AT::is_exact) {
      out << "AT::is_exact is wrong" << endl;
      success = 0;
    }
    // Call the base class' implementation.  Every subclass'
    // implementation of testHostImpl() should (must) do this, in
    // order to include the parent class' tests.
    const int parentSuccess = base_type::testHostImpl (out);
    success = success && parentSuccess;

    return success;
  }
};


/// \class ArithTraitsTester
/// \brief Tests for Kokkos::Details::ArithTraits
/// \tparam ScalarType Any type for which Kokkos::Details::ArithTraits
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
/// operator for testing Kokkos::Details::ArithTraits.  This test
/// works for any type <tt>ScalarType</tt> for which
/// Kokkos::Details::ArithTraits has a specialization, and which can
/// be executed on the parallel device.
///
/// The tests include those suitable for execution on the parallel
/// device (operator()) and those suitable for execution on the host
/// (testHost()).  The device-based test is a reduction over redundant
/// executions of the test.  All redundant executions must return
/// '1' (passed).
template<class ScalarType, class DeviceType>
class ArithTraitsTester :
  public ArithTraitsTesterFloatingPointBase<ScalarType, DeviceType> {
public:
  typedef DeviceType execution_space;
  typedef typename execution_space::size_type size_type;
  //! Type of the result of the reduction.
  typedef int value_type;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTester () {}
};


/// \brief Run the Kokkos::Details::ArithTraits tests on the parallel device.
/// \tparam ScalarType Any type for which Kokkos::Details::ArithTraits
///   has a specialization, and which can be executed on the parallel
///   device.
/// \tparam DeviceType A Kokkos parallel device type.
///
/// This function must be called on the host, but it executes on the
/// device with (redundant) parallelism.
///
/// \return \c 1 if all redundant executions pass, else \c 0.
template<class ScalarType, class DeviceType>
int testArithTraitsOnDevice (std::ostream& out, const int verbose)
{
  using std::endl;
  typedef ArithTraitsTester<ScalarType, DeviceType> functor_type;
  int success = 1; // output argument of parallel_reduce
  Kokkos::parallel_reduce (1, functor_type (), success);
  if (success) {
    if (verbose)
      out << Kokkos::Details::ArithTraits<ScalarType>::name () << " passed" << endl;
  } else {
    out << Kokkos::Details::ArithTraits<ScalarType>::name () << " FAILED" << endl;
  }
  return success;
}

/// \brief Run the Kokkos::Details::ArithTraits tests on the host.
/// \tparam ScalarType Any type for which Kokkos::Details::ArithTraits
///   has a specialization.
/// \tparam DeviceType A Kokkos parallel device type.
///
/// This function must be called on the host, and executes on the host.
///
/// \return \c 1 if all tests pass, else \c 0.
template<class ScalarType, class DeviceType>
int testArithTraitsOnHost (std::ostream& out, const int verbose)
{
  using std::endl;
  ArithTraitsTester<ScalarType, DeviceType> f;
  const int localSuccess = f.testHost (out);

  if (localSuccess) {
    if (verbose)
      out << Kokkos::Details::ArithTraits<ScalarType>::name () << " passed" << endl;
  } else {
    out << Kokkos::Details::ArithTraits<ScalarType>::name () << " FAILED" << endl;
  }
  return localSuccess;
}


/// \brief Run the Kokkos::Details::ArithTraits tests for all (valid)
///   scalar types, on the given parallel device.
/// \tparam DeviceType A Kokkos parallel device type.
///
/// This is the "outward-facing" function meant to be called by
/// main().  This function must be called on the host, but it executes
/// on the device with (redundant) parallelism.
///
/// \return \c 1 if all tests pass, else \c 0.
template<class DeviceType>
int runAllArithTraitsDeviceTests (std::ostream& out, const int verbose)
{
  int success = 1;
  int curSuccess = 1;
  //
  // Built-in char(acter) types
  //

  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<char, DeviceType> (out, verbose);
  // Interestingly enough, char and int8_t are different types, but
  // signed char and int8_t are the same (on my system).
  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<signed char, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<unsigned char, DeviceType> (out, verbose);

  //
  // Built-in integer types
  //

  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<short, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<unsigned short, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<int8_t, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<uint8_t, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<int16_t, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<uint16_t, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<int32_t, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<uint32_t, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<int, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<unsigned int, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<int64_t, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<uint64_t, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<long, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<unsigned long, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<long long, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<unsigned long long, DeviceType> (out, verbose);

  //
  // Built-in real floating-point types
  //

  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<float, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<double, DeviceType> (out, verbose);

  //
  // Kokkos' complex floating-point types
  //

  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<Kokkos::complex<float>, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnDevice<Kokkos::complex<double>, DeviceType> (out, verbose);

  return success;
}


/// \brief Run the Kokkos::Details::ArithTraits tests for all scalar
///   types, on the host.
/// \tparam DeviceType A Kokkos parallel device type.
///
/// This is the "outward-facing" function meant to be called by
/// main().  This function must be called on the host, and executes on
/// the host.
///
/// \return \c 1 if all tests pass, else \c 0.
template<class DeviceType>
int runAllArithTraitsHostTests (std::ostream& out, const int verbose)
{
  int success = 1;
  int curSuccess = 1;

  //
  // Built-in char(acter) types
  //

  success = success && curSuccess; curSuccess = testArithTraitsOnHost<char, DeviceType> (out, verbose);
  // Interestingly enough, char and int8_t are different types, but
  // signed char and int8_t are the same (on my system).
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<signed char, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<unsigned char, DeviceType> (out, verbose);

  //
  // Built-in integer types
  //

  success = success && curSuccess; curSuccess = testArithTraitsOnHost<short, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<unsigned short, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<int8_t, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<uint8_t, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<int16_t, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<uint16_t, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<int32_t, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<uint32_t, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<int, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<unsigned int, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<int64_t, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<uint64_t, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<long, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<unsigned long, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<long long, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<unsigned long long, DeviceType> (out, verbose);

  //
  // Built-in real and complex floating-point types
  //

  success = success && curSuccess; curSuccess = testArithTraitsOnHost<float, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<double, DeviceType> (out, verbose);
#ifndef KOKKOS_ENABLE_CUDA
  // This would spill tons of warnings about host device stuff otherwise
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<long double, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<std::complex<float>, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<std::complex<double>, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<std::complex<long double>, DeviceType> (out, verbose);
#endif
  //
  // Kokkos' complex floating-point types
  //

  success = success && curSuccess; curSuccess = testArithTraitsOnHost<Kokkos::complex<float>, DeviceType> (out, verbose);
  success = success && curSuccess; curSuccess = testArithTraitsOnHost<Kokkos::complex<double>, DeviceType> (out, verbose);
  //success = success && curSuccess; curSuccess = testArithTraitsOnHost<Kokkos::complex<long double>, DeviceType> (out, verbose);

  return success;
}

template <typename device>
void test_ArithTraits ()
{
  using std::endl;

  class NullBuffer : public std::streambuf
  {
  public:
    int overflow(int c) { return c; }
  };
  NullBuffer null_buffer;
  //std::ostream &out = std::cout;
  std::ostream out(&null_buffer);


  bool success = true;
  success = runAllArithTraitsDeviceTests <device>(out, 0);
  EXPECT_TRUE( success);
  success = runAllArithTraitsHostTests <device>(out, 0);
  EXPECT_TRUE( success);
}
TEST_F( TestCategory, common_ArithTraits) {
  test_ArithTraits<TestExecSpace>();
}

#endif // KOKKOS_ARITHTRAITSTEST_HPP
