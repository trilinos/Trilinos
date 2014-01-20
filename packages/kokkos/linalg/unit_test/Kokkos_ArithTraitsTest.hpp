/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
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

#include <Kokkos_ParallelReduce.hpp>
#include "Kokkos_ArithTraits.hpp"
#include <limits> // std::numeric_limits
#include <typeinfo> // typeid (T)

// If KOKKOS_INLINE_FUNCTION isn't already defined from Kokkos, define
// it here.  If compiling with CUDA, the macro includes both the
// __host__ and __device__ attributes.  Whether or not compiling with
// CUDA, the macro also includes the inline attribute.
#ifndef KOKKOS_INLINE_FUNCTION
#  ifdef __CUDA_ARCH__
#    define KOKKOS_INLINE_FUNCTION inline __host__ __device__
#  else
#    define KOKKOS_INLINE_FUNCTION inline
#  endif // __CUDA_ARCH__
#endif // KOKKOS_INLINE_FUNCTION

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
/// 'true' (passed).
template<class ScalarType, class DeviceType>
class ArithTraitsTesterBase {
public:
  typedef DeviceType device_type;
  typedef typename device_type::size_type size_type;
  //! Type of the result of the reduction.
  typedef bool value_type;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterBase () {}

  /// \brief Set the initial value (\c true) of the reduction.
  ///
  /// Subclasses need not and must not override this method.
  KOKKOS_INLINE_FUNCTION void init (value_type& dst) const {
    dst = true;
  }

  /// \brief Combine two intermediate reduction results into \c dst.
  ///
  /// Subclasses need not and must not override this method.
  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& dst,
        const volatile value_type& src) const
  {
    dst = dst && src;
  }

  /// \brief The "parallel for" part of the reduction.
  ///
  /// This is the method that actually runs the tests on the device.
  /// It runs through a sequence of tests, and produces a \c true
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
    bool success = true;

    // Make sure that the typedef exists.
    typedef typename AT::mag_type mag_type;

    // ArithTraits should not even compile if it's not specialized for
    // T, but we check for this bool constant for compatibility with
    // std::numeric_limits.
    if (! AT::is_specialized) {
      success = false;
    }

    // It's OK to refer to std::numeric_limits constants in a device
    // function, just not to its class methods (which are not marked
    // as device functions).
    if (AT::is_integer != std::numeric_limits<ScalarType>::is_integer) {
      success = false;
    }
    if (AT::is_exact != std::numeric_limits<ScalarType>::is_exact) {
      success = false;
    }

    const ScalarType zero = AT::zero ();
    const ScalarType one = AT::one ();

    // Test properties of the arithmetic and multiplicative identities.
    if (zero + zero != zero) {
      success = false;
    }
    if (zero + one != one) {
      success = false;
    }
    if (one - one != zero) {
      success = false;
    }
    // This is technically true even of Z_2, since in that field, one
    // is its own inverse (so -one == one).
    if ((one + one) - one != one) {
      success = false;
    }

    if (AT::abs (zero) != zero) {
      success = false;
    }
    if (AT::abs (one) != one) {
      success = false;
    }
    if (AT::is_signed && AT::abs (-one) != one) {
      success = false;
    }
    // Need enable_if to test whether T can be compared using <=.
    // However, mag_type should always be comparable using <=.
    //
    // These are very mild ordering properties.
    // They should work even for a set only containing zero.
    if (AT::abs (zero) > AT::abs (AT::max ())) {
      success = false;
    }

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
        success = false;
      }
    }
    if (AT::pow (three, zero) != one) {
      success = false;
    }
    if (AT::pow (three, one) != three) {
      success = false;
    }
    if (AT::pow (three, two) != nine) {
      success = false;
    }

    // This fails inexplicably for complex numbers on gcc 4.2.1 on Mac.
    if (! AT::is_complex) {
      result = AT::pow (three, three);
      if (result != twentySeven) {
        success = false;
      }
    }

    // These fail inexplicably for complex numbers on gcc 4.2.1 on Mac.
    if (AT::is_signed && ! AT::is_complex) {
      result = AT::pow (-three, one);
      if (result != -three) {
        success = false;
      }
      result = AT::pow (-three, two);
      if (result != nine) {
        success = false;
      }
      result = AT::pow (-three, three);
      if (result != -twentySeven) {
        success = false;
      }
    }

    if (AT::sqrt (zero) != zero) {
      success = false;
    }
    if (AT::sqrt (one) != one) {
      success = false;
    }
    if (AT::sqrt (thirtySix) != six) {
      success = false;
    }
    if (AT::sqrt (sixtyFour) != eight) {
      success = false;
    }
    if (AT::is_integer) {
      success = success && (AT::sqrt (fortyTwo) == six);
      success = success && (AT::sqrt (oneTwentySeven) == eleven);
    }

    if (AT::log (one) != zero) {
      success = false;
    }
    if (AT::log10 (one) != zero) {
      success = false;
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
  /// \return \c true if all tests succeeded, else \c false.
  bool testHostImpl (std::ostream& out) const {
    return true; // there are no tests, so trivially, all the tests pass
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
  /// \return \c true if all the tests pass, else \c false.
  bool testHost (std::ostream& out) const {
    typedef Kokkos::Details::ArithTraits<ScalarType> AT;
    using std::endl;
    bool success = true;

    // Make sure that the typedef exists.
    typedef typename AT::mag_type mag_type;

    // ArithTraits should not even compile if it's not specialized for
    // T, but we check for this bool constant for compatibility with
    // std::numeric_limits.
    if (! AT::is_specialized) {
      out << "ArithTraits is not specialized for T" << endl;
      success = false;
    }

    if (AT::is_integer != std::numeric_limits<ScalarType>::is_integer) {
      out << "AT::is_integer != std::numeric_limits<ScalarType>::is_integer" << endl;
      success = false;
    }

    if (AT::is_exact != std::numeric_limits<ScalarType>::is_exact) {
      out << "AT::is_exact != std::numeric_limits<ScalarType>::is_exact" << endl;
      success = false;
    }

    const ScalarType zero = AT::zero ();
    const ScalarType one = AT::one ();
    // Test properties of the arithmetic and multiplicative identities.

    if (zero + zero != zero) {
      out << "zero + zero != zero" << endl;
      success = false;
    }
    if (zero + one != one) {
      out << "zero + one != one" << endl;
      success = false;
    }
    if (one - one != zero) {
      out << "one - one != zero" << endl;
      success = false;
    }
    // This is technically true even of Z_2, since in that field, one
    // is its own inverse (so -one == one).
    if ((one + one) - one != one) {
      out << "(one + one) - one != one" << endl;
      success = false;
    }

    if (AT::abs (zero) != zero) {
      out << "AT::abs (zero) != zero" << endl;
      success = false;
    }
    if (AT::abs (one) != one) {
      out << "AT::abs (one) != one" << endl;
      success = false;
    }
    if (AT::is_signed) {
      if (AT::abs (-one) != one) {
        out << "AT::abs (-one) != one" << endl;
        success = false;
      }
    }
    // Need enable_if to test whether T can be compared using <=.
    // However, mag_type should always be comparable using <=.
    //
    // // These are very mild ordering properties.
    // // They should work even for a set only containing zero.
    if (AT::abs (zero) > AT::abs (AT::max ())) {
      out << "AT::abs (zero) > AT::abs (AT::max ())" << endl;
      success = false;
    }
    //success = success && (AT::abs (AT::min ()) <= AT::abs (AT::max ()));

    // Need enable_if to do a complex test.
    // if (AT::is_complex) {

    // }
    // else {

    // }

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
        success = false;
      }
    }
    if (AT::pow (three, zero) != one) {
      out << "AT::pow (three, zero) != one" << endl;
      success = false;
    }
    if (AT::pow (three, one) != three) {
      out << "AT::pow (three, one) != three" << endl;
      success = false;
    }
    if (AT::pow (three, two) != nine) {
      out << "AT::pow (three, two) != nine" << endl;
      success = false;
    }

    // This fails inexplicably for complex numbers on gcc 4.2.1 on Mac.
    if (! AT::is_complex) {
      result = AT::pow (three, three);
      if (result != twentySeven) {
        out << "AT::pow (three, three) = " << result
            << " != twentySeven = " << twentySeven << endl;
        success = false;
      }
    }

    // These fail inexplicably for complex numbers on gcc 4.2.1 on Mac.
    if (AT::is_signed && ! AT::is_complex) {
      result = AT::pow (-three, one);
      if (result != -three) {
        out << "AT::pow (-three, one) = " << result
            << " != -three = " << -three << endl;
        success = false;
      }
      result = AT::pow (-three, two);
      if (result != nine) {
        out << "AT::pow (-three, two) = " << result
            << " != nine = " << nine << endl;
        success = false;
      }
      result = AT::pow (-three, three);
      if (result != -twentySeven) {
        out << "AT::pow (-three, three) = " << result
            << " != -twentySeven = " << twentySeven << endl;
        success = false;
      }
    }

    if (AT::sqrt (zero) != zero) {
      out << "AT::sqrt (zero) != zero" << endl;
      success = false;
    }
    if (AT::sqrt (one) != one) {
      out << "AT::sqrt (one) != one" << endl;
      success = false;
    }
    if (AT::sqrt (thirtySix) != six) {
      out << "AT::sqrt (thirtySix) != six" << endl;
      success = false;
    }
    if (AT::sqrt (sixtyFour) != eight) {
      out << "AT::sqrt (sixtyFour) != eight" << endl;
      success = false;
    }
    if (AT::is_integer) {
      success = success && (AT::sqrt (fortyTwo) == six);
      success = success && (AT::sqrt (oneTwentySeven) == eleven);
    }

    if (AT::log (one) != zero) {
      out << "AT::log (one) != zero" << endl;
      success = false;
    }
    if (AT::log10 (one) != zero) {
      out << "AT::log10 (one) != zero" << endl;
      success = false;
    }

    // Run the subclass' remaining tests, if any.
    success = success && testHostImpl (out);

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
         const bool is_complex = Kokkos::Details::ArithTraits<ScalarType>::is_complex>
class ArithTraitsTesterComplexBase : public ArithTraitsTesterBase<ScalarType, DeviceType> {
private:
  //! The base class of this class.
  typedef ArithTraitsTesterBase<ScalarType, DeviceType> base_type;

public:
  typedef DeviceType device_type;
  typedef typename device_type::size_type size_type;
  //! Type of the result of the reduction.
  typedef bool value_type;

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
  virtual bool testHostImpl (std::ostream& out) const;
};

//
// Specialization of ArithTraitsTesterComplexBase for real T.
//
template<class ScalarType,
         class DeviceType>
class ArithTraitsTesterComplexBase<ScalarType, DeviceType, false> :
  public ArithTraitsTesterBase<ScalarType, DeviceType> {
private:
  //! The base class of this class.
  typedef ArithTraitsTesterBase<ScalarType, DeviceType> base_type;

public:
  typedef DeviceType device_type;
  typedef typename device_type::size_type size_type;
  //! Type of the result of the reduction.
  typedef bool value_type;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterComplexBase () {}

  KOKKOS_INLINE_FUNCTION void
  operator () (size_type iwork, value_type& dst) const {
    typedef Kokkos::Details::ArithTraits<ScalarType> AT;
    (void) iwork; // forestall compiler warning for unused variable
    bool success = true;

    // Apparently, std::numeric_limits<ScalarType>::is_signed is true
    // only for real numbers.
    if (AT::is_signed != std::numeric_limits<ScalarType>::is_signed) {
      success = false;
    }
    if (AT::is_complex) {
      success = false;
    }

    // Call the base class' implementation.  Every subclass'
    // implementation of operator() must do this, in order to include
    // the parent class' tests.
    bool baseResult = true;
    base_type::operator () (iwork, baseResult);
    success = success && baseResult;

    dst = dst && success;
  }

protected:
  virtual bool testHostImpl (std::ostream& out) const {
    using std::endl;
    typedef Kokkos::Details::ArithTraits<ScalarType> AT;

    bool success = true;
    // Apparently, std::numeric_limits<ScalarType>::is_signed is true only for real numbers.
    if (AT::is_signed != std::numeric_limits<ScalarType>::is_signed) {
      out << "ArithTraits<T>::is_signed != std::numeric_limits<ScalarType>::is_signed" << endl;
      success = false;
    }
    if (AT::is_complex) {
      out << "ArithTraits<T>::is_complex is wrong" << endl;
      success = false;
    }
    // Call the base class' implementation.  Every subclass'
    // implementation of testHostImpl() should (must) do this, in
    // order to include the parent class' tests.  In the case of this
    // particular class, the base class' implementation doesn't do
    // anything, but that's OK.
    success = success && base_type::testHostImpl (out);

    return success;
  }
};

// Specialization for complex T.
template<class ScalarType,
         class DeviceType>
class ArithTraitsTesterComplexBase<ScalarType, DeviceType, true> :
  public ArithTraitsTesterBase<ScalarType, DeviceType> {
private:
  //! The base class of this class.
  typedef ArithTraitsTesterBase<ScalarType, DeviceType> base_type;

public:
  typedef DeviceType device_type;
  typedef typename device_type::size_type size_type;
  //! Type of the result of the reduction.
  typedef bool value_type;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterComplexBase () {}

  KOKKOS_INLINE_FUNCTION void
  operator () (size_type iwork, value_type& dst) const {
    typedef Kokkos::Details::ArithTraits<ScalarType> AT;
    (void) iwork; // forestall compiler warning for unused variable
    bool success = true;

    if (! AT::is_complex) {
      success = false;
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
      success = false;
    }

    // Call the base class' implementation.  Every subclass'
    // implementation of operator() must do this, in order to include
    // the parent class' tests.
    bool baseResult = true;
    base_type::operator () (iwork, baseResult);
    success = success && baseResult;

    dst = dst && success;
  }

protected:
  virtual bool testHostImpl (std::ostream& out) const {
    using std::endl;
    typedef Kokkos::Details::ArithTraits<ScalarType> AT;
    bool success = true;

    if (! AT::is_complex) {
      out << "ArithTraits<T>::is_complex is wrong" << endl;
      success = false;
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
      success = false;
    }
    if (AT::conj (onePlusOne) != oneMinusOne) {
      out << "AT::conj ((1, 1)) != (1, -1)" << endl;
      success = false;
    }
    // Call the base class' implementation.  Every subclass'
    // implementation of testHostImpl() should (must) do this, in
    // order to include the parent class' tests.  In the case of this
    // particular class, the base class' implementation doesn't do
    // anything, but that's OK.
    success = success && base_type::testHostImpl (out);

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
/// 'true' (passed).
template<class ScalarType,
         class DeviceType,
         const bool is_exact = Kokkos::Details::ArithTraits<ScalarType>::is_exact>
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
  typedef DeviceType device_type;
  typedef typename device_type::size_type size_type;
  //! Type of the result of the reduction.
  typedef bool value_type;

  /// \brief The "parallel for" part of the reduction.
  ///
  /// See comments of ArithTraitsTesterBase's operator().
  KOKKOS_INLINE_FUNCTION void
  operator () (size_type iwork, value_type& dst) const;

protected:
  virtual bool testHostImpl (std::ostream& out) const;
};


//
// Specialization for is_exact = false (i.e., ScalarType is a
// floating-point type).
//
template<class ScalarType, class DeviceType>
class ArithTraitsTesterFloatingPointBase<ScalarType, DeviceType, false> :
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
  typedef DeviceType device_type;
  typedef typename device_type::size_type size_type;
  //! Type of the result of the reduction.
  typedef bool value_type;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterFloatingPointBase () {}

  KOKKOS_INLINE_FUNCTION void
  operator () (size_type iwork, value_type& dst) const {
    typedef Kokkos::Details::ArithTraits<ScalarType> AT;
    (void) iwork; // forestall compiler warning for unused variable
    bool success = true;

    if (AT::is_exact) {
      success = false;
    }
    if (! AT::isNan (AT::nan ())) {
      success = false;
    }

    const ScalarType zero = AT::zero ();
    const ScalarType one = AT::one ();

    if (AT::isInf (zero)) {
      success = false;
    }
    if (AT::isInf (one)) {
      success = false;
    }
    if (AT::isNan (zero)) {
      success = false;
    }
    if (AT::isNan (one)) {
      success = false;
    }

    // Call the base class' implementation.  Every subclass'
    // implementation of operator() must do this, in order to include
    // the parent class' tests.
    bool baseResult = true;
    base_type::operator () (iwork, baseResult);
    success = success && baseResult;

    dst = dst && success;
  }

protected:
  virtual bool testHostImpl (std::ostream& out) const {
    typedef Kokkos::Details::ArithTraits<ScalarType> AT;
    using std::endl;
    bool success = true;

    if (AT::is_exact) {
      out << "AT::is_exact is wrong" << endl;
      success = false;
    }

    //if (std::numeric_limits<ScalarType>::is_iec559) {
    //success = success && AT::isInf (AT::inf ());
    success = success && AT::isNan (AT::nan ());
    if (! success) {
      out << "isNaN or nan failed" << endl;
    }
    //}

    const ScalarType zero = AT::zero ();
    const ScalarType one = AT::one ();

    if (AT::isInf (zero)) {
      out << "isInf(zero) is true" << endl;
      success = false;
    }
    if (AT::isInf (one)) {
      out << "isInf(one) is true" << endl;
      success = false;
    }
    if (AT::isNan (zero)) {
      out << "isNan(zero) is true" << endl;
      success = false;
    }
    if (AT::isNan (one)) {
      out << "isNan(one) is true" << endl;
      success = false;
    }

    // Call the base class' implementation.  Every subclass'
    // implementation of testHostImpl() should (must) do this, in
    // order to include the parent class' tests.
    success = success && base_type::testHostImpl (out);

    return success;
  }
};


//
// Specialization for is_exact = true (i.e., ScalarType is <i>not</i>
// a floating-point type).
//
template<class ScalarType, class DeviceType>
class ArithTraitsTesterFloatingPointBase<ScalarType, DeviceType, true> :
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
  typedef DeviceType device_type;
  typedef typename device_type::size_type size_type;
  //! Type of the result of the reduction.
  typedef bool value_type;

  //! Constructor (does nothing, but marked as device function).
  KOKKOS_INLINE_FUNCTION ArithTraitsTesterFloatingPointBase () {}

  KOKKOS_INLINE_FUNCTION void
  operator () (size_type iwork, value_type& dst) const {
    typedef Kokkos::Details::ArithTraits<ScalarType> AT;
    (void) iwork; // forestall compiler warning for unused variable
    bool success = true;

    if (! AT::is_exact) {
      success = false;
    }

    // Call the base class' implementation.  Every subclass'
    // implementation of operator() must do this, in order to include
    // the parent class' tests.
    bool baseResult = true;
    base_type::operator () (iwork, baseResult);
    success = success && baseResult;

    dst = dst && success;
  }

protected:
  virtual bool testHostImpl (std::ostream& out) const {
    typedef Kokkos::Details::ArithTraits<ScalarType> AT;
    using std::endl;
    bool success = true;

    if (! AT::is_exact) {
      out << "AT::is_exact is wrong" << endl;
      success = false;
    }
    // Call the base class' implementation.  Every subclass'
    // implementation of testHostImpl() should (must) do this, in
    // order to include the parent class' tests.
    success = success && base_type::testHostImpl (out);

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
/// 'true' (passed).
template<class ScalarType, class DeviceType>
class ArithTraitsTester :
  public ArithTraitsTesterFloatingPointBase<ScalarType, DeviceType> {
public:
  typedef DeviceType device_type;
  typedef typename device_type::size_type size_type;
  //! Type of the result of the reduction.
  typedef bool value_type;

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
/// \return \c true if all redundant executions pass, else \c false.
template<class ScalarType, class DeviceType>
bool testArithTraitsOnDevice (std::ostream& out, const bool verbose)
{
  using std::endl;
  typedef ArithTraitsTester<ScalarType, DeviceType> functor_type;
  bool success = true; // output argument of parallel_reduce
  Kokkos::parallel_reduce (1, functor_type (), success);
  if (success) {
    if (verbose) {
      out << typeid (ScalarType).name () << " passed" << endl;
    }
  } else {
    out << typeid (ScalarType).name () << " FAILED" << endl;
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
/// \return \c true if all tests pass, else \c false.
template<class ScalarType, class DeviceType>
bool testArithTraitsOnHost (std::ostream& out, const bool verbose)
{
  using std::endl;
  ArithTraitsTester<ScalarType, DeviceType> f;
  const bool localSuccess = f.testHost (out);
  if (localSuccess) {
    if (verbose) {
      out << typeid (ScalarType).name () << " passed" << endl;
    }
  } else {
    out << typeid (ScalarType).name () << " FAILED" << endl;
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
/// \return \c true if all tests pass, else \c false.
template<class DeviceType>
bool runAllArithTraitsDeviceTests (std::ostream& out, const bool verbose)
{
  bool success = true;

  //
  // Built-in char(acter) types
  //

  success = success && testArithTraitsOnDevice<char, DeviceType> (out, verbose);
  // Interestingly enough, char and int8_t are different types, but
  // signed char and int8_t are the same (on my system).
  success = success && testArithTraitsOnDevice<signed char, DeviceType> (out, verbose);
  success = success && testArithTraitsOnDevice<unsigned char, DeviceType> (out, verbose);

  //
  // Built-in integer types
  //

  success = success && testArithTraitsOnDevice<short, DeviceType> (out, verbose);
  success = success && testArithTraitsOnDevice<unsigned short, DeviceType> (out, verbose);
  success = success && testArithTraitsOnDevice<int8_t, DeviceType> (out, verbose);
  success = success && testArithTraitsOnDevice<uint8_t, DeviceType> (out, verbose);
  success = success && testArithTraitsOnDevice<int16_t, DeviceType> (out, verbose);
  success = success && testArithTraitsOnDevice<uint16_t, DeviceType> (out, verbose);
  success = success && testArithTraitsOnDevice<int32_t, DeviceType> (out, verbose);
  success = success && testArithTraitsOnDevice<uint32_t, DeviceType> (out, verbose);
  success = success && testArithTraitsOnDevice<int, DeviceType> (out, verbose);
  success = success && testArithTraitsOnDevice<unsigned int, DeviceType> (out, verbose);
  success = success && testArithTraitsOnDevice<int64_t, DeviceType> (out, verbose);
  success = success && testArithTraitsOnDevice<uint64_t, DeviceType> (out, verbose);
  success = success && testArithTraitsOnDevice<long, DeviceType> (out, verbose);
  success = success && testArithTraitsOnDevice<unsigned long, DeviceType> (out, verbose);
  success = success && testArithTraitsOnDevice<long long, DeviceType> (out, verbose);
  success = success && testArithTraitsOnDevice<unsigned long long, DeviceType> (out, verbose);

  //
  // Built-in real floating-point types
  //

  success = success && testArithTraitsOnDevice<float, DeviceType> (out, verbose);
  success = success && testArithTraitsOnDevice<double, DeviceType> (out, verbose);

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
/// \return \c true if all tests pass, else \c false.
template<class DeviceType>
bool runAllArithTraitsHostTests (std::ostream& out, const bool verbose)
{
  bool success = true;

  //
  // Built-in char(acter) types
  //

  success = success && testArithTraitsOnHost<char, DeviceType> (out, verbose);
  // Interestingly enough, char and int8_t are different types, but
  // signed char and int8_t are the same (on my system).
  success = success && testArithTraitsOnHost<signed char, DeviceType> (out, verbose);
  success = success && testArithTraitsOnHost<unsigned char, DeviceType> (out, verbose);

  //
  // Built-in integer types
  //

  success = success && testArithTraitsOnHost<short, DeviceType> (out, verbose);
  success = success && testArithTraitsOnHost<unsigned short, DeviceType> (out, verbose);
  success = success && testArithTraitsOnHost<int8_t, DeviceType> (out, verbose);
  success = success && testArithTraitsOnHost<uint8_t, DeviceType> (out, verbose);
  success = success && testArithTraitsOnHost<int16_t, DeviceType> (out, verbose);
  success = success && testArithTraitsOnHost<uint16_t, DeviceType> (out, verbose);
  success = success && testArithTraitsOnHost<int32_t, DeviceType> (out, verbose);
  success = success && testArithTraitsOnHost<uint32_t, DeviceType> (out, verbose);
  success = success && testArithTraitsOnHost<int, DeviceType> (out, verbose);
  success = success && testArithTraitsOnHost<unsigned int, DeviceType> (out, verbose);
  success = success && testArithTraitsOnHost<int64_t, DeviceType> (out, verbose);
  success = success && testArithTraitsOnHost<uint64_t, DeviceType> (out, verbose);
  success = success && testArithTraitsOnHost<long, DeviceType> (out, verbose);
  success = success && testArithTraitsOnHost<unsigned long, DeviceType> (out, verbose);
  success = success && testArithTraitsOnHost<long long, DeviceType> (out, verbose);
  success = success && testArithTraitsOnHost<unsigned long long, DeviceType> (out, verbose);

  //
  // Built-in real and complex floating-point types
  //

  success = success && testArithTraitsOnHost<float, DeviceType> (out, verbose);
  success = success && testArithTraitsOnHost<double, DeviceType> (out, verbose);
  success = success && testArithTraitsOnHost<long double, DeviceType> (out, verbose);
  success = success && testArithTraitsOnHost<std::complex<float>, DeviceType> (out, verbose);
  success = success && testArithTraitsOnHost<std::complex<double>, DeviceType> (out, verbose);
  success = success && testArithTraitsOnHost<std::complex<long double>, DeviceType> (out, verbose);

  return success;
}

#endif // KOKKOS_ARITHTRAITSTEST_HPP
