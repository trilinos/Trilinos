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

#ifndef KOKKOS_ARITHTRAITSTEST_HPP
#define KOKKOS_ARITHTRAITSTEST_HPP

#include "Kokkos_ArithTraits.hpp"
#include <limits> // std::numeric_limits

#ifndef KOKKOS_DEVICE_FUNCTION
#  ifdef __CUDA_ARCH__
#    define KOKKOS_DEVICE_FUNCTION inline __host__ __device__
#  else
#    define KOKKOS_DEVICE_FUNCTION
#  endif // __CUDA_ARCH__
#endif // KOKKOS_DEVICE_FUNCTION

template<class T>
class ArithTraitsTesterBase {
protected:
  static bool testFloatingPoint (std::ostream& out) {
    using Kokkos::Details::ArithTraits;
    using std::endl;
    bool success = true;

    //if (std::numeric_limits<T>::is_iec559) {
    //success = success && ArithTraits<T>::isInf (ArithTraits<T>::inf ());
    success = success && ArithTraits<T>::isNan (ArithTraits<T>::nan ());
    if (! success) {
      out << "isNaN or nan failed" << endl;
    }
    //}

    const T zero = ArithTraits<T>::zero ();
    const T one = ArithTraits<T>::one ();

    if (ArithTraits<T>::isInf (zero)) {
      out << "isInf(zero) is true" << endl;
      success = false;
    }
    if (ArithTraits<T>::isInf (one)) {
      out << "isInf(one) is true" << endl;
      success = false;
    }
    if (ArithTraits<T>::isNan (zero)) {
      out << "isNan(zero) is true" << endl;
      success = false;
    }
    if (ArithTraits<T>::isNan (one)) {
      out << "isNan(one) is true" << endl;
      success = false;
    }
    return success;
  }

  static bool test (std::ostream& out) {
    using Kokkos::Details::ArithTraits;
    using std::endl;
    bool success = true;

    // Make sure that the typedef exists.
    typedef typename ArithTraits<T>::mag_type mag_type;

    // ArithTraits should not even compile if it's not specialized for
    // T, but we check for this bool constant for compatibility with
    // std::numeric_limits.
    if (! ArithTraits<T>::is_specialized) {
      out << "ArithTraits is not specialized for T" << endl;
      success = false;
    }

    if (ArithTraits<T>::is_integer != std::numeric_limits<T>::is_integer) {
      out << "ArithTraits<T>::is_integer != std::numeric_limits<T>::is_integer" << endl;
      success = false;
    }

    if (ArithTraits<T>::is_exact != std::numeric_limits<T>::is_exact) {
      out << "ArithTraits<T>::is_exact != std::numeric_limits<T>::is_exact" << endl;
      success = false;
    }

    const T zero = ArithTraits<T>::zero ();
    const T one = ArithTraits<T>::one ();
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

    if (ArithTraits<T>::abs (zero) != zero) {
      out << "ArithTraits<T>::abs (zero) != zero" << endl;
      success = false;
    }
    if (ArithTraits<T>::abs (one) != one) {
      out << "ArithTraits<T>::abs (one) != one" << endl;
      success = false;
    }
    if (ArithTraits<T>::is_signed) {
      if (ArithTraits<T>::abs (-one) != one) {
        out << "ArithTraits<T>::abs (-one) != one" << endl;
        success = false;
      }
    }
    // Need enable_if to test whether T can be compared using <=.
    // However, mag_type should always be comparable using <=.
    //
    // // These are very mild ordering properties.
    // // They should work even for a set only containing zero.
    if (ArithTraits<T>::abs (zero) > ArithTraits<T>::abs (ArithTraits<T>::max ())) {
      out << "ArithTraits<T>::abs (zero) > ArithTraits<T>::abs (ArithTraits<T>::max ())" << endl;
      success = false;
    }
    //success = success && (ArithTraits<T>::abs (ArithTraits<T>::min ()) <= ArithTraits<T>::abs (ArithTraits<T>::max ()));

    // Need enable_if to do a complex test.
    // if (ArithTraits<T>::is_complex) {

    // }
    // else {

    // }

    const T two = one + one;
    const T three = one + one + one;
    const T four = two * two;
    const T five = four + one;
    const T six = three * two;
    const T seven = four + three;
    const T eight = four * two;
    const T nine = eight + one;
    const T eleven = five + six;
    const T twentySeven = nine * three;
    const T thirtySix = six * six;
    const T fortyTwo = six * seven;
    const T sixtyThree = eight * eight - one;
    const T sixtyFour = eight * eight;
    // max char value, for 8-bit char
    const T oneTwentySeven = sixtyFour + sixtyThree;

    T result;

    // This fails inexplicably for complex numbers on gcc 4.2.1 on Mac.
    if (! ArithTraits<T>::is_complex) {
      result = ArithTraits<T>::pow (two, three);
      if (result != eight) {
        out << "ArithTraits<T>::pow (two, three) != eight" << endl;
        success = false;
      }
    }
    if (ArithTraits<T>::pow (three, zero) != one) {
      out << "ArithTraits<T>::pow (three, zero) != one" << endl;
      success = false;
    }
    if (ArithTraits<T>::pow (three, one) != three) {
      out << "ArithTraits<T>::pow (three, one) != three" << endl;
      success = false;
    }
    if (ArithTraits<T>::pow (three, two) != nine) {
      out << "ArithTraits<T>::pow (three, two) != nine" << endl;
      success = false;
    }

    // This fails inexplicably for complex numbers on gcc 4.2.1 on Mac.
    if (! ArithTraits<T>::is_complex) {
      result = ArithTraits<T>::pow (three, three);
      if (result != twentySeven) {
        out << "ArithTraits<T>::pow (three, three) = " << result
            << " != twentySeven = " << twentySeven << endl;
        success = false;
      }
    }

    // These fail inexplicably for complex numbers on gcc 4.2.1 on Mac.
    if (ArithTraits<T>::is_signed && ! ArithTraits<T>::is_complex) {
      result = ArithTraits<T>::pow (-three, one);
      if (result != -three) {
        out << "ArithTraits<T>::pow (-three, one) = " << result
            << " != -three = " << -three << endl;
        success = false;
      }
      result = ArithTraits<T>::pow (-three, two);
      if (result != nine) {
        out << "ArithTraits<T>::pow (-three, two) = " << result
            << " != nine = " << nine << endl;
        success = false;
      }
      result = ArithTraits<T>::pow (-three, three);
      if (result != -twentySeven) {
        out << "ArithTraits<T>::pow (-three, three) = " << result
            << " != -twentySeven = " << twentySeven << endl;
        success = false;
      }
    }

    if (ArithTraits<T>::sqrt (zero) != zero) {
      out << "ArithTraits<T>::sqrt (zero) != zero" << endl;
      success = false;
    }
    if (ArithTraits<T>::sqrt (one) != one) {
      out << "ArithTraits<T>::sqrt (one) != one" << endl;
      success = false;
    }
    if (ArithTraits<T>::sqrt (thirtySix) != six) {
      out << "ArithTraits<T>::sqrt (thirtySix) != six" << endl;
      success = false;
    }
    if (ArithTraits<T>::sqrt (sixtyFour) != eight) {
      out << "ArithTraits<T>::sqrt (sixtyFour) != eight" << endl;
      success = false;
    }
    if (ArithTraits<T>::is_integer) {
      success = success && (ArithTraits<T>::sqrt (fortyTwo) == six);
      success = success && (ArithTraits<T>::sqrt (oneTwentySeven) == eleven);
    }

    if (ArithTraits<T>::log (one) != zero) {
      out << "ArithTraits<T>::log (one) != zero" << endl;
      success = false;
    }
    if (ArithTraits<T>::log10 (one) != zero) {
      out << "ArithTraits<T>::log10 (one) != zero" << endl;
      success = false;
    }

    return success;
  }
};


template<class T, const bool is_complex = Kokkos::Details::ArithTraits<T>::is_complex>
class ArithTraitsTester : public ArithTraitsTesterBase<T> {
public:
  static bool testFloatingPoint (std::ostream& out);
  static bool test (std::ostream& out);
};

// Specialization for real T.
template<class T>
class ArithTraitsTester<T, false> : public ArithTraitsTesterBase<T> {
public:
  static bool testFloatingPoint (std::ostream& out) {
    return ArithTraitsTesterBase<T>::testFloatingPoint (out);
  }
  static bool test (std::ostream& out) {
    using Kokkos::Details::ArithTraits;

    bool success = ArithTraitsTesterBase<T>::test (out);
    // Apparently, std::numeric_limits<T>::is_signed is true only for real numbers.
    if (ArithTraits<T>::is_signed != std::numeric_limits<T>::is_signed) {
      out << "ArithTraits<T>::is_signed != std::numeric_limits<T>::is_signed" << std::endl;
      success = false;
    }
    if (ArithTraits<T>::is_complex) {
      out << "ArithTraits<T>::is_complex is wrong" << std::endl;
      success = false;
    }
    return success;
  }
};

// Specialization for complex T.
template<class T>
class ArithTraitsTester<T, true> : public ArithTraitsTesterBase<T> {
public:
  static bool testFloatingPoint (std::ostream& out) {
    return ArithTraitsTesterBase<T>::testFloatingPoint (out);
  }
  static bool test (std::ostream& out) {
    using Kokkos::Details::ArithTraits;

    bool success = ArithTraitsTesterBase<T>::test (out);
    if (! ArithTraits<T>::is_complex) {
      out << "ArithTraits<T>::is_complex is wrong" << std::endl;
      success = false;
    }
    return success;
  }
};



#endif // KOKKOS_ARITHTRAITSTEST_HPP
