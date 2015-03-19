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

/// \file float128.cpp
/// \brief Test for \c __float128 with Kokkos.
///
/// This file exercises the GCC extension \c __float128 with Kokkos.
/// This only works with GCC.  It also requires the "quadmath" TPL
/// (that provides \c __float128 library functions and macros) as well
/// as the following build flags:
///
/// <ul>
/// <li> <tt>-std=gnu++11</tt> (NOT c++11, else warnings) </li>
/// <li> <li>-fext-numeric-literals</tt> (else build errors) </li>
/// </ul>
///
/// Be sure to turn off "<tt>-ansi -pedantic</tt>", otherwise you'll
/// get more warnings (it's not ANSI C++ any more).

#include <iostream>
#include <sstream>
#include <quadmath.h>
#include <stdexcept>
#include <Kokkos_ArithTraits.hpp>
#include <Kokkos_Core.hpp>

namespace Kokkos {
namespace Details {

template<>
class ArithTraits<__float128> {
public:
  typedef __float128 val_type;
  typedef val_type mag_type;

  static const bool is_specialized = true;
  static const bool is_signed = true;
  static const bool is_integer = false;
  static const bool is_exact = false;
  static const bool is_complex = false;

  static bool isInf (const __float128 x) {
    return isinfq (x);
  }
  static bool isNan (const __float128 x) {
    return isnanq (x);
  }
  static mag_type abs (const __float128 x) {
    return fabsq (x);
  }
  static __float128 zero () {
    return 0.0;
  }
  static __float128 one () {
    return 1.0;
  }
  static __float128 min () {
    return FLT128_MIN;
  }
  static __float128 max () {
    return FLT128_MAX;
  }
  static mag_type real (const __float128 x) {
    return x;
  }
  static mag_type imag (const __float128 /* x */) {
    return 0.0;
  }
  static __float128 conj (const __float128 x) {
    return x;
  }
  static __float128 pow (const __float128 x, const __float128 y) {
    return powq (x, y);
  }
  static __float128 sqrt (const __float128 x) {
    return sqrtq (x);
  }
  static __float128 log (const __float128 x) {
    return logq (x);
  }
  static __float128 log10 (const __float128 x) {
    return log10q (x);
  }
  static mag_type epsilon () {
    return FLT128_EPSILON;
  }

  // Backwards compatibility with Teuchos::ScalarTraits.
  typedef mag_type magnitudeType;
  typedef double halfPrecision;
  // Unfortunately, we can't rely on a standard __float256 type.
  typedef __float128 doublePrecision;

  static const bool isComplex = false;
  static const bool isOrdinal = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = true;
  static bool isnaninf (const __float128 x) {
    return isNan (x) || isInf (x);
  }
  static magnitudeType magnitude (const __float128 x) {
    return abs (x);
  }
  static __float128 conjugate (const __float128 x) {
    return conj (x);
  }
  static std::string name () {
    return "__float128";
  }
  static __float128 squareroot (const __float128 x) {
    return sqrt (x);
  }
  static __float128 nan () {
    return strtoflt128 ("NAN()", NULL); // ???
  }
  static mag_type eps () {
    return epsilon ();
  }
  static mag_type sfmin () {
    return FLT128_MIN; // ???
  }
  static int base () {
    return 2;
  }
  static mag_type prec () {
    return eps () * static_cast<mag_type> (base ());
  }
  static int t () {
    return FLT_MANT_DIG;
  }
  static mag_type rnd () {
    return 1.0;
  }
  static int emin () {
    return FLT128_MIN_EXP;
  }
  static mag_type rmin () {
    return FLT128_MIN; // ??? // should be base^(emin-1)
  }
  static int emax () {
    return FLT128_MAX_EXP;
  }
  static mag_type rmax () {
    return FLT128_MAX; // ??? // should be (base^emax)*(1-eps)
  }
};

} // namespace Kokkos
} // namespace Details

using std::cout;
using std::endl;

// GCC / libquadmath doesn't implement an std::ostream operator<< for
// __float128, so we have to write our own.  At least libquadmath
// provides a printing function specifically for __float128.
std::ostream&
operator<< (std::ostream& out, const __float128& x)
{
  const size_t bufSize = 128;
  char buf[128];

  const int numCharPrinted = quadmath_snprintf (buf, bufSize, "%.30Qe", x);
  if (static_cast<size_t> (numCharPrinted) >= bufSize) {
    std::ostringstream os;
    os << "Failed to print __float128 value: buffer has " << bufSize
       << " characters, but quadmath_snprintf wanted " << numCharPrinted
       << " characters!";
    throw std::runtime_error (os.str ());
  }
  out << buf;
  return out;
}

int
main (int argc, char* argv[])
{
  Kokkos::initialize (argc, argv);

  __float128 x = 1.0;
  __float128 y = strtoflt128 ("1.111112222233333", NULL);
  __float128 z = strtoflt128 ("1.111112222233333444445555566666", NULL);

  // Make sure that all the digits print.
  cout << "x = " << x << endl
       << "y = " << y << endl
       << "z = " << z << endl
       << "(double) z = " << static_cast<double> (z) << endl
       << "z - (double) z = " << (z - static_cast<__float128> (static_cast<double> (z))) << endl;

  // Create a Kokkos::View on the host (CUDA doesn't work yet, since
  // __float128 is a GCC extension not available in CUDA).
  Kokkos::View<__float128*, Kokkos::HostSpace> view ("view", 20);

  // Increment the first entry, nonatomically.
  view(0)++;
  cout << "view(0) after increment = " << view(0) << endl;

  // Increment the first entry, atomically.
  Kokkos::atomic_add (&view(0), x);
  cout << "view(0) after atomic_add (x) = " << view(0) << endl;

  // Assign to the first entry, atomically.
  Kokkos::atomic_assign (&view(0), z);
  cout << "view(0) after atomic_assign (z) = " << view(0) << endl;

  cout << "Test PASSED" << endl;

  Kokkos::finalize ();
  return 0;
}
