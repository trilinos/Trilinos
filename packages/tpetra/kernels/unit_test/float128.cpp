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


// GCC / libquadmath doesn't implement an std::ostream operator<< for
// __float128, so we have to write our own.  At least libquadmath
// provides a printing function specifically for __float128.
//
// FIXME (mfh 19 Mar 2015) This will break if users have already
// defined their own operator<< in the global namespace.  Note that we
// already have implemented this in Teuchos_ScalarTraits.hpp, which is
// why we enclose this in an anonymous namespace.
namespace {
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
}

using std::cout;
using std::endl;

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
