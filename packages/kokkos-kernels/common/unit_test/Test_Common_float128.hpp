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
// #include <cstdlib>
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
std::ostream& operator<<(std::ostream& out, const __float128& x) {
  const size_t bufSize = 128;
  char buf[128];

  const int numCharPrinted = quadmath_snprintf(buf, bufSize, "%.30Qe", x);
  if (static_cast<size_t>(numCharPrinted) >= bufSize) {
    std::ostringstream os;
    os << "Failed to print __float128 value: buffer has " << bufSize << " characters, but quadmath_snprintf wanted "
       << numCharPrinted << " characters!";
    throw std::runtime_error(os.str());
  }
  out << buf;
  return out;
}
}  // namespace

using std::cout;
using std::endl;

void testfloat128() {
  bool success = true;
  __float128 x = 1.0;
  __float128 y = strtoflt128("1.111112222233333", NULL);
  __float128 z = strtoflt128("1.111112222233333444445555566666", NULL);

  // Make sure that all the digits print.
  cout << "x = " << x << endl
       << "y = " << y << endl
       << "z = " << z << endl
       << "(double) z = " << static_cast<double>(z) << endl
       << "z - (double) z = " << (z - static_cast<__float128>(static_cast<double>(z))) << endl;

  // FIXME (mfh 04 Sep 2015) The results of printing could depend on
  // the locale.  This works fine for the default locale on my system.
  {
    std::ostringstream os;
    os << x;
    if (os.str() != "1.000000000000000000000000000000e+00") {
      success = false;
      cout << "'_float128 x = 1.0' does not print correctly!  It prints as " << os.str() << "." << endl;
    }
  }
  {
    std::ostringstream os;
    os << y;
    if (os.str() != "1.111112222233333000000000000000e+00") {
      success = false;
      cout << "'__float128 y = strtoflt128 (\"1.111112222233333\", NULL);' "
              "does not print correctly!  It prints as "
           << os.str() << "." << endl;
    }
  }
  {
    std::ostringstream os;
    os << z;
    if (os.str() != "1.111112222233333444445555566666e+00") {
      success = false;
      cout << "'__float128 z = strtoflt128 "
              "(\"1.111112222233333444445555566666\", NULL);' "
              "does not print correctly!  It prints as "
           << os.str() << "." << endl;
    }
  }

  // Create a Kokkos::View on the host (CUDA doesn't work yet, since
  // __float128 is a GCC extension not available in CUDA).
  Kokkos::View<__float128*, Kokkos::HostSpace> view("view", 20);

  // Increment the first entry, nonatomically.
  view(0)++;
  cout << "view(0) after increment = " << view(0) << endl;
  if (view(0) != static_cast<__float128>(1.0)) {
    success = false;
  }

  // Increment the first entry, atomically.
  Kokkos::atomic_add(&view(0), x);
  cout << "view(0) after atomic_add (x) = " << view(0) << endl;
  if (view(0) != static_cast<__float128>(2.0)) {
    success = false;
  }

  // Assign to the first entry, atomically.
  Kokkos::atomic_assign(&view(0), z);
  cout << "view(0) after atomic_assign (z) = " << view(0) << endl;
  if (view(0) != z) {
    success = false;
  }
  EXPECT_TRUE((success));
}

TEST_F(TestCategory, common_float128) { testfloat128(); }
