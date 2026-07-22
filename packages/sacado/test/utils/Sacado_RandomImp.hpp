// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstdlib>

template <typename ScalarT>
Sacado::Random<ScalarT>::
Random() :
  a(0.0),
  b(1.0),
  seed(static_cast<ScalarT>(rand()))
{
  // rand() can return 0 or 2147483647, so adjust seed if that happens
  if ((seed == 0.0) || (seed == 2147483647.0))
    seed = 1.0;
}

template <typename ScalarT>
Sacado::Random<ScalarT>::
Random(ScalarT a_, ScalarT b_) :
  a(a_),
  b(b_),
  seed(static_cast<ScalarT>(rand()))
{
  // rand() can return 0 or 2147483647, so adjust seed if that happens
  if ((seed == 0.0) || (seed == 2147483647.0))
    seed = 1.0;
}

template <typename ScalarT>
Sacado::Random<ScalarT>::
Random(ScalarT a_, ScalarT b_, int s) :
  a(a_),
  b(b_),
  seed(0.0)
{
  setSeed(s);
}

template <typename ScalarT>
Sacado::Random<ScalarT>::
~Random()
{
}

template <typename ScalarT>
void
Sacado::Random<ScalarT>::
setSeed(int s) {
  int ss = checkSeed("setSeed", s);
  srand(ss);
  seed = static_cast<ScalarT>(s);
}

template <typename ScalarT>
ScalarT
Sacado::Random<ScalarT>::
number() {
  const ScalarT A = 16807.0;
  const ScalarT bigInt = 2147483647.0;
      
  seed = std::fmod(A*seed, bigInt);
  return (b-a)*(seed/bigInt) + a;
}

template <typename ScalarT>
int
Sacado::Random<ScalarT>::
checkSeed(const std::string& func, int s) {
  if ((s < 1) || (s > 2147483646)) {
    std::cerr << "Error in Sacado::Random::" << s << "():  " 
	      << "supplied seed " 
	      << s << " is not an integer between 1 and 2147483646." 
	      << std::endl << "Using a seed of 1 instead." << std::endl;
    return 1;
  }
  else
    return s;
}

#ifdef HAVE_SACADO_COMPLEX

template <typename T>
Sacado::Random< std::complex<T> >::
Random() :
  rand_real(0.0, 1.0),
  rand_imag(0.0, 1.0)
{
}

template <typename T>
Sacado::Random< std::complex<T> >::
Random(const std::complex<T>& a, const std::complex<T>& b) :
  rand_real(a.real(), b.real()),
  rand_imag(a.imag(), b.imag())
{
}

template <typename T>
Sacado::Random< std::complex<T> >::
Random(const std::complex<T>& a, const std::complex<T>& b, int s) :
  rand_real(a.real(), b.real(), s),
  rand_imag(a.imag(), b.imag(), s+1)
{
}

template <typename T>
Sacado::Random< std::complex<T> >::
~Random()
{
}

template <typename T>
void
Sacado::Random< std::complex<T> >::
setSeed(int s) {
  rand_real.setSeed(s);
  rand_imag.setSeed(s+1);
}

template <typename T>
std::complex<T>
Sacado::Random< std::complex<T> >::
number() {
  return std::complex<T>(rand_real.number(), rand_imag.number());
}

#endif // HAVE_SACADO_COMPLEX
