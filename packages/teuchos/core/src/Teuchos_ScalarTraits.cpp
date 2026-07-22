// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Assert.hpp"
#include <limits>

// Define this to throw exceptions when any Teuchos::ScalarTraits function
// encounters a NaN or an Inf.
//#define TEUCHOS_SCALAR_TRAITS_THROW_NAN_INF_ERR

#ifdef HAVE_TEUCHOSCORE_QUADMATH
namespace std {

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

istream&
operator>> (std::istream& in, __float128& x)
{
  std::string tmpStr;
  in >> tmpStr;
  // FIXME (mfh 10 Sep 2015) I don't think this routine does any error
  // checking.
  x = strtoflt128 (tmpStr.c_str (), NULL);
  return in;
}

} // namespace std
#endif // HAVE_TEUCHOSCORE_QUADMATH

void Teuchos::throwScalarTraitsNanInfError( const std::string &errMsg )
{
  (void)errMsg;
#ifdef TEUCHOS_SCALAR_TRAITS_THROW_NAN_INF_ERR
  TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, errMsg );
#endif
}

#ifdef HAVE_TEUCHOS_GNU_MP
gmp_randclass Teuchos::gmp_rng ( gmp_randinit_default );
#endif

#ifdef HAVE_TEUCHOS_QD
bool Teuchos::operator&&(const dd_real &a, const dd_real &b) {
  return !a.is_zero() && !b.is_zero();
}
bool Teuchos::operator&&(const qd_real &a, const qd_real &b) {
  return !a.is_zero() && !b.is_zero();
}
#endif

#ifndef __sun
// This is an intentional computation of NaN.
namespace Teuchos {
  const float  flt_nan = std::numeric_limits<float>::quiet_NaN();
  const double dbl_nan = std::numeric_limits<double>::quiet_NaN();
}
#endif
