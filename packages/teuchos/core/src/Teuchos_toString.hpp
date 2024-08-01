// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_TO_STRING_HPP
#define TEUCHOS_TO_STRING_HPP

#include "Teuchos_ConfigDefs.hpp"
#ifdef HAVE_TEUCHOSCORE_QUADMATH
#  include <quadmath.h> // __float128 functions
#endif // HAVE_TEUCHOSCORE_QUADMATH

namespace Teuchos {

/** \brief Default traits class for converting objects into strings.
 *
 * NOTE: This default implementation relies on operator<<(std::ostream&, ...)
 * being defined for the data type T.
 *
 * \ingroup teuchos_language_support_grp
 */
template<typename T>
class ToStringTraits {
public:
  static std::string toString( const T &t )
    {
      std::ostringstream oss;
      oss << t;
      return oss.str();
    }
};


/** \brief Utility function for returning a pretty string representation of
 * a object of type T.
 *
 * NOTE: This helper function simply returns ToStringTraits<T>::toString(t)
 * and the right way to speicalize the behavior is to specialize
 * ToStringTraits.
 *
 * \ingroup teuchos_language_support_grp
 */
template<typename T>
inline
std::string toString(const T& t)
{
  return ToStringTraits<T>::toString(t);
}


/** \brief Specialization for bool. */
template<>
class ToStringTraits<bool> {
public:
  static std::string toString( const bool &t )
    {
      if (t)
        return "true";
      return "false";
    }
};


/** \brief Specialization for std::string. */
template<>
class ToStringTraits<std::string> {
public:
  static std::string toString( const std::string &t )
    {
      return t;
    }
};

/** \brief Specialization for double. */
template<>
class ToStringTraits<double> {
public:
  static std::string toString (const double& t) {
    std::ostringstream os;
    os.setf (std::ios::scientific);
    // 17 = round(52 * log10(2)) + 1.  That's one decimal digit more
    // than the binary precision justifies, which should be plenty.
    // Guy Steele et al. have a better algorithm for floating-point
    // I/O, but using a lot of digits is the lazy approach.
    os.precision (17);
    os << t;
    return os.str();
  }
};

#ifdef HAVE_TEUCHOS_LONG_DOUBLE
/** \brief Specialization for long double. */
template<>
class ToStringTraits<long double> {
public:
  static std::string toString (const long double& t) {
    std::ostringstream os;
    os.setf (std::ios::scientific);
    // 26 = round(80 * log10(2)) + 1.  That's one decimal digit more
    // than the binary precision justifies, which should be plenty.
    // Guy Steele et al. have a better algorithm for floating-point
    // I/O, but using a lot of digits is the lazy approach.
    os.precision (26);
    os << t;
    return os.str();
  }
};
#endif

/** \brief Specialization for float. */
template<>
class ToStringTraits<float> {
public:
  static std::string toString (const float& t) {
    std::ostringstream os;
    os.setf (std::ios::scientific);
    // 8 = round(23 * log10(2)) + 1.  That's one decimal digit more
    // than the binary precision justifies, which should be plenty.
    // Guy Steele et al. have a better algorithm for floating-point
    // I/O, but using a lot of digits is the lazy approach.
    os.precision (8);
    os << t;
    return os.str();
  }
};


#ifdef HAVE_TEUCHOSCORE_QUADMATH
/// \brief Partial specialization for \c __float128.
///
/// \c __float128 is a GCC language extension.  It requires linking
/// with libquadmath.
template<>
class ToStringTraits<__float128> {
public:
  static std::string toString (const __float128& val)
  {
    // libquadmath doesn't implement operator<< (std::ostream&,
    // __float128), but it does have a print function.
    const size_t bufSize = 128;
    char buf[128];

    // FIXME (mfh 04 Sep 2015) We should test the returned int value
    // to make sure that it is < bufSize.  On the other hand, I do not
    // want to add more header dependencies to this file, and
    // TEUCHOS_TEST_FOR_EXCEPTION is not already included.
#if 0
    const int numCharPrinted = quadmath_snprintf (buf, bufSize, "%.30Qe", val);
    TEUCHOS_TEST_FOR_EXCEPTION
      (static_cast<size_t> (numCharPrinted) >= bufSize, std::runtime_error,
       "Teuchos::toString: Failed to print __float128 value: buffer has "
       << bufSize << " characters, but quadmath_snprintf wanted "
       << numCharPrinted << " characters!");
#else
    (void) quadmath_snprintf (buf, bufSize, "%.30Qe", val);
#endif
    return std::string (buf);
  }
};
#endif // HAVE_TEUCHOSCORE_QUADMATH


/// \brief Partial specialization for std::pair<T1, T2>.
///
/// \note This relies on operator<< working for T1 and T2.
template<typename T1, typename T2>
class ToStringTraits<std::pair<T1, T2> > {
public:
  static std::string toString (const std::pair<T1, T2>& t) {
    std::ostringstream oss;
    oss << "(" << t.first << "," << t.second << ")";
    return oss.str();
  }
};

} // end namespace Teuchos


#endif // TEUCHOS_TO_STRING_HPP

