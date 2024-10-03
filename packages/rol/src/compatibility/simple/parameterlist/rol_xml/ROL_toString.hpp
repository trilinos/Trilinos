// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TO_STRING_HPP
#define ROL_TO_STRING_HPP

namespace ROL {

/** \brief Default traits class for converting objects into strings.
 *
 * NOTE: This default implementation relies on operator<<(std::ostream&, ...)
 * being defined for the data type T.
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

} // end namespace ROL


#endif // ROL_TO_STRING_HPP
