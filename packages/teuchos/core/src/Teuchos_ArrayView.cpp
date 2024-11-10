// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ArrayView.hpp>

namespace Teuchos {

template<>
TEUCHOSCORE_LIB_DLL_EXPORT std::string
ArrayView<float>::toString() const
{
  using Teuchos::as;
  std::ostringstream ss;

  debug_assert_valid_ptr();

  ss.setf (std::ios::scientific);
  // 8 = round(23 * log10(2)) + 1.  That's one decimal digit more
  // than the binary precision justifies, which should be plenty.
  // Guy Steele et al. have a better algorithm for floating-point
  // I/O, but using a lot of digits is the lazy approach.
  ss.precision (8);
  ss << "{";
  for (size_type i = 0; i < size (); ++i) {
    ss << operator[] (i);
    if (i + 1 < size ()) {
      ss << ", ";
    }
  }
  ss << "}";
  return ss.str ();
}

template<>
TEUCHOSCORE_LIB_DLL_EXPORT std::string
ArrayView<const float>::toString() const
{
  using Teuchos::as;
  std::ostringstream ss;

  debug_assert_valid_ptr();

  ss.setf (std::ios::scientific);
  // 8 = round(23 * log10(2)) + 1.  That's one decimal digit more
  // than the binary precision justifies, which should be plenty.
  // Guy Steele et al. have a better algorithm for floating-point
  // I/O, but using a lot of digits is the lazy approach.
  ss.precision (8);
  ss << "{";
  for (size_type i = 0; i < size (); ++i) {
    ss << operator[] (i);
    if (i + 1 < size ()) {
      ss << ", ";
    }
  }
  ss << "}";
  return ss.str ();
}

template<>
TEUCHOSCORE_LIB_DLL_EXPORT std::string
ArrayView<double>::toString() const
{
  using Teuchos::as;
  std::ostringstream ss;

  debug_assert_valid_ptr();

  ss.setf (std::ios::scientific);
  // 17 = round(52 * log10(2)) + 1.  That's one decimal digit more
  // than the binary precision justifies, which should be plenty.  Guy
  // Steele et al. have a better algorithm for floating-point I/O, but
  // using a lot of digits is the lazy approach.
  ss.precision (17);
  ss << "{";
  for (size_type i = 0; i < size (); ++i) {
    ss << operator[] (i);
    if (i + 1 < size ()) {
      ss << ", ";
    }
  }
  ss << "}";
  return ss.str ();
}

template<>
TEUCHOSCORE_LIB_DLL_EXPORT std::string
ArrayView<const double>::toString() const
{
  using Teuchos::as;
  std::ostringstream ss;

  debug_assert_valid_ptr();

  ss.setf (std::ios::scientific);
  // 17 = round(52 * log10(2)) + 1.  That's one decimal digit more
  // than the binary precision justifies, which should be plenty.  Guy
  // Steele et al. have a better algorithm for floating-point I/O, but
  // using a lot of digits is the lazy approach.
  ss.precision (17);
  ss << "{";
  for (size_type i = 0; i < size (); ++i) {
    ss << operator[] (i);
    if (i + 1 < size ()) {
      ss << ", ";
    }
  }
  ss << "}";
  return ss.str ();
}

} // namespace Teuchos
