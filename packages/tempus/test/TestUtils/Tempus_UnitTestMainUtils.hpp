//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_UNIT_TEST_MAIN_UTILS_HPP
#define TEMPUS_UNIT_TEST_MAIN_UTILS_HPP

#if defined(__linux__) && defined(__GNUC__) && !defined(__INTEL_COMPILER)
#include <fenv.h>
#elif defined(__APPLE__) && defined(__GNUC__) && defined(__SSE__)
#include <xmmintrin.h>
#endif

namespace Tempus_Test {

/** \brief Enable Floating Point Exceptions.
 *
 *  This enables/disables float point exceptions, Divide-by-Zero,
 *  Overflow, and Invalid for Linux platforms with GCC and
 *  Mac OS with GCC.  Otherwise, floating point exceptions are
 *  not used.
 */
void enableFPE(bool enableFPE)
{
#if defined(__APPLE__) && defined(__GNUC__) && defined(__SSE__)
  static int eMask = _MM_GET_EXCEPTION_MASK();
#endif

  if (enableFPE) {
#if defined(__linux__) && defined(__GNUC__) && !defined(__INTEL_COMPILER)
    feenableexcept(FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID);
#elif defined(__APPLE__) && defined(__GNUC__) && defined(__SSE__)
    eMask = _MM_GET_EXCEPTION_MASK();  // Save current eMask so we can disable.
    _MM_SET_EXCEPTION_MASK(eMask & ~_MM_MASK_DIV_ZERO & ~_MM_MASK_OVERFLOW &
                           ~_MM_MASK_INVALID);
#endif
  }
  else {
#if defined(__linux__) && defined(__GNUC__) && !defined(__INTEL_COMPILER)
    fedisableexcept(FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID);
#elif defined(__APPLE__) && defined(__GNUC__) && defined(__SSE__)
    _MM_SET_EXCEPTION_MASK(eMask);
#endif
  }
}

}  // namespace Tempus_Test

#endif  // TEMPUS_UNIT_TEST_MAIN_UTILS_HPP
