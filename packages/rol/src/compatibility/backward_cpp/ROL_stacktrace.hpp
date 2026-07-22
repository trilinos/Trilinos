// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_STACKTRACE_HPP
#define ROL_STACKTRACE_HPP

#include <sstream>
#include "backward.hpp"

/* \file  ROL_StackTrace.hpp
 * \brief Defines ROL_TEST_FOR_EXCEPTION using backward-cpp
 *
 * https://github.com/bombela/backward-cpp
 * 
 */


#define ROL_TEST_FOR_EXCEPTION(throw_exception_test, Exception, msg) \
{ \
  const bool throw_exception = (throw_exception_test); \
  if(throw_exception) { \
    std::ostringstream omsg; \
    omsg \
      << __FILE__ << ":" << __LINE__ << ":\n\n" \
      << "\n\n" \
      << "Throw test that evaluated to true: "#throw_exception_test \
      << "\n\n" \
      << msg; \
    using namespace backward;  \
    StackTrace st; st.load_here(32); \
    Printer p;  \
    p.object = true; \
    p.color_mode = ColorMode::always; \
    p.address = true;  \
    p.print(st, omsg); \
    const std::string &omsgstr = omsg.str(); \
    throw Exception(omsgstr); \
   } \
}


#endif // ROL_STACKTRACE_HPP

