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

/* \file  ROL_StackTrace.hpp
 * \brief Defines the "RELEASE" version of ROL_TEST_FOR_EXCEPTION 
 *        which simply throws an exception if the test returns true,
 *        bit does not create a stack trace.
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
    const std::string &omsgstr = omsg.str(); \
    throw Exception(omsgstr); \
  }\
}

#endif // ROL_STACKTRACE_HPP

