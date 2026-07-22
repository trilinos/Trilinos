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

#include <dlfcn.h>
#include <boost/stacktrace.hpp>
#include <boost/exception/all.hpp>

/* \file  ROL_StackTrace.hpp
 * \brief Defines ROL_TEST_FOR_EXCEPTION using boost::stacktrace
 *
 * Minimally requires boost::stacktrace, boost::config, and boost::container_hash
 */


namespace ROL {

using traced = boost::error_info<struct tag_stacktrace, boost::stacktrace::stacktrace>;

template <class E>
void throw_with_trace(const E& e) {
    throw boost::enable_error_info(e)
        << traced(boost::stacktrace::stacktrace());
}

} // namespace ROL

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
    ROL::throw_with_trace<Exception>(Exception(omsgstr)); \
  } \
}


#endif // ROL_STACKTRACE_HPP

