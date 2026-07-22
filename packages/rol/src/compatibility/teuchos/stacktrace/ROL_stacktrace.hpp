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

#include "Teuchos_TestForException.hpp"

/* \file  ROL_StackTrace.hpp
 * \brief Defines ROL_TEST_FOR_EXCEPTION using the Teuchos implementation
 */

#define ROL_TEST_FOR_EXCEPTION(throw_exception_test, Exception, msg) \
  TEUCHOS_TEST_FOR_EXCEPTION(throw_exception_test, Exception, msg)

#endif // ROL_STACKTRACE_HPP

