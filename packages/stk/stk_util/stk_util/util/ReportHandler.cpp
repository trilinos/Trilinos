// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include "stk_util/stk_config.h"

#include "stk_util/util/ReportHandler.hpp"
#include "Kokkos_Core.hpp"
#include <iostream>   // for cout
#include <sstream>    // for ostringstream, operator<<, basic_ostream, basic_ostream::operator<<
#include <stdexcept>  // for runtime_error, logic_error, invalid_argument

#ifdef STK_HAVE_BOOST
#include "boost/version.hpp"
#if BOOST_VERSION >= 106500
#ifndef STK_NO_BOOST_STACKTRACE
#define STK_HAVE_BOOST_STACKTRACE
#endif
#endif
#endif

#ifdef STK_HAVE_BOOST_STACKTRACE
#ifdef __INTEL_COMPILER
#include "boost/stacktrace.hpp"
#else
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wattributes"
#include "boost/stacktrace.hpp"
#pragma GCC diagnostic pop
#endif
#endif

namespace stk {
namespace {

// Global variables used to store the current handlers. We do not want direct
// external access to this state, so we put these in an empty namespace.
ErrorHandler s_assert_handler = &default_assert_handler;
ErrorHandler s_error_handler = &default_error_handler;
ErrorHandler s_invalid_arg_handler = &default_invalid_arg_handler;
REH s_reportHandler = &default_report_handler;

template<class EXCEPTION>
void default_handler_req(const char* expr,
                         const std::string& location,
                         std::ostringstream& message)
{
  std::string error_msg = "";
  if (message.str() != "") {
    error_msg = std::string("Error: ") + message.str() + "\n";
  }

  throw EXCEPTION(
    std::string("Requirement( ") + expr + " ) FAILED\n" +
    "Error occurred at: " + location + "\n" + error_msg);
}


template<class EXCEPTION>
void default_handler_exc(const char* expr,
                         const std::string& location,
                         std::ostringstream& message)
{
  std::string error_msg = "";
  if (message.str() != "") {
    error_msg = std::string("Error: ") + message.str() + "\n";
  }

  std::string expr_msg = "";
  if (expr != std::string("")) {
    expr_msg = std::string("Expr '") + expr + "' eval'd to true, throwing.\n";
  }

  throw EXCEPTION(
    expr_msg + "Error occurred at: " + location + "\n" + error_msg);
}

}

void clean_error_handler(const char* expr,
                         const std::string& location,
                         std::ostringstream& message)
{
  std::string error_msg = "";
  if (message.str() != "") {
    error_msg = std::string("Error: ") + message.str() + "\n";
  }

  throw std::logic_error(
    std::string("Requirement( ") + expr + " ) FAILED\n" +
    error_msg);
}


void
default_report_handler(
  const char *		message,
  int                   type)
{
  std::cout << "Message type " << type << ": " << message << std::endl;
}

void
report(
  const char *		message,
  int                   type)
{
    (*s_reportHandler)(message, type);
}


REH
set_report_handler(
  REH		        reh)
{
  if (!reh)
    throw std::runtime_error("Cannot set report handler to NULL");

  REH prev_reh = s_reportHandler;
  s_reportHandler = reh;

  return prev_reh;
}


std::string
source_relative_path(
  const std::string &	path)
{
  static const char *prefix[] = {"/src/", "/include/", "/Apps_", "/stk_"};

  for (unsigned int i = 0; i < sizeof(prefix)/sizeof(prefix[0]); ++i) {
    std::string::size_type j = path.rfind(prefix[i], path.length());
    if (j != std::string::npos) {
      j = path.rfind("/", j - 1);
      if (j != std::string::npos)
	return path.substr(j + 1, path.length());
      else
	return path;
    }
  }
  return path;
}

void default_assert_handler(const char* expr,
                            const std::string& location,
                            std::ostringstream& message)
{
  default_handler_req<std::logic_error>(expr, location, message);
}

void default_error_handler(const char* expr,
                           const std::string& location,
                           std::ostringstream& message)
{
  default_handler_exc<std::runtime_error>(expr, location, message);
}

void default_invalid_arg_handler(const char* expr,
                                 const std::string& location,
                                 std::ostringstream& message)
{
  default_handler_exc<std::invalid_argument>(expr, location, message);
}

ErrorHandler set_assert_handler(ErrorHandler handler)
{
  if (!handler)
    throw std::runtime_error("Cannot set assert handler to NULL");

  ErrorHandler prev_handler = s_assert_handler;
  s_assert_handler = handler;

  return prev_handler;
}

ErrorHandler set_error_handler(ErrorHandler handler)
{
  if (!handler)
    throw std::runtime_error("Cannot set error handler to NULL");

  ErrorHandler prev_handler = s_error_handler;
  s_error_handler = handler;

  return prev_handler;
}

ErrorHandler set_invalid_arg_handler(ErrorHandler handler)
{
  if (!handler)
    throw std::runtime_error("Cannot set invalid_arg handler to NULL");

  ErrorHandler prev_handler = s_invalid_arg_handler;
  s_invalid_arg_handler = handler;

  return prev_handler;
}

void handle_assert(const char* expr,
                   const std::string& location,
                   std::ostringstream& message)
{
  (*s_assert_handler)(expr, location, message);
}

void handle_error(const char* expr,
                  const std::string& location,
                  std::ostringstream& message)
{
  (*s_error_handler)(expr, location, message);
}

void handle_invalid_arg(const char* expr,
                        const std::string& location,
                        std::ostringstream& message)
{
  (*s_invalid_arg_handler)(expr, location, message);
}

#ifdef STK_HAVE_BOOST_STACKTRACE
std::ostream & output_stacktrace(std::ostream & os)
{
#ifdef __INTEL_COMPILER
#pragma warning( disable: 2196 )
#endif
  os << boost::stacktrace::stacktrace();
  return os;
}
#else
std::ostream & output_stacktrace(std::ostream & os)
{
  return os;
}
#endif

} // namespace stk

#ifndef STK_ENABLE_GPU_BUT_NO_RDC

STK_FUNCTION void ThrowMsgDevice(const char * message)
{
  Kokkos::abort(message);
}

STK_FUNCTION void ThrowErrorMsgDevice(const char * message)
{
  Kokkos::abort(message);
}

#endif

