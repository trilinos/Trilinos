/*------------------------------------------------------------------------*/
/*                 Copyright (c) 2013, Sandia Corporation.
/*                 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*                 the U.S. Governement retains certain rights in this software.
/*                 
/*                 Redistribution and use in source and binary forms, with or without
/*                 modification, are permitted provided that the following conditions are
/*                 met:
/*                 
/*                     * Redistributions of source code must retain the above copyright
/*                       notice, this list of conditions and the following disclaimer.
/*                 
/*                     * Redistributions in binary form must reproduce the above
/*                       copyright notice, this list of conditions and the following
/*                       disclaimer in the documentation and/or other materials provided
/*                       with the distribution.
/*                 
/*                     * Neither the name of Sandia Corporation nor the names of its
/*                       contributors may be used to endorse or promote products derived
/*                       from this software without specific prior written permission.
/*                 
/*                 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*                 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*                 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*                 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*                 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*                 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*                 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*                 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*                 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*                 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*                 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*                 
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/environment/ReportHandler.hpp>
#include <iostream>                     // for cout
#include <sstream>                      // for operator<<, basic_ostream, etc
#include <stdexcept>                    // for runtime_error, logic_error, etc


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
    "Error occured at: " + source_relative_path(location) + "\n" +
    error_msg);
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
    expr_msg +
    "Error occured at: " + source_relative_path(location) + "\n" +
    error_msg);
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
  /* %TRACE[ON]% */  /* %TRACE% */
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

} // namespace stk
