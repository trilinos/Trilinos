/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/environment/ReportHandler.hpp>

#include <iostream>
#include <stdexcept>

namespace stk {

namespace {

// Global variables used to store the current handlers. We do not want direct
// external access to this state, so we put these in an empty namespace.
ErrorHandler s_assert_handler = &default_assert_handler;
REH s_reportHandler = &default_report_handler;

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
  std::string error_msg = "";
  if (message.str() != "") {
    error_msg = std::string("Error message: ") + message.str() + "\n";
  }
  throw std::logic_error(
    std::string("Requirement( ") + expr + " ) FAILED\n" +
    "Error occured at: " + source_relative_path(location) + "\n" +
    error_msg);
}

ErrorHandler set_assert_handler(ErrorHandler error_handler)
{
  if (!error_handler)
    throw std::runtime_error("Cannot set assert handler to NULL");

  ErrorHandler prev_error_handler = s_assert_handler;
  s_assert_handler = error_handler;

  return prev_error_handler;
}

void handle_assert(const char* expr,
                   const std::string& location,
                   std::ostringstream& message)
{
  (*s_assert_handler)(expr, location, message);
}

} // namespace stk
