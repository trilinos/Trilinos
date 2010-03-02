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

void
default_report_handler(
  const char *		message,
  int                   type)
{
  std::cout << "Message type " << type << ": " << message << std::endl;
}

REH
s_reportHandler = &default_report_handler;


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

} // namespace stk
