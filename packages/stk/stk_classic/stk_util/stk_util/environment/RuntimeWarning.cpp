/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/environment/RuntimeWarning.hpp>

namespace stk_classic {

unsigned 
get_warning_count()
{
  return get_message_count(MSG_WARNING);
}


void
reset_warning_count()
{
  reset_message_count(MSG_WARNING);
}


void
set_max_warning_count(
  unsigned int	max_warnings)
{
  set_max_message_count(MSG_WARNING, max_warnings);
}


unsigned 
get_max_warning_count()
{
  return   get_max_message_count(MSG_WARNING);
}


void
report_warning(
  const char *          message,
  const MessageCode &   message_code)
{
  report_message(message, MSG_WARNING, message_code);
}


void
report_symmetric_warning(
  const char *          message,
  const MessageCode &   message_code)
{
  report_message(message, MSG_SYMMETRIC | MSG_WARNING, message_code);
}


void
report_deferred_warning(
  const char *          message,
  const char *          aggregate,
  const MessageCode &   message_code)
{
  add_deferred_message(MSG_WARNING, message_code.m_id, message_code.m_throttle.m_cutoff, message_code.m_throttle.m_group, message, aggregate);
}


RuntimeWarningAdHoc::RuntimeWarningAdHoc(
  MessageCode &   message_code)
  : m_messageCode(message_code)
{}


RuntimeWarningAdHoc::~RuntimeWarningAdHoc()
{
  try {
    report_warning(message.str().c_str(), m_messageCode);
  }
  catch (std::exception &)
  {}
}


RuntimeWarningSymmetric::RuntimeWarningSymmetric(
  MessageCode &   message_code)
  : m_messageCode(message_code)
{}


RuntimeWarningSymmetric::~RuntimeWarningSymmetric()
{
  try {
    report_symmetric_warning(message.str().c_str(), m_messageCode);
  }
  catch (std::exception &)
  {}
}


RuntimeWarningDeferred::RuntimeWarningDeferred(
  const MessageCode &   message_code)
  : m_messageCode(message_code)
{}


RuntimeWarningDeferred::~RuntimeWarningDeferred()
{
  try {
    report_deferred_warning(message.str().c_str(), aggregate.str().c_str(), m_messageCode);
  }
  catch (std::exception &)
  {}
}

} // namespace stk_classic
