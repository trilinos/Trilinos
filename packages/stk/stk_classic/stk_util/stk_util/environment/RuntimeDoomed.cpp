/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <vector>
#include <stdexcept>

#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/environment/RuntimeDoomed.hpp>

namespace stk_classic {

unsigned int
get_doomed_count()
{
  return get_message_count(MSG_DOOMED);
}


void
reset_doomed_count()
{
  reset_message_count(MSG_DOOMED);
}


void
set_max_doomed_count(
  unsigned int	max_messages)
{
  set_max_message_count(MSG_DOOMED, max_messages);
}


unsigned
get_max_doomed_count()
{
  return get_max_message_count(MSG_DOOMED);
}


void
report_doomed(
  const char *          message,
  const MessageCode &   message_code)
{
  report_message(message, MSG_DOOMED, message_code);
}


void
report_symmetric_doomed(
  const char *          message,
  const MessageCode &   message_code)
{
  report_message(message, MSG_SYMMETRIC | MSG_DOOMED, message_code);
}
 

void
report_deferred_doomed(
  const char *          message,
  const char *          aggregate,
  const MessageCode &   message_code)
{
  add_deferred_message(MSG_DOOMED, message_code.m_id, message_code.m_throttle.m_cutoff, message_code.m_throttle.m_group, message, aggregate);
}


RuntimeDoomedAdHoc::RuntimeDoomedAdHoc(
  const MessageCode & message_code)
  : m_messageCode(message_code)
{}


RuntimeDoomedAdHoc::~RuntimeDoomedAdHoc()
{
  try {
    report_doomed(message.str().c_str(), m_messageCode);
  }
  catch (std::exception &)
  {}
}


RuntimeDoomedSymmetric::RuntimeDoomedSymmetric(
  const MessageCode & message_code)
  : m_messageCode(message_code)
{}


RuntimeDoomedSymmetric::~RuntimeDoomedSymmetric()
{
  try {
    report_symmetric_doomed(message.str().c_str(), m_messageCode);
  }
  catch (std::exception &)
  {}
}


RuntimeDoomedDeferred::RuntimeDoomedDeferred(
  const MessageCode & message_code)
  : m_messageCode(message_code)
{}


RuntimeDoomedDeferred::~RuntimeDoomedDeferred()
{
  try {
    report_deferred_doomed(message.str().c_str(), aggregate.str().c_str(), m_messageCode);
  }
  catch (std::exception &)
  {}
}

} // namespace stk_classic

