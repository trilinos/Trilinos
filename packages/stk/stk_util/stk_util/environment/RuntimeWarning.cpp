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

#include "stk_util/environment/RuntimeWarning.hpp"
#include "stk_util/environment/RuntimeMessage.hpp"  // for MSG_WARNING, report_message, MessageCode
#include <exception>                                // for exception
#include <string>                                   // for basic_string

namespace stk {

unsigned 
get_warning_count()
{
  return get_message_count(MSG_WARNING);
}

unsigned 
get_warning_printed_count()
{
  return get_message_printed_count(MSG_WARNING);
}

unsigned 
get_warning_printed_count(const MessageCode &messageCode)
{
  return get_message_printed_count(messageCode);
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

} // namespace stk
