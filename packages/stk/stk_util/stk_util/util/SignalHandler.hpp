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

#ifndef STK_UTIL_DIAG_SignalHandler_hpp
#define STK_UTIL_DIAG_SignalHandler_hpp

#include "stk_util/util/Callback.hpp"  // for CallbackBase
#include <map>                         // for multimap, multimap<>::value_compare
#include <string>                      // for string

struct sigaction;

/**
 * @file
 *
 * The signal handler
 */

namespace sierra {

/**
 * @brief Class <b>SignalHandler</b> ...
 *
 */
class SignalHandler
{
public:
  /**
   * @brief Member function <b>instance</b> ...
   *
   * @return a <b>Handler</b> ...
   */
  static SignalHandler &instance();

  static bool check_signal_name(const std::string &signal);

  /**
   * @brief Member function <b>handle_signal</b> ...
   *
   * @param signal	an <b>int</b> variable ...
   */
  void handle_signal(int signal);

  /**
   * @brief Member function <b>add_handler</b> ...
   *
   * @param signal	an <b>int</b> variable ...
   * @param callback	a <b>CallbackBase</b> variable ...
   */
  void add_handler(int signal, CallbackBase &callback);

  /**
   * @brief Member function <b>add_handler</b> ...
   *
   * @param signal_name	a <b>String</b> variable ...
   * @param callback	a <b>CallbackBase</b> variable ...
   */
  void add_handler(const std::string &signal_name, CallbackBase &callback);

  /**
   * @brief Member function <b>remove_handler</b> ...
   *
   * @param signal	an <b>int</b> variable ...
   * @param callback	a <b>CallbackBase</b> variable ...
   */
  void remove_handler(int signal, CallbackBase &callback);

  /**
   * @brief Member function <b>remove_handler</b> ...
   *
   * @param signal_name	a <b>String</b> variable ...
   * @param callback	a <b>CallbackBase</b> variable ...
   */
  void remove_handler(const std::string &signal_name, CallbackBase &callback);

  /**
   * @brief Member function <b>remove_all_handlers</b> ...
   *
   */
  void remove_all_handlers();

private:
  typedef std::multimap<int, CallbackBase *> HandlerMap;
  typedef std::multimap<int, struct sigaction *> OldActionMap;

  HandlerMap		m_handlerMap;
  OldActionMap		m_oldActionMap;
};

} // namespace sierra

#endif // STK_UTIL_DIAG_SignalHandler_h
