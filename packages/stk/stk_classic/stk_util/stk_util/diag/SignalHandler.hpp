/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_UTIL_DIAG_SignalHandler_hpp
#define STK_UTIL_DIAG_SignalHandler_hpp

#include <map>

#include <stk_util/diag/String.hpp>
#include <stk_util/util/Callback.hpp>

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

  static bool check_signal_name(const sierra::String& signal);

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
  void add_handler(const String &signal_name, CallbackBase &callback);

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
  void remove_handler(const String &signal_name, CallbackBase &callback);

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
