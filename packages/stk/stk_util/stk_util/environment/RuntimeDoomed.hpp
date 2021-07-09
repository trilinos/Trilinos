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

#ifndef STK_UTIL_ENVIRONMENT_RUNTIMEDOOMED_HPP
#define STK_UTIL_ENVIRONMENT_RUNTIMEDOOMED_HPP

#include "stk_util/environment/RuntimeMessage.hpp"  // for MessageCode, MessageCode::s_defaultMe...
#include <sstream>                                  // for ostringstream, ios_base, ostream


namespace stk {

typedef std::ostream &(*OStreamFunctionPtr)(std::ostream &);
typedef std::ios_base &(*IOSBaseFunctionPtr)(std::ios_base &);

///
/// @addtogroup runtime_message_detail
/// @{
///

/**
 * @brief Function <b>get_doomed_count</b> returns the accumulated doomed count.
 *
 */
unsigned get_doomed_count();

/**
 * @brief Function <b>is_doomed</b> returns true if get_doomed_count() > 0.
 *
 * @return			a <b>bool</b> value of true if get_doomed_count() > 0.
 */
inline bool is_doomed() {
  return get_doomed_count() > 0;
}

/**
 * @brief Function <b>reset_doomed_count</b> sets the accumulated doomed count to zero.
 *
 */
void reset_doomed_count();

/**
 * @brief Function <b>set_max_messages</b> sets the maximum number of doomed before no more doomed
 * will be displayed.
 */
void set_max_doomed_count(unsigned int max_messages);

/**
 * @brief Function <b>report_symmetric_doomed</b> sends a doomed message to the reporter.
 *
 * Since the issue causing the doomed error may occur thousands of times during a run, it is
 * desirable to throttle the number of times that a message is displayed.  If you desire to limit
 * the number of times a message is displayed, obtain a unique id via the get_next_message_id()
 * function.  This function allows you to assign the maximum display count.  When the count of
 * messages reported with this message id exceeds the display count, it is no longer displayed.  The
 * default id zero (0) is assigned which has an extremely large display count.
 *
 * @param message		a <b>char</b> const pointer to a message that is to be displayed.
 *
 * @param message_code		a <b>MessageCode</b> const pointer to a cross-processor unique
 *                              identification of this message.
 *
 */
void report_doomed(const char *message, const MessageCode &message_code = MessageCode::s_defaultMessageCode);

/**
 * @brief Function <b>report_symmetric_doomed</b> sends a doomed message to the reporter.
 *
 * Since the issue causing the doomed error may occur thousands of times during a run, it is
 * desirable to throttle the number of times that a message is displayed.  If you desire to limit
 * the number of times a message is displayed, obtain a unique id via the get_next_message_id()
 * function.  This function allows you to assign the maximum display count.  When the count of
 * messages reported with this message id exceeds the display count, it is no longer displayed.  The
 * default id zero (0) is assigned which has an extremely large display count.
 *
 * @param message		a <b>char</b> const pointer to a message that is to be displayed.
 *
 * @param message_code		a <b>MessageCode</b> const pointer to a cross-processor unique
 *                              identification of this message.
 *
 */
void report_symmetric_doomed(const char *message, const MessageCode &message_code = MessageCode::s_defaultMessageCode);

void report_deferred_doomed(const char *message, const char *aggregate, const MessageCode &message_code);

/**
 * @brief Class <b>RuntimeDoomedAdHoc</b> reports an ad hoc doomed message to the report system.
 *
 * For example:
 *
 * <PRE>
 *     if (adhoc_runtime_doomed_condition)
 *       RuntimeDoomedAdHoc() << "My useful message about " << some_data;
 *
 *     if (adhoc_runtime_doomed_condition) {
 *       static MessageCode mc;
 *       RuntimeDoomedAdHoc(mc) << "My useful message about " << some_data;
 *     }
 * </PRE>
 */
class RuntimeDoomedAdHoc
{
public:
  /**
   * @brief Creates a new <b>RuntimeDoomedAdHoc</b> instance, setting the message code.
   *
   * @param message_code	an <b>MessageCode</b> const reference to the message code associated
   *                            with this message.
   *
   */
  explicit RuntimeDoomedAdHoc(const MessageCode &message_code = MessageCode::s_defaultMessageCode);

  /**
   * @brief Destroys a <b>RuntimeDoomedAdHoc</b> instance.
   *
   * The message is displayed by calling the report_doomed() function.  However, if the count of
   * remaining messages for this message id is zero, the message is not displayed.
   *
   */
  ~RuntimeDoomedAdHoc();

private:
  /**
   * @brief Make copy of <b>RuntimeDoomedAdHoc</b> invalid.
   *
   */
  RuntimeDoomedAdHoc(const RuntimeDoomedAdHoc &);

  /**
   * @brief Make assignment of <b>RuntimeDoomedAdHoc</b> invalid.
   *
   */
  RuntimeDoomedAdHoc &operator=(const RuntimeDoomedAdHoc &);

public:
  /**
   * @brief Member function <b>operator&lt;&lt;</b> passes the output stream manipulator to the
   * output stream.
   *
   * @return			a <b>RuntimeDoomedAdHoc</b> reference to this object
   */
  RuntimeDoomedAdHoc &operator<<(OStreamFunctionPtr f) {
    f(message);
    return *this;
  }

  /**
   * @brief Member function <b>operator&lt;&lt;</b> passes the ios_base manipulator to the output
   * stream.
   *
   * @return			a <b>RuntimeDoomedAdHoc</b> reference to this object
   */
  RuntimeDoomedAdHoc &operator<<(IOSBaseFunctionPtr f) {
    f(message);
    return *this;
  }

  /**
   * @brief Member function <b>operator&lt;&lt;</b> passes any data type to the
   * exception string class for conversion to a string.
   *
   * @param t			a <b>T</b> const reference that is to be converted
   *				to a string.
   *
   * @return			a <b>RuntimeDoomedAdHoc</b> reference to this object;
   */
  template <class T>
  RuntimeDoomedAdHoc &operator<<(const T &t) {
    message << t;
    return *this;
  }

public:
  std::ostringstream    message;                ///< Stream to receive message content

private:
  const MessageCode     m_messageCode;          ///< Message id and uninitialized throttle
};


/**
 * @brief Class <b>RuntimeDoomedSymmetric</b> reports a fatal error message to the report system.
 *
 * For example:
 *
 * <PRE>
 *     if (symmetric_runtime_doomed_condition)
 *       RuntimeDoomedSymmetric() << "My useful message about " << some_data;
 *
 *     if (symmetric_runtime_doomed_condition) {
 *       static MessageCode mc;
 *       RuntimeDoomedSymmetric(mc) << "My useful message about " << some_data;
 *     }
 * </PRE>
 */
class RuntimeDoomedSymmetric
{
public:
  /**
   * @brief Creates a new <b>RuntimeDoomedSymmetric</b> instance, setting the message code.
   *
   * @param message_code	an <b>MessageCode</b> const reference to the message code associated
   *                            with this message.
   *
   */
  explicit RuntimeDoomedSymmetric(const MessageCode &message_code = MessageCode::s_defaultMessageCode);

  /**
   * @brief Destroys a <b>RuntimeDoomedSymmetric</b> instance.
   *
   * The message is displayed by calling the report_doomed() function.  However, if the count of
   * remaining messages for this message id is zero, the message is not displayed.
   *
   */
  ~RuntimeDoomedSymmetric();

private:
  /**
   * @brief Make copy of <b>RuntimeDoomedSymmetric</b> invalid.
   *
   */
  RuntimeDoomedSymmetric(const RuntimeDoomedSymmetric &);

  /**
   * @brief Make assignment of <b>RuntimeDoomedSymmetric</b> invalid.
   *
   */
  RuntimeDoomedSymmetric &operator=(const RuntimeDoomedSymmetric &);

public:
  /**
   * @brief Member function <b>operator&lt;&lt;</b> passes the ios_base manipulator
   * to the output stream.
   *
   * @return			a <b>RuntimeDoomedSymmetric</b> reference to this object
   */
  RuntimeDoomedSymmetric &operator<<(OStreamFunctionPtr f) {
    f(message);
    return *this;
  }

  /**
   * @brief Member function <b>operator&lt;&lt;</b> passes the ios_base manipulator
   * to the output stream.
   *
   * @return			a <b>RuntimeDoomedSymmetric</b> reference to this object
   */
  RuntimeDoomedSymmetric &operator<<(IOSBaseFunctionPtr f) {
    f(message);
    return *this;
  }

  /**
   * @brief Member function <b>operator&lt;&lt;</b> passes any data type to the
   * exception string class for conversion to a string.
   *
   * @param t			a <b>T</b> const reference that is to be converted
   *				to a string.
   *
   * @return			a <b>RuntimeDoomedSymmetric</b> reference to this object;
   */
  template <class T>
  RuntimeDoomedSymmetric &operator<<(const T &t) {
    message << t;
    return *this;
  }

public:
  std::ostringstream    message;                ///< Stream to receive message content

private:
  const MessageCode     m_messageCode;          ///< Message id and uninitialized throttle
};


/**
 * @brief Class <b>RuntimeDoomedDeferred</b> reports a deferred fatal error message to the report
 * system.
 *
 * For example:
 *
 * <PRE>
 *     if (deferred_runtime_doomed_condition) {
 *       static MessageCode mc;
 *       RuntimeDoomedDeferred(mc) << "My useful message about " << some_data;
 *     }
 *
 *     if (deferred_runtime_doomed_condition) {
 *       static MessageCode mc;
 *       RuntimeDoomedDeferred x;
 *       x << "My useful message about " << some_data;
 *       x.aggregate << proc_specific_data;
 *     }
 * </PRE>
 */
class RuntimeDoomedDeferred
{
public:
  /**
   * @brief Creates a new <b>RuntimeDoomedDeferred</b> instance, setting the message code.
   *
   * @param message_code	an <b>MessageCode</b> const reference to the message code associated
   *                            with this message.
   *
   */
  explicit RuntimeDoomedDeferred(const MessageCode &message_code);

  /**
   * @brief Destroys a <b>RuntimeDoomed</b> instance.
   *
   * The message is displayed by calling the add_deferred_message() function.
   *
   */
  ~RuntimeDoomedDeferred();

private:
  /**
   * @brief Make copy of <b>RuntimeDoomedDeferred</b> invalid.
   *
   */
  RuntimeDoomedDeferred(const RuntimeDoomedDeferred &);

  /**
   * @brief Make assignment of <b>RuntimeDoomedDeferred</b> invalid.
   *
   */
  RuntimeDoomedDeferred &operator=(const RuntimeDoomedDeferred &);

public:
  /**
   * @brief Member function <b>operator&lt;&lt;</b> passes the ios_base manipulator
   * to the output stream.
   *
   * @return			a <b>RuntimeDoomedDeferred</b> reference to this object
   */
  RuntimeDoomedDeferred &operator<<(OStreamFunctionPtr f) {
    f(message);
    return *this;
  }

  /**
   * @brief Member function <b>operator&lt;&lt;</b> passes the ios_base manipulator
   * to the output stream.
   *
   * @return			a <b>RuntimeDoomedDeferred</b> reference to this object
   */
  RuntimeDoomedDeferred &operator<<(IOSBaseFunctionPtr f) {
    f(message);
    return *this;
  }

  /**
   * @brief Member function <b>operator&lt;&lt;</b> passes any data type to the
   * exception string class for conversion to a string.
   *
   * @param t			a <b>T</b> const reference that is to be converted
   *				to a string.
   *
   * @return			a <b>RuntimeDoomedDeferred</b> reference to this object;
   */
  template <class T>
  RuntimeDoomedDeferred &operator<<(const T &t) {
    message << t;
    return *this;
  }

public:
  std::ostringstream    message;                ///< Stream to receive message header content
  std::ostringstream    aggregate;              ///< Stream to receive message aggregate content

private:
  const MessageCode     m_messageCode;          ///< Message id and uninitialized throttle
};

using RuntimeDoomed = stk::RuntimeDoomedAdHoc;

///
/// @}
///

} // namespace stk

#endif // STK_UTIL_ENVIRONMENT_RUNTIMEDOOMED_HPP
