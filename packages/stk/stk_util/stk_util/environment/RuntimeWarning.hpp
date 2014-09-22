// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#ifndef STK_UTIL_ENVIRONMENT_RUNTIMEWARNING_HPP
#define STK_UTIL_ENVIRONMENT_RUNTIMEWARNING_HPP

#include <sstream>

#include <stk_util/environment/RuntimeMessage.hpp>

namespace stk {

typedef std::ostream &(*OStreamFunctionPtr)(std::ostream &);
typedef std::ios_base &(*IOSBaseFunctionPtr)(std::ios_base &);

///
/// @addtogroup runtime_message_detail
/// @{
///

/**
 * @brief Function <b>get_warning_count</b> returns the accumulated warning count.
 *
 */
unsigned get_warning_count();

/**
 * @brief Function <b>reset_warning_count</b> sets the accumulated warning count to zero.
 *
 */
void reset_warning_count();

/**
 * @brief Function <b>set_max_messages</b> sets the maximum number of warning before no more warning
 * will be displayed.
 *
 * @param max_messages	an <b>int</b> variable ...
 *
 */
void set_max_warning_count(unsigned int max_messages);

/**
 * @brief Function <b>set_max_messages</b> sets the maximum number of warning
 * and doomed messages displayed before the message is thrown as a RuntimeError exception.
 *
 */
unsigned get_max_warning_count();

/**
 * @brief Member function <b>report_warning</b> ...
 *
 * @param message		a <b>char</b> const pointer ...
 *
 * @param message_code		a <b>MessageCode</b> const ...
 *
 */
void report_warning(const char *message, const MessageCode &message_code = MessageCode::s_defaultMessageCode);

/**
 * @brief Function <b>report_symmetric_warning</b> sends a warning message to the reporter.
 *
 * Since the issue causing the warnign may occure thousands of times during a run, it is desirable
 * to throttle the number of times that a message is displayed.  If you desire to limit the number
 * of times a message is displayed, obtain a unique id via the get_next_message_id() function.  This
 * function allows you to assign the maximum display count.  When the count of messages reported
 * with this message id exceeds the display count, it is no longer displayed.  The default id zero
 * (0) is assigned which has an extremely large display count.
 *
 * @param message		a <b>char</b> const pointer to a message that is to be displayed.
 *
 * @param message_code		a <b>size_t</b> value of the message id for the message.
 *
 */
void report_symmetric_warning(const char *message, const MessageCode &message_code = MessageCode::s_defaultMessageCode);

/**
 * @brief Member function <b>report_deferred_warning</b> ...
 *
 * @param message		a <b>std::string</b> const ...
 *
 * @param aggregate		a <b>std::string</b> const ...
 *
 * @param message_code		a <b>MessageCode</b> const ...
 *
 */
void report_deferred_warning(const char *message, const char *aggregate, const MessageCode &message_code);

/**
 * @brief Class <b>RuntimeWarningAdHoc</b> reports an ad hoc warning message to the
 * report system.
 *
 *
 * For example:
 *
 * <PRE>
 *     if (adhoc_runtime_warning_condition)
 *       RuntimeWarningAdHoc() << "My useful message about " << some_data;
 *
 *     if (adhoc_runtime_warning_condition) {
 *       static MessageCode mc;
 *       RuntimeWarningAdHoc(mc) << "My useful message about " << some_data;
 *     }
 * </PRE>
 */
class RuntimeWarningAdHoc
{
public:
  /**
   * @brief Creates a new <b>RuntimeWarningAdHoc</b> instance, setting the message code.
   *
   * @param message_code	an <b>MessageCode</b> const reference to the message code associated
   *                            with this message.
   *
   */
  explicit RuntimeWarningAdHoc(MessageCode &message_code = MessageCode::s_defaultMessageCode);

  /**
   * @brief Destroys a <b>RuntimeWarningAdHoc</b> instance.
   *
   * The message is displayed by calling the report_warning() function.  However, if the count of
   * remaining messages for this message id is zero, the message is not displayed.
   *
   */
  ~RuntimeWarningAdHoc();

private:
  /**
   * @brief Make copy of <b>RuntimeWarningAdHoc</b> invalid.
   *
   */
  RuntimeWarningAdHoc(const RuntimeWarningAdHoc &);

  /**
   * @brief Make assignment of <b>RuntimeWarningAdHoc</b> invalid.
   *
   */
  RuntimeWarningAdHoc &operator=(const RuntimeWarningAdHoc &);
  
public:
  /**
   * @brief Member function <b>operator&lt;&lt;</b> passes the ios_base manipulator
   * to the output stream.
   *
   * @return			a <b>RuntimeWarningAdHoc</b> reference to this object
   */
  RuntimeWarningAdHoc &operator<<(OStreamFunctionPtr f) {
    f(message);
    return *this;
  }

  /**
   * @brief Member function <b>operator&lt;&lt;</b> passes the ios_base manipulator
   * to the output stream.
   *
   * @return			a <b>RuntimeWarningAdHoc</b> reference to this object
   */
  RuntimeWarningAdHoc &operator<<(IOSBaseFunctionPtr f) {
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
   * @return			a <b>RuntimeWarningAdHoc</b> reference to this object;
   */
  template <class T>
  RuntimeWarningAdHoc &operator<<(const T &t) {
    message << t;
    return *this;
  }

public:
  std::ostringstream    message;                ///< Stream to receive message content

private:
  const MessageCode     m_messageCode;          ///< Message id and uninitialized throttle
};


/**
 * @brief Class <b>RuntimeWarningSymmetric</b> reports a symmetric warning message to the report system.
 *
 * For example:
 *
 * <PRE>
 *     if (symmetric_runtime_warning_condition)
 *       RuntimeWarningSymmetric() << "My useful message about " << some_data;
 *
 *     if (symmetric_runtime_warning_condition) {
 *       static MessageCode mc;
 *       RuntimeWarningSymmetric(mc) << "My useful message about " << some_data;
 *     }
 * </PRE>
 */
class RuntimeWarningSymmetric
{
public:
  /**
   * @brief Creates a new <b>RuntimeWarning</b> instance, setting the message code.
   *
   * @param message_code	an <b>MessageCode</b> const reference to the message code associated
   *                            with this message.
   *
   */
  explicit RuntimeWarningSymmetric(MessageCode &message_code = MessageCode::s_defaultMessageCode);

  /**
   * @brief Destroys a <b>RuntimeWarningSymmetric</b> instance.
   *
   * The message is displayed by calling the report_symmetric_warning() function.  However, if the count of
   * remaining messages for this message id is zero, the message is not displayed.
   *
   */
  ~RuntimeWarningSymmetric();

private:
  /**
   * @brief Make copy of <b>RuntimeWarningSymmetric</b> invalid.
   *
   */
  RuntimeWarningSymmetric(const RuntimeWarningSymmetric &);

  /**
   * @brief Make assignment of <b>RuntimeWarningSymmetric</b> invalid.
   *
   */
  RuntimeWarningSymmetric &operator=(const RuntimeWarningSymmetric &);
  
public:
  /**
   * @brief Member function <b>operator&lt;&lt;</b> passes the ios_base manipulator
   * to the output stream.
   *
   * @return			a <b>RuntimeWarningSymmetric</b> reference to this object
   */
  RuntimeWarningSymmetric &operator<<(OStreamFunctionPtr f) {
    f(message);
    return *this;
  }

  /**
   * @brief Member function <b>operator&lt;&lt;</b> passes the ios_base manipulator
   * to the output stream.
   *
   * @return			a <b>RuntimeWarningSymmetric</b> reference to this object
   */
  RuntimeWarningSymmetric &operator<<(IOSBaseFunctionPtr f) {
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
   * @return			a <b>RuntimeWarningSymmetric</b> reference to this object;
   */
  template <class T>
  RuntimeWarningSymmetric &operator<<(const T &t) {
    message << t;
    return *this;
  }

public:
  std::ostringstream    message;                ///< Stream to receive message content

private:
  const MessageCode     m_messageCode;          ///< Message id and uninitialized throttle
};


/**
 * @brief Class <b>RuntimeWarningDeferred</b> reports a deferred warning message to the report
 * system.
 *
 * For example:
 *
 * <PRE>
 *     if (deferred_runtime_warning_condition) {
 *       static MessageCode mc;
 *       RuntimeWarningDeferred(mc) << "My useful message about " << some_data;
 *     }
 *
 *     if (deferred_runtime_warning_condition) {
 *       static MessageCode mc;
 *       RuntimeWarningDeferred x;
 *       x << "My useful message about " << some_data;
 *       x.aggregate << proc_specific_data;
 *     }
 * </PRE>
 */
class RuntimeWarningDeferred
{
public:
  /**
   * @brief Creates a new <b>RuntimeWarningDeferred</b> instance, setting the message code.
   *
   * @param message_code	an <b>MessageCode</b> const reference to the message code associated
   *                            with this message.
   *
   */
  explicit RuntimeWarningDeferred(const MessageCode &message_code);

  /**
   * @brief Destroys a <b>RuntimeWarning</b> instance.
   *
   * The message is displayed by calling the add_deferred_message() function.
   */
  ~RuntimeWarningDeferred();

private:
  /**
   * @brief Make copy of <b>RuntimeWarningDeferred</b> invalid.
   */
  RuntimeWarningDeferred(const RuntimeWarningDeferred &);

  /**
   * @brief Make assignment of <b>RuntimeWarningDeferred</b> invalid.
   */
  RuntimeWarningDeferred &operator=(const RuntimeWarningDeferred &);

public:
  /**
   * @brief Member function <b>operator&lt;&lt;</b> passes the ios_base manipulator
   * to the output stream.
   *
   * @return			a <b>RuntimeWarningDeferred</b> reference to this object
   */
  RuntimeWarningDeferred &operator<<(OStreamFunctionPtr f) {
    f(message);
    return *this;
  }

  /**
   * @brief Member function <b>operator&lt;&lt;</b> passes the ios_base manipulator
   * to the output stream.
   *
   * @return			a <b>RuntimeWarningDeferred</b> reference to this object
   */
  RuntimeWarningDeferred &operator<<(IOSBaseFunctionPtr f) {
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
   * @return			a <b>RuntimeWarningDeferred</b> reference to this object;
   */
  template <class T>
  RuntimeWarningDeferred &operator<<(const T &t) {
    message << t;
    return *this;
  }

public:
  std::ostringstream    message;                ///< Stream to receive message header content
  std::ostringstream    aggregate;              ///< Stream to receive message aggregate content

private:
  const MessageCode     m_messageCode;          ///< Message id and uninitialized throttle
};

///
/// @}
///

} // namespace stk

#endif // STK_UTIL_ENVIRONMENT_RUNTIMEWARNING_HPP
