#ifndef STK_UTIL_ENVIRONMENT_REPORTHANDLER_HPP
#define STK_UTIL_ENVIRONMENT_REPORTHANDLER_HPP

#include <iosfwd>
#include <string>

namespace stk {

///
/// @addtogroup runtime_message_detail
/// @{
///

/**
 * @brief Type definition REH is a pointer to a function of type void that takes a const
 * std::exception reference as a parameter.
 *
 */
typedef void (*REH)(const char *message, int type);

/**
 * @brief Function <b>default_report_handler</b> is the default
 * error reporter for sierra exceptions.  Note that it is implemented in Fmwk_sierra.C so
 * that it can participate
 *
 * @param message		a <b>char</b> const pointer to the message to be
 *				displayed.
 *
 * @param type			an <b>int</b> value of the type of message from the
 *				enumeration <b>type</b>
 *
 */
void default_report_handler(const char *message, int type);

/**
 * @brief Function <b>set_report_handler</b> sets the exception report function to be called when an
 * report_exception() is called.
 *
 * @param reh		a <b>REH</b> of the new exception reporter.
 *
 * @return		a <b>REH</b> to the previous exception reporter.
 */
REH set_report_handler(REH reh);

/**
 * @brief Function <b>report</b> calls the current exception reporter to report the message in
 * <b>x</b>.
 *
 * @param message	a <b>char</b> const pointer to the message to report.
 *
 * @param type		an <b>int</b> value of the type of message.
 *
 */
void report(const char *message, int type);

/**
 * @brief Function <b>source_relative_path</b> strips everything through "/src/",
 * "/include/" or "/App_" so that error message output doesn't mention names.
 *
 * @param s		a <b>std::string</b> const reference to the original path.
 *
 * @return		a <b>std::string</b> value of the stripped path.
 */
std::string source_relative_path(const std::string &path);

///
/// @}
///

} // namespace stk

///
/// @addtogroup runtime_message_detail
/// @{
///

/**
 * @brief Define statements to add __FILE___ and __LINE__.
 *
 * These just make it a little easier to add the __FILE__ and __LINE__ macros into the traceback
 * functions.  Inline functions do not work because the __FILE__ and __LINE__ expansions take place
 * before the compiler inlining.
 */
#define XSTR_TRACE_LINE(s) STR_TRACE_LINE(s)
#define STR_TRACE_LINE(s) #s

#ifdef __PRETTY_FUNCTION__

#define COUT_TRACE  " Function::Line="<<__PRETTY_FUNCTION__<<":"<<__LINE__
#define STR_TRACE   (std::string(__FILE__) +  ":" + XSTR_TRACE_LINE(__LINE__) + " in " + std::string(__PRETTY_FUNCTION__))

#else

#define COUT_TRACE  " File::Line="<<__FILE__<<":"<<__LINE__
#define STR_TRACE   (std::string(__FILE__) +  ":" + XSTR_TRACE_LINE(__LINE__))

#endif

#define StackTrace std::string(std::string("  exception thrown from ") + stk::source_relative_path(STR_TRACE))

#ifdef NDEBUG
#  define ThrowAssert(expr)		((void) (0))
#  define ThrowAssertMsg(expr,message)	((void) (0))
#else
#  define ThrowAssert(expr)		((expr) ? (void) 0 : throw std::runtime_error(std::string("Assertion ") + #expr + " failed\n" + StackTrace))
#  define ThrowAssertMsg(expr,message)	((expr) ? (void) 0 : throw std::runtime_error(std::string(message) + ", assertion " + #expr + " failed\n" + StackTrace))
#endif

#define ThrowRequire(expr)		((expr) ? (void) 0 : throw std::runtime_error(std::string("Requirement ") + #expr + " failed\n" + StackTrace))
#define ThrowRequireMsg(expr,message)	((expr) ? (void) 0 : throw std::runtime_error(std::string(message) + ", requirement " + #expr + " failed\n" + StackTrace))

///
/// @}
///

#endif // STK_UTIL_ENVIRONMENT_REPORTHANDLER_HPP
