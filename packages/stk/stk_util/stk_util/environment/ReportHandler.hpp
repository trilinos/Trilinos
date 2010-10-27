/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_UTIL_ENVIRONMENT_REPORTHANDLER_HPP
#define STK_UTIL_ENVIRONMENT_REPORTHANDLER_HPP

#include <iosfwd>
#include <string>
#include <sstream>

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
 * @brief ErrorHandler defines the signature of functions that can be used
 * to handle errors. expr is the expression of the failing error-check,
 * location is a raw code location (something like file:line, no prose),
 * and message is the error message.
 */
typedef void (*ErrorHandler)(const char* expr, const std::string& location, std::ostringstream& message);

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

/**
 * A function used to create and throw nice-looking exceptions (used by
 * the Throw* macros).
 */
void default_assert_handler(const char* expr,
                            const std::string& location,
                            std::ostringstream& message);

/**
 * Change the error handler for ThrowAssert and ThrowRequire.
 *
 * @return The previous error handler (useful if you want to restore it later)
 */
ErrorHandler set_assert_handler(ErrorHandler error_handler);

/**
 * Makes the call to the current assert handler
 */
void handle_assert(const char* expr,
                   const std::string& location,
                   std::ostringstream& message);

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

// The do-while is necessary to prevent usage of this macro from changing
// program semantics (e.g. dangling-else problem). The obvious implementation:
// if (expr) ; else throw ...
// is not adequate because it causes ambiguous else statements in this context:
// if (something)
//   ThrowRequire(foo);
// The compiler does not know whether the else statement that the macro inserts
// applies to the "if (something) " or the "if (expr)".
#define ThrowRequireMsg(expr, message)                                  \
  do {                                                                  \
    if ( !(expr) ) {                                                    \
      std::ostringstream stk_util_internal_throw_require_oss;           \
      stk_util_internal_throw_require_oss << message;                   \
      stk::handle_assert( #expr,                                        \
                          STR_TRACE,                                    \
                          stk_util_internal_throw_require_oss );        \
    }                                                                   \
  } while (false)

#define ThrowRequire(expr) ThrowRequireMsg(expr, "")

#ifdef NDEBUG
#  define ThrowAssert(expr)            ((void) (0))
#  define ThrowAssertMsg(expr,message) ((void) (0))
#else
#  define ThrowAssert(expr)            ThrowRequire(expr)
#  define ThrowAssertMsg(expr,message) ThrowRequireMsg(expr,message)
#endif

///
/// @}
///

#endif // STK_UTIL_ENVIRONMENT_REPORTHANDLER_HPP
