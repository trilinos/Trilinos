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

#ifndef STK_UTIL_UTIL_REPORTHANDLER_HPP
#define STK_UTIL_UTIL_REPORTHANDLER_HPP

#include <sstream>    // for ostringstream
#include <stdexcept>  // for logic_error, runtime_error
#include <string>     // for operator+, allocator, string, char_traits
#include <type_traits>

#include "stk_util/diag/String.hpp"
#include "stk_util/stk_kokkos_macros.h"  // for STK_FUNCTION

#ifdef STK_ENABLE_GPU_BUT_NO_RDC
#include <Kokkos_Core.hpp>
#endif

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
 * @brief Function <b>source_relative_path</b> strips everything through
 * "/src/", "/include/", "/App_", or "/stk_" so that error message output
 * doesn't mention names.
 *
 * @param path  a <b>std::string</b> const reference to the original path.
 *
 * @return      a <b>std::string</b> value of the stripped path.
 */
std::string source_relative_path(const std::string &path);

/**
 * A function used to create and throw nice-looking exceptions for assertion
 * failures.
 */
void default_assert_handler(const char* expr,
                            const std::string& location,
                            std::ostringstream& message);

/**
 * A function used to create and throw nice-looking exceptions for runtime
 * errors.
 */
void default_error_handler(const char* expr,
                           const std::string& location,
                           std::ostringstream& message);

/**
 * This is the same as default_error_handler but does not output traceback
 * */
void clean_error_handler(const char* expr,
                           const std::string& location,
                           std::ostringstream& message);

/**
 * A function used to create and throw nice-looking exceptions for invalid
 * argument errors.
 */
void default_invalid_arg_handler(const char* expr,
                                 const std::string& location,
                                 std::ostringstream& message);

/**
 * Change the error handler for ThrowAssert and ThrowRequire.
 *
 * @return The previous error handler (useful if you want to restore it later)
 */
ErrorHandler set_assert_handler(ErrorHandler handler);

/**
 * Change the error handler for ThrowError.
 *
 * @return The previous error handler (useful if you want to restore it later)
 */
ErrorHandler set_error_handler(ErrorHandler handler);

/**
 * Change the error handler for ThrowInvalidArg.
 *
 * @return The previous error handler (useful if you want to restore it later)
 */
ErrorHandler set_invalid_arg_handler(ErrorHandler handler);

/**
 * Makes the call to the current assert handler
 */
void handle_assert(const char* expr,
                   const std::string& location,
                   std::ostringstream& message);

/**
 * Makes the call to the current error handler
 */
void handle_error(const char* expr,
                   const std::string& location,
                   std::ostringstream& message);

/**
 * Makes the call to the current invalid_arg handler
 */
void handle_invalid_arg(const char* expr,
                        const std::string& location,
                        std::ostringstream& message);


///
/// @}
///

std::ostream & output_stacktrace(std::ostream & os);

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
#define XSTK_STR_TRACE_LINE(s) STK_STR_TRACE_LINE(s)
#define STK_STR_TRACE_LINE(s) #s

#ifdef __PRETTY_FUNCTION__

#define COUT_TRACE  " Function::Line=" __PRETTY_FUNCTION__ ":" __LINE__
#define STK_STR_TRACE   (__FILE__ ":" XSTK_STR_TRACE_LINE(__LINE__) " in " __PRETTY_FUNCTION__)

#else

#define COUT_TRACE  " File::Line=" __FILE__ ":" __LINE__
#define STK_STR_TRACE   (__FILE__ ":" XSTK_STR_TRACE_LINE(__LINE__))

#endif

#define StackTrace std::string(std::string("  exception thrown from ") + stk::source_relative_path(STK_STR_TRACE))

#define STRINGIZE(x) STRINGIZE2(x)
#define STRINGIZE2(x) #x
#define LINE_STRING STRINGIZE(__LINE__)

// The do-while is necessary to prevent usage of this macro from changing
// program semantics (e.g. dangling-else problem). The obvious implementation:
// if (expr) ; else throw ...
// is not adequate because it causes ambiguous else statements in this context:
// if (something)
//   STK_ThrowRequire(foo);
// The compiler does not know whether the else statement that the macro inserts
// applies to the "if (something) " or the "if (expr)".

namespace stk::impl
{
template <typename T>
class is_string_literal
{
  using literal_t = const char (&)[sizeof(T)];

 public:
  static constexpr bool value = std::is_same_v<T, literal_t>;
};

template <typename T, typename... Args>
struct is_same_as_any {
  static constexpr bool value = false;
};
template <typename T, typename First, typename... Args>
struct is_same_as_any<T, First, Args...> {
  static constexpr bool value = std::is_same_v<T, First> || is_same_as_any<T, Args...>::value;
};

template <typename T>
class is_valid_throw_condition
{
  using raw_t = std::decay_t<T>;

 public:
  static constexpr bool value =
      !is_same_as_any<raw_t, sierra::String, std::string, const char*, char*>::value && !is_string_literal<T>::value;
};
template <typename T>
inline auto eval_test_condition(const T& val)
{
  static_assert(stk::impl::is_valid_throw_condition<T>::value,
      "Cannot use string type as the condition in STK_ThrowRequire*. "
      "Use (ptr != nullptr) instead.");
  return !val;
}
}  // namespace stk::impl

#define STK_ThrowGenericCond(expr, message, handler)                                                           \
  do {                                                                                                         \
    if (stk::impl::eval_test_condition(expr)) {                                                                \
      std::ostringstream stk_util_internal_throw_require_oss;                                                  \
      stk_util_internal_throw_require_oss << message;                                                          \
      std::ostringstream stk_util_internal_throw_require_loc_oss;                                              \
      stk_util_internal_throw_require_loc_oss << stk::source_relative_path(STK_STR_TRACE) << "\n";             \
      stk::output_stacktrace(stk_util_internal_throw_require_loc_oss);                                         \
      stk::handler(#expr, stk_util_internal_throw_require_loc_oss.str(), stk_util_internal_throw_require_oss); \
    }                                                                                                          \
  } while (false)

inline void ThrowMsgHost(bool /* expr */, const char* exprString, const char* message, const std::string& location)
{
  std::ostringstream stk_util_internal_throw_require_loc_oss;
  stk_util_internal_throw_require_loc_oss << stk::source_relative_path(location) << "\n";
  stk::output_stacktrace(stk_util_internal_throw_require_loc_oss);
  throw std::logic_error(
    std::string("Requirement( ") + exprString + " ) FAILED\n" +
    "Error occurred at: " + stk_util_internal_throw_require_loc_oss.str() + "\n" +
    "Error: " + message + "\n");
}

#ifdef STK_ENABLE_GPU_BUT_NO_RDC
STK_INLINE_FUNCTION void ThrowMsgDevice(const char * message)
{ 
  Kokkos::abort(message);
}
#else
STK_FUNCTION void ThrowMsgDevice(const char * message);
#endif

inline void ThrowHost(bool /* expr */, const char* exprString, const std::string& location)
{
  std::ostringstream stk_util_internal_throw_require_loc_oss;
  stk_util_internal_throw_require_loc_oss << stk::source_relative_path(location) << "\n";
  stk::output_stacktrace(stk_util_internal_throw_require_loc_oss);
  throw std::logic_error(
    std::string("Requirement( ") + exprString + " ) FAILED\n" +
    "Error occurred at: " + stk_util_internal_throw_require_loc_oss.str() + "\n");
}

inline void ThrowErrorMsgHost(const char * message, const std::string & location)
{
  std::ostringstream stk_util_internal_throw_require_loc_oss;
  stk_util_internal_throw_require_loc_oss << stk::source_relative_path(location) << "\n";
  stk::output_stacktrace(stk_util_internal_throw_require_loc_oss);
  throw std::runtime_error(
    std::string("Error occurred at: ") + stk_util_internal_throw_require_loc_oss.str() + "\n" +
    "Error: " + message + "\n");
}

#ifdef STK_ENABLE_GPU_BUT_NO_RDC
STK_INLINE_FUNCTION void ThrowErrorMsgDevice(const char * message)
{ 
  Kokkos::abort(message);
}
#else
STK_FUNCTION void ThrowErrorMsgDevice(const char * message);
#endif

// This generic macro is for unconditional throws. We pass "" as the expr
// string, the handler should be smart enough to realize that this means there
// was not expression checked, AKA, this throw was unconditional.
#define STK_ThrowGeneric(message, handler)                                 \
  do {                                                                     \
    std::ostringstream stk_util_internal_throw_require_oss;                \
    stk_util_internal_throw_require_oss << message;                        \
    std::ostringstream stk_util_internal_throw_require_loc_oss;            \
    stk_util_internal_throw_require_loc_oss <<                             \
      stk::source_relative_path(STK_STR_TRACE) << "\n";                    \
    stk::output_stacktrace(stk_util_internal_throw_require_loc_oss);       \
    stk::handler( "",                                                      \
                  stk_util_internal_throw_require_loc_oss.str(),           \
                  stk_util_internal_throw_require_oss );                   \
} while (false)

// The macros below define the exceptions that we want to support within
// STK. The intent is that STK developers will never call throw XXX
// directly. This will give us full control over the exceptions being generated
// by STK and how they are handled. These macros are also designed to make
// it as easy as possible for developers to throw exceptions that have good
// error messages and to reduce the volume of coded needed for error handling.
//
// We currently support the following exceptions in STK:
//   logic_error <-> STK_ThrowAssert, STK_ThrowAsserMsg, STK_ThrowRequire,
//                   STK_ThrowRequireMsg
//   runtime_error <-> STK_ThrowErrorMsgIf, STK_ThrowErrorMsg
//   invalid_argument <-> STK_ThrowInvalidArgMsgIf, STK_ThrowInvalidArgIf
//
// Please note the logic of the errors is the opposite of the asserts. The
// asserts will throw exceptions if the given expression is false; for the
// error macros, exceptions will be thrown if the given expression is true.
//
// USE:
//     All of the following have versions that do not require a message, but
//     we strongly encourage developers to use the versions that take the
//     message.
//
//   ASSERTS:
//     STK_ThrowAssertMsg(expr, message);
//       If NDEBUG is not defined, throw a logic error if expr evaluates to
//       false, adding message to the error message . Use this for expensive
//       logic-mistake checks that could impact performance.
//
//     STK_ThrowRequireMsg(code, message);
//       Always throw a logic error if expr evaluates to false, adding message
//       to error message. Use this for inexpensive logic-mistake checks
//       that do not impact performance.
//
//   ERRORS:
//     STK_ThrowErrorMsgIf(expr, message);
//       Throw a runtime error if expr evaluates to true, adding message to
//       the error message. Use this to generate errors dealing with system
//       errors or other errors that do not involve invalid parameters being
//       passed to functions.
//
//     STK_ThrowInvalidArgMsgIf(expr, message);
//       Throw an invalid_argument error if expr evaluates to true, adding
//       message to the error message. Use this to generate errors dealing with
//       users passing invalid arguments to functions in the API.
//
// EXAMPLES:
//
// 1) Require that i equals j, demonstate use of put-tos in the message arg
//   STK_ThrowRequireMsg(i == j, "i(" << i << ") != j(" << j << ")");
//
// 2) Check method argument foo is not NULL
//   STK_ThrowInvalidArgMsgIf(foo != NULL, "Arg foo is NULL");

#define STK_ThrowRequireWithSierraHelpMsg(expr) STK_ThrowGenericCond(expr, "Program error. Contact sierra-help@sandia.gov for support.", handle_assert)
#define STK_ThrowRequireMsg(expr, message) STK_ThrowGenericCond(expr, message, handle_assert)
#define STK_ThrowRequire(expr) STK_ThrowRequireMsg(expr, "")

#ifndef __HIP_DEVICE_COMPILE__

#ifdef NDEBUG
#  define STK_ThrowAssert(expr)                                              (static_cast<void>(0))
#  define STK_ThrowAssertMsg(expr,message)                                   (static_cast<void>(0))
#else 
#  define STK_ThrowAssert(expr)                                              STK_ThrowRequire(expr)
#  define STK_ThrowAssertMsg(expr,message)                                   STK_ThrowRequireMsg(expr,message)
#endif

#define STK_ThrowErrorMsgIf(expr, message)                              STK_ThrowGenericCond( !(expr), message, handle_error)
#define STK_ThrowErrorIf(expr)                                          STK_ThrowErrorMsgIf(expr, "")
#define STK_ThrowErrorMsg(message)                                      STK_ThrowGeneric( message, handle_error )

#define STK_ThrowInvalidArgMsgIf(expr, message)                                    STK_ThrowGenericCond( !(expr), message, handle_invalid_arg)
#define STK_ThrowInvalidArgIf(expr)                                                STK_ThrowInvalidArgMsgIf(expr, "")

#else
// FIXME: unsupported indirect call to function on HIP-Clang
#define STK_ThrowAssert(expr)
#define STK_ThrowAssertMsg(expr,message)

#define STK_ThrowErrorMsgIf(expr, message)
#define STK_ThrowErrorIf(expr)
#define STK_ThrowErrorMsg(message)

#define STK_ThrowInvalidArgMsgIf(expr, message)
#define STK_ThrowInvalidArgIf(expr)
#endif

#if ((defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0)) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__))
#define STK_NGP_ThrowRequireMsg(expr, message)               \
  do {                                                       \
    const bool __stk_expr_res = bool(expr);                  \
    if (!__stk_expr_res) {                                   \
      ThrowMsgDevice(message ": " __FILE__ ":" LINE_STRING); \
    }                                                        \
  } while (false);
#else
#define STK_NGP_ThrowRequireMsg(expr, message)                     \
  do {                                                             \
    const bool __stk_expr_res = bool(expr);                        \
    if (!__stk_expr_res) {                                         \
      ThrowMsgHost(__stk_expr_res, #expr, message, STK_STR_TRACE); \
    }                                                              \
  } while (false);
#endif

#if ((defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0)) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__))
#define STK_NGP_ThrowRequire(expr)                              \
  do {                                                          \
    const bool __stk_expr_res = bool(expr);                     \
    if (!__stk_expr_res) {                                      \
      ThrowMsgDevice("(" #expr "): " __FILE__ ":" LINE_STRING); \
    }                                                           \
  } while (false);
#else
#define STK_NGP_ThrowRequire(expr)                     \
  do {                                                 \
    const bool __stk_expr_res = bool(expr);            \
    if (!__stk_expr_res) {                             \
      ThrowHost(__stk_expr_res, #expr, STK_STR_TRACE); \
    }                                                  \
  } while (false);
#endif

#ifdef NDEBUG
#  define STK_NGP_ThrowAssert(expr)                                                  (static_cast<void>(0))
#  define STK_NGP_ThrowAssertMsg(expr,message)                                       (static_cast<void>(0))
#else
#  define STK_NGP_ThrowAssert(expr)                                                  STK_NGP_ThrowRequire(expr)
#  define STK_NGP_ThrowAssertMsg(expr,message)                                       STK_NGP_ThrowRequireMsg(expr, message)
#endif

#if ((defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0)) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__))
#define STK_NGP_ThrowErrorMsgIf(expr, message)                                        STK_NGP_ThrowRequireMsg(!(expr), message);
#else
#define STK_NGP_ThrowErrorMsgIf(expr, message)                              \
  do {                                                                      \
    const bool __stk_expr_res = bool(expr);                                 \
    if (__stk_expr_res) {                                                   \
      ThrowMsgHost(__stk_expr_res, "!(" #expr ")", message, STK_STR_TRACE); \
    }                                                                       \
  } while (false);
#endif

#if ((defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0)) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__))
#define STK_NGP_ThrowErrorIf(expr)                                     STK_NGP_ThrowRequireMsg(!(expr), "!(" #expr ")");
#else
#define STK_NGP_ThrowErrorIf(expr)                              \
  do {                                                          \
    const bool __stk_expr_res = bool(expr);                     \
    if (__stk_expr_res) {                                       \
      ThrowHost(__stk_expr_res, "!(" #expr ")", STK_STR_TRACE); \
    }                                                           \
  } while (false);
#endif

#if ((defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0)) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__))
#define STK_NGP_ThrowErrorMsg(message)                                      ThrowErrorMsgDevice(message ": " __FILE__ ":" LINE_STRING);
#else
#define STK_NGP_ThrowErrorMsg(message)                                      ThrowErrorMsgHost(message, STK_STR_TRACE);
#endif


/**
 * @ingroup Exception
 * @brief Macro WarnTrace makes a pretty warning message with file and line number.
 */
#define WarnTrace std::string(std::string("  warning at ") + stk::source_relative_path(STK_STR_TRACE))

/**
 * @ingroup Exception
 * @brief Macro ErrorTrace makes a pretty error message with file and line number.
 */
#define ErrorTrace std::string(std::string("  error thrown from ") + stk::source_relative_path(STK_STR_TRACE))

///
/// @}
///

#endif // STK_UTIL_UTIL_REPORTHANDLER_HPP
