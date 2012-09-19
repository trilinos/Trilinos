#ifndef STK_UTIL_PARALLEL_ExceptionReport_hpp
#define STK_UTIL_PARALLEL_ExceptionReport_hpp

#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/environment/RuntimeWarning.hpp>
#include <stk_util/environment/RuntimeDoomed.hpp>

#include <stk_util/util/Fortran.hpp>
#include <stk_util/parallel/Exception.hpp>


/**
 * @ingroup Exception
 * @brief Macro WarnTrace makes a pretty warning message with file and line number.
 *
 */
#define WarnTrace std::string(std::string("  warning at ") + stk::source_relative_path(STR_TRACE))

/**
 * @ingroup Exception
 * @brief Macro ErrorTrace makes a pretty error message with file and line number.
 *
 */
#define ErrorTrace std::string(std::string("  error thrown from ") + stk::source_relative_path(STR_TRACE))

namespace sierra {

  enum ErrorDieEnum{DIE_ON_WARN=0, DIE_ON_ERROR=1, DIE_ON_MESSAGE=2};

typedef stk::MessageCode MessageCode;

typedef stk::RuntimeWarningAdHoc RuntimeWarning;                ///< Deprecated
typedef stk::RuntimeWarningSymmetric RuntimeWarningP0;          ///< Deprecated

typedef stk::RuntimeWarningAdHoc RuntimeWarningAdHoc;
typedef stk::RuntimeWarningSymmetric RuntimeWarningSymmetric;
typedef stk::RuntimeWarningDeferred RuntimeWarningDeferred;

typedef stk::RuntimeDoomedAdHoc RuntimeDoomed;                  ///< Deprecated
typedef stk::RuntimeDoomedSymmetric RuntimeDoomedP0;            ///< Deprecated

typedef stk::RuntimeDoomedAdHoc RuntimeDoomedAdHoc;
typedef stk::RuntimeDoomedSymmetric RuntimeDoomedSymmetric;
typedef stk::RuntimeDoomedDeferred RuntimeDoomedDeferred;

void set_test_error_messages_file(const std::string &test_error_messages_path);

std::ofstream *get_test_error_messages_file();

 void set_test_error_messages_die_on_first_message(std::vector<ErrorDieEnum> errorTypes);

bool get_test_error_messages_die_on_first_warning();
bool get_test_error_messages_die_on_first_error();


///
/// @}
///

} // namespace sierra

/**
 * @ingroup Exception
 * @brief Function <b>report_error</b> reports a message from a Fortran function to
 * the exception reporting system.  If the <b>int_val</b> is 1, the message is printed
 * to sierra::Env::outputP0().  If it is 2, a sierra::RuntimeWarning is issued.  And, if
 * it is 3, a sierra::RuntimeError is thrown.
 *
 * @param int_val		a <b>int</b> value of the type of message to
 *				report.
 *
 * @param message		a <b>char</b> const pointer to the start of the
 *				message.
 *
 * @param message_length	a <b>int</b> const value the length of the message.
 *
 */
extern "C"
  void SIERRA_FORTRAN(report_error)(int &int_val, const char *message, const int message_length);

#endif // STK_UTIL_PARALLEL_ExceptionReport_hpp
