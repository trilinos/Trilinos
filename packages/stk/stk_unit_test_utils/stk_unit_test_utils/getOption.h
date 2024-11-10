#ifndef UNITTESTUTILS_OPTIONS_PARSING
#define UNITTESTUTILS_OPTIONS_PARSING

#include <sstream>
#include "stk_util/stk_config.h"
#include <stk_unit_test_utils/CommandLineArgs.hpp>
#include <string>

namespace stk
{
namespace unit_test_util
{

inline bool has_option(const std::string& option)
{
  stk::unit_test_util::GlobalCommandLineArguments& args = stk::unit_test_util::GlobalCommandLineArguments::self();
  if (args.get_argv() != nullptr) {
    for (int i = 0; i < args.get_argc(); i++) {
      std::string input_argv(args.get_argv()[i]);
      if (option == input_argv) {
        return true;
      }
    }
  }
    return false;
}

inline std::string get_option(const std::string& option, const std::string defaultString="no")
{
  stk::unit_test_util::GlobalCommandLineArguments& args = stk::unit_test_util::GlobalCommandLineArguments::self();
  std::string returnValue = defaultString;
  if (args.get_argv() != nullptr) {
    for (int i = 0; i < args.get_argc(); i++) {
      std::string input_argv(args.get_argv()[i]);
      if (option == input_argv) {
        if ((i + 1) < args.get_argc()) {
          returnValue = std::string(args.get_argv()[i + 1]);
        }
        break;
      }
    }
  }
    return returnValue;
}

template <typename T>
T get_command_line_option(const std::string &option, const T &defaultValue)
{
    std::ostringstream os;
    os << defaultValue;
    std::string str = get_option(option, os.str());
    std::istringstream ss(str);
    T val(defaultValue);
    ss >> val;
    return val;
}

namespace simple_fields {

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
bool has_option(const std::string& option);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
std::string get_option(const std::string& option, const std::string defaultString="no");

template <typename T>
STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
T get_command_line_option(const std::string &option, const T &defaultValue) {
  return stk::unit_test_util::get_command_line_option(option, defaultValue);
}

} // namespace simple_fields

}
}

#endif
