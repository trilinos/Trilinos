#include <stk_unit_test_utils/getOption.h>

namespace stk
{
namespace unit_test_util
{

namespace simple_fields {

bool has_option(const std::string& option) {
  return stk::unit_test_util::has_option(option);
}

std::string get_option(const std::string& option, const std::string defaultString) {
  return stk::unit_test_util::get_option(option, defaultString);
}

} // namespace simple_fields

}
}
