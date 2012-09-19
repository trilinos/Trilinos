#include <stdexcept>
#include <sstream>
#include <iostream>

namespace stk {
namespace sddm {

std::runtime_error
cast_error(
  const char *  	from_type_name,
  const char *          to_type_name)
{
  std::ostringstream oss;
  oss << "Cannot cast property of type " << from_type_name
      << " to type " << to_type_name;

  return std::runtime_error(oss.str());
}

} // namespace sddm
} // namespace stk
