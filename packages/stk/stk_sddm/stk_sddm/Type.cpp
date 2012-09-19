#include <iomanip>
#include <sstream>
#include <iostream>

#include <stk_sddm/DataTypeTraits.hpp>
#include <stk_sddm/Type.hpp>

namespace stk {
namespace sddm {

#define INST(TYPE)                              \
template<>                                      \
AnyType *                                       \
Type<TYPE>::instance()                          \
{                                               \
  static Type<TYPE>     s_type;                 \
                                                \
  return &s_type;                               \
}

INST(int)
INST(double)
INST(std::string)
INST(std::vector<int>)
INST(std::vector<double>)
INST(std::vector<std::string>)

std::vector<const AnyType *> &
defaultTypes(
  std::vector<const AnyType *> &        default_type_vector)
{
  default_type_vector.push_back(Type<int>::instance());
  default_type_vector.push_back(Type<double>::instance());
  default_type_vector.push_back(Type<std::string>::instance());
  default_type_vector.push_back(Type<std::vector<int> >::instance());
  default_type_vector.push_back(Type<std::vector<double> >::instance());
  default_type_vector.push_back(Type<std::vector<std::string> >::instance());

  return default_type_vector;
}

} // namespace sddm
} // namespace stk
