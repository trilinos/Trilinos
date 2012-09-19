#ifndef STK_SDDM_PROPERTY_IO_HPP
#define STK_SDDM_PROPERTY_IO_HPP

#include <stk_sddm/Property.hpp>

namespace stk {
namespace sddm {

std::ostream &print_property(std::ostream &os, bool feature_coverage, const Property &property);

std::ostream &print_xml(std::ostream &os, const Property &property, size_t indent = 0);
  

inline
std::ostream &operator<<(std::ostream &os, const AnyValue &any_value) {
  return any_value.dump(os);
}

std::ostream &operator<<(std::ostream &os, const Property &property);

std::ostream &operator>>(std::ostream &is, AnyValue &any_value);

template <class T>
inline
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
  for (typename std::vector<T>::const_iterator it = v.begin(); it != v.end(); ++it) {
    if (it != v.begin())
      os << " ";    
    os << (*it);
  }
  
  return os;
}

} // namespace sddm
} // namespace stk

#endif // STK_SDDM_PROPERTY_IO_HPP

