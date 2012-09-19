#include <ostream>
#include <iomanip>

#include <stk_sddm/DefaultAttribute.hpp>

namespace stk {
namespace sddm {

const char *
DefaultAttribute::name() const {
  return "Default";
}
  
std::ostream &
DefaultAttribute::dump(std::ostream &os) const {
  os << "@" << name() << "()";
  return os;
}
  
std::ostream &
DefaultAttribute::xml(std::ostream &os, size_t indent) const {
  os << std::setw(indent*2) << "" << "<Attribute type=\"" << name() << "\"/>" << std::endl;

  return os;
}  

} // namespace sddm
} // namespace stk
