#ifndef STK_SDDM_PRESTO_TRANSLATOR_HPP
#define STK_SDDM_PRESTO_TRANSLATOR_HPP

#include <iosfwd>

namespace stk {
namespace sddm {

class Property;

namespace presto {

void write(std::ostream &        os, const Property &      property);

} // namespace presto
} // namespace sddm
} // namespace stk

#endif // STK_SDDM_PRESTO_TRANSLATOR_HPP
