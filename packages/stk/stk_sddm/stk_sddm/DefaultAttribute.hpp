#ifndef SDDM_DEFAULT_ATTRIBUTE_HPP
#define SDDM_DEFAULT_ATTRIBUTE_HPP

#include <iosfwd>

#include <stk_sddm/AttributeInterface.hpp>

namespace stk {
namespace sddm {

class DefaultAttribute : public AttributeInterface<DefaultAttribute>
{
  virtual const char *name() const;
  
  virtual std::ostream &dump(std::ostream &os) const;
  
  virtual std::ostream &xml(std::ostream &os, size_t indent) const;
};


} // namespace sddm
} // namespace stk

#endif // SDDM_DEFAULT_ATTRIBUTE_HPP
