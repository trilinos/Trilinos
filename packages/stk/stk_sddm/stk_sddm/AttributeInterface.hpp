#ifndef SDDM_ATTRIBUTE_INTERFACE_HPP
#define SDDM_ATTRIBUTE_INTERFACE_HPP

#include <typeinfo>
#include <iosfwd>
#include <stddef.h>

namespace stk {
namespace sddm {

class AttributeInterfaceBase
{
protected:
  AttributeInterfaceBase()
  {}

public:
  virtual ~AttributeInterfaceBase()
  {}

  virtual const char *name() const = 0;
  virtual std::ostream &dump(std::ostream &os) const = 0;
  virtual std::ostream &xml(std::ostream &os, size_t indent) const = 0;
};


template <class T>
class AttributeInterface : public AttributeInterfaceBase {
protected:
  AttributeInterface()
  {}

public:
  typedef T Type;
  
  virtual ~AttributeInterface()
  {}

  static const std::type_info &type() {
    return typeid(Type);
  }
  
  virtual const char *name() const = 0;

  virtual std::ostream &dump(std::ostream &os) const = 0;
  
  virtual std::ostream &xml(std::ostream &os, size_t indent) const = 0;
};

} // namespace sddm
} // namespace stk

#endif // SDDM_ATTRIBUTE_INTERFACE_HPP
