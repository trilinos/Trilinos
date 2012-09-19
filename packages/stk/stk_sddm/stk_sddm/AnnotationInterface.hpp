#ifndef SDDM_ANNOTATION_INTERFACE_HPP
#define SDDM_ANNOTATION_INTERFACE_HPP

#include <typeinfo>
#include <iosfwd>

namespace stk {
namespace sddm {

class Property;

class AnnotationInterfaceBase
{
protected:
  AnnotationInterfaceBase()
  {}
public:
  virtual ~AnnotationInterfaceBase()
  {}

  virtual const char *name() const = 0;
  virtual std::ostream &dump(std::ostream &os) const = 0;
  virtual std::ostream &xml(std::ostream &os) const = 0;
  virtual bool validate(const Property &property) const = 0;
};


template <class T>
class AnnotationInterface : public AnnotationInterfaceBase
{
protected:
  AnnotationInterface()
  {}

public:
  virtual ~AnnotationInterface()
  {}

  static const std::type_info &type() {
    return typeid(T);
  }
  
  virtual const char *name() const = 0;

  virtual std::ostream &dump(std::ostream &os) const = 0;
  
  virtual std::ostream &xml(std::ostream &os) const = 0;

  virtual bool validate(const Property &property) const = 0;
};

} // namespace sddm
} // namespace stk

#endif // SDDM_ANNOTATION_INTERFACE_HPP
