#ifndef SDDM_SOMEOF_ANNOTATION_HPP
#define SDDM_SOMEOF_ANNOTATION_HPP

#include <ostream>

#include <stk_sddm/AnnotationInterface.hpp>

namespace stk {
namespace sddm {

class SomeOfAnnotation : public AnnotationInterface<SomeOfAnnotation>
{
protected:
  SomeOfAnnotation(size_t min_count, size_t max_count)
    : m_minCount(min_count),
      m_maxCount(max_count)
  {}

  virtual const char *name() const {
    return "SomeOf";
  }
  
  virtual std::ostream &dump(std::ostream &os) const {
    return os;
  }
  
  virtual std::ostream &xml(std::ostream &os) const {
    return os;
  }
  
//  virtual std::ostream &describe(std::ostream &oss) const;

  virtual bool validate() const {
    return true;
  }
  
protected:
  size_t        m_minCount;
  size_t        m_maxCount;
};

class OneOfAnnotation : public SomeOfAnnotation
{
protected:
  OneOfAnnotation()
    : SomeOfAnnotation(1, 1)
  {}

  virtual const char *name() const {
    return "OneOf";
  }
  
  virtual std::ostream &dump(std::ostream &os) const {
    return os;
  }
  
  virtual std::ostream &xml(std::ostream &os) const {
    return os;
  }
  
//  virtual std::ostream &describe(std::ostream &oss) const;

  virtual bool validate() const {
    return true;
  }
};

} // namespace sddm
} // namespace stk

#endif // SDDM_SOMEOF_ANNOTATION_HPP
