#ifndef SDDM_OCCURRENCE_ANNOTATION_HPP
#define SDDM_OCCURRENCE_ANNOTATION_HPP

#include <ostream>

#include <stk_sddm/AnnotationInterface.hpp>

namespace stk {
namespace sddm {

class OccurrenceAnnotation : public AnnotationInterface<OccurrenceAnnotation>
{
public:
  OccurrenceAnnotation(size_t min_count, size_t max_count)
    : m_minCount(min_count),
      m_maxCount(max_count)
  {}

  size_t getMinCount() const {
    return m_minCount;
  }
  
  size_t getMaxCount() const {
    return m_maxCount;
  }
  
  virtual const char *name() const;
  
  virtual std::ostream &dump(std::ostream &os) const;
  
  virtual std::ostream &xml(std::ostream &os) const;

  virtual bool validate(const Property &property) const;
  
protected:
  size_t        m_minCount;
  size_t        m_maxCount;
};

} // namespace sddm
} // namespace stk

#endif // SDDM_OCCURRENCE_ANNOTATION_HPP
