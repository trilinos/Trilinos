#include <limits>

#include <stk_sddm/OccurrenceAnnotation.hpp>
#include <stk_sddm/Property.hpp>
#include <stk_sddm/Algorithm.hpp>

namespace stk {
namespace sddm {

const char *
OccurrenceAnnotation::name() const
{
  return "SomeOf";
}


std::ostream &
OccurrenceAnnotation::dump(
  std::ostream &        os) const
{
  return os;
}


std::ostream &
OccurrenceAnnotation::xml(
  std::ostream &        os) const
{
  os << "    <Annotation type=\"OccurrenceCount\">" << std::endl
     << "      <Min>" << m_minCount << "</Min>" << std::endl;
  if (m_maxCount == std::numeric_limits<size_t>::max())
    os << "      <Max>unlimited</Max>" << std::endl;
  else
    os << "      <Max>" << m_maxCount << "</Max>" << std::endl;
  os << "    </Annotation>" << std::endl;

  return os;
}


bool
OccurrenceAnnotation::validate(
  const Property &      property) const
{
  const Property *parent = property.getParent();
  if (parent) {
    size_t count = 0;;
    stk::sddm::query(parent->begin(), parent->end(), counter(count), property.getName());

    if (count < m_minCount || count > m_maxCount) {
      std::ostringstream oss;
      
      if (count < m_minCount)
        oss << " must occur at least " << m_minCount << " times";
      else if (count > m_maxCount)
        oss << " must not occur more than " << m_maxCount << " times";
        
      property.reportValidationError(oss.str());
      return false;
    }
  }
  
  return true;
}

} // namespace sddm
} // namespace stk
