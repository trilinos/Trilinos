#include <stk_sddm/AnnotationMap.hpp>

namespace stk {
namespace sddm {

AnnotationMap::~AnnotationMap()
{
  for (AnnotationVector::iterator it = m_annotationVector.begin(); it != m_annotationVector.end(); ++it)
    delete *it;
}

void
AnnotationMap::addAnnotation(
  const std::type_info *        type,
  AnnotationInterfaceBase *     annotation)
{
  TypeIndexVector::iterator it = findAnnotation(type);

  if (it == m_typeVector.end()) {
    m_typeVector.push_back(TypeIndex(type, m_annotationVector.size(), m_annotationVector.size() + 1));
    m_annotationVector.push_back(annotation);
  }
  else {
    AnnotationVector::iterator annotation_it = m_annotationVector.begin() + (*it).m_end;
    m_annotationVector.insert(annotation_it, annotation);
    ++(*it).m_end;

    ++it;
    
    while (it != m_typeVector.end()) {
      ++(*it).m_begin;
      ++(*it).m_end;
      ++it;
    }
  }
}


AnnotationMap::TypeIndexVector::const_iterator
AnnotationMap::findAnnotation(
  const std::type_info *        type) const
{
  TypeIndexVector::const_iterator it = m_typeVector.begin();
  
  while (it != m_typeVector.end() && (*(*it).m_type).before(*type))
    ++it;

  return it;
}

AnnotationMap::TypeIndexVector::iterator
AnnotationMap::findAnnotation(
  const std::type_info *        type)
{
  TypeIndexVector::iterator it = m_typeVector.begin();
  
  while (it != m_typeVector.end() && (*(*it).m_type).before(*type))
    ++it;

  return it;
}

} // namespace sddm
} // namespace stk

