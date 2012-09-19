#ifndef SDDM_ANNOTATION_MAP_HPP
#define SDDM_ANNOTATION_MAP_HPP

#include <vector>
#include <stdexcept>
#include <sstream>
#include <typeinfo>

#include <stk_sddm/AnnotationInterface.hpp>

namespace stk {
namespace sddm {

class AnnotationMap {
  struct TypeIndex {
    TypeIndex(const std::type_info *type, size_t begin, size_t end)
      : m_type(type),
        m_begin(begin),
        m_end(end)
    {}

    const std::type_info *      m_type;
    size_t                      m_begin;
    size_t                      m_end;
  };

  typedef std::vector<TypeIndex> TypeIndexVector;
  typedef std::vector<AnnotationInterfaceBase *>  AnnotationVector;

public:
  typedef AnnotationVector::const_iterator const_iterator;

  AnnotationMap()
    : m_typeVector(),
      m_annotationVector()
  {}
  
  ~AnnotationMap();
  
  const_iterator begin() const {
    return m_annotationVector.begin();
  }
  
  const_iterator end() const {
    return m_annotationVector.end();
  }
  
  bool empty() const {
    return m_annotationVector.empty();
  }
  
  template <class T>
  void addAnnotation(AnnotationInterface<T> *annotation) {
    addAnnotation(&annotation->type(), annotation);
  }

  template <class T, class Out>
  void findAnnotation(Out out) const {
    TypeIndexVector::const_iterator type_it = findAnnotation(&typeid(AnnotationInterface<T>));

    if (type_it != m_typeVector.end())
      for (size_t i = (*type_it).m_begin; i != (*type_it).m_end; ++i, ++out)
        *out = static_cast<const T *>(m_annotationVector[i]);
  }
  
private:
  void addAnnotation(const std::type_info *type, AnnotationInterfaceBase *annotation);
  TypeIndexVector::const_iterator findAnnotation(const std::type_info *type) const;
  TypeIndexVector::iterator findAnnotation(const std::type_info *type);

private:
  TypeIndexVector       m_typeVector;
  AnnotationVector      m_annotationVector;
};

} // namespace sddm
} // namespace stk

#endif // SDDM_ANNOTATION_MAP_HPP
