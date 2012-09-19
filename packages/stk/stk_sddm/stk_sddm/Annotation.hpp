#ifndef SDDM_ANNOTATIONS_HPP
#define SDDM_ANNOTATIONS_HPP

#include <map>
#include <stdexcept>
#include <string>
#include <typeinfo>

namespace stk {
namespace sddm {

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
  virtual bool validate() const = 0;
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

  virtual bool validate() const = 0;
};


class AnnotationMap {
  typedef std::map<const std::type_info *, AnnotationInterfaceBase *> Map;
public:
  template <class T>
  void addAnnotation(AnnotationInterface<T> *annotation) {
    m_annotationMap[&annotation->type()] = annotation;
  }
    
  template <class T>
  T *findAnnotation() const {
    T *annotation = 0;
  
    Map::const_iterator it = m_annotationMap.find(&AnnotationInterface<T>::type());

    if (it != m_annotationMap.end())
      annotation = static_cast<T *>((*it).second);
    return annotation;
  }

  template <class T>
  bool exists() const {
    return findAnnotation<T>(m_annotationMap);
  }
  
  template <class T>
  T &getAnnotation() const {
    T *annotation = findAnnotation<T>(m_annotationMap);

    if (!annotation) {
      std::ostringstream oss;

      oss << "Annotation " << T::name() << " was not added to annotation map";
    
      throw std::runtime_error(oss.str());
    }
    return *annotation;
  }
  
private:
  Map                   m_annotationMap;
};


class SomeOfAnnotation : public AnnotationInterface<SomeOfAnnotation>
{
protected:
  SomeOfAnnotation(size_t min_count, size_t max_count)
    : m_min(min),
      m_max(max)
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
  
  virtual std::ostream &describe(std::ostream &oss) const;

  virtual bool validate() const {
    return true;
  }
  
protected:
  size_t        m_min;
  size_t        m_max;
};

} // namespace sddm
} // namespace stk

#endif // SDDM_ANNOTATIONS_HPP
