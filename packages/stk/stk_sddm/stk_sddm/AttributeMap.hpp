#ifndef SDDM_ATTRIBUTE_MAP_HPP
#define SDDM_ATTRIBUTE_MAP_HPP

#include <map>
#include <typeinfo>

#include <stk_sddm/AttributeInterface.hpp>

namespace stk {
namespace sddm {

class AttributeMap {
  typedef std::map<const std::type_info *, AttributeInterfaceBase *> Map;

public:
  typedef Map::const_iterator const_iterator;

  AttributeMap()
    : m_attributeMap()
  {}

  ~AttributeMap();

  const_iterator begin() const {
    return m_attributeMap.begin();
  }

  const_iterator end() const {
    return m_attributeMap.end();
  }

  bool empty() const {
    return m_attributeMap.empty();
  }

  template <class T>
  void addAttribute(AttributeInterface<T> *attribute) {
    addAttribute(&AttributeInterface<T>::type(), attribute);
  }

  template <class T>
  T *findAttribute() const {
    T *attribute = 0;

    Map::const_iterator it = m_attributeMap.find(&AttributeInterface<T>::type());

    if (it != m_attributeMap.end())
      attribute = static_cast<T *>((*it).second);
    return attribute;
  }

private:
  void addAttribute(const std::type_info *type, AttributeInterfaceBase *attribute);

private:
  Map                   m_attributeMap;
};

} // namespace sddm
} // namespace stk

#endif // SDDM_ATTRIBUTE_MAP_HPP
