#ifndef STK_SDDM_PROPERTY_HPP
#define STK_SDDM_PROPERTY_HPP

#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>
#include <map>

#include <stk_sddm/ci_traits.hpp>
#include <stk_sddm/DataTypeTraits.hpp>
#include <stk_sddm/Type.hpp>
#include <stk_sddm/AttributeMap.hpp>
#include <stk_sddm/AttributeInterface.hpp>
#include <stk_sddm/DefaultAttribute.hpp>
#include <stk_sddm/SourceAttribute.hpp>
#include <stk_sddm/Value.hpp>

namespace stk {
namespace sddm {

class Property;
class AnyValue;
class PropertyRepository;
class Taxon;

typedef std::vector<Property *> PropertyVector;
typedef std::vector<const Property *> ConstPropertyVector;

/**
 * @brief Function <b>cast_error</b> creates a <b>bad_any_data_cast</b>, assembling a
 * message describing the name and invalid conversion.
 *
 * @param property	        a <b>Property</b> const ...
 *
 * @param of_type		a <b>AnyValue</b> const pointer
 *
 * @param to_type		a <b>std::type_info</b> const ...
 *
 * @return			a <b>bad_any_data_cast</b> ...
 */
std::runtime_error cast_error(const Property &property, const AnyValue *data, const char *to_type_name);


class Property
{
  friend class PropertyRepository;

public:
  typedef std::vector<Property *> Children;
  typedef Children::iterator iterator;
  typedef Children::const_iterator const_iterator;

private:
  Property(Property &parent, const std::string &name, AnyValue *any_value);

  Property(const std::string &name, AnyValue *any_value, const PropertyRepository &property_repository);

  Property(const Property &);
  Property &operator=(const Property &);

  Property &add(const std::string &name, AnyValue *any_value);

public:
  virtual ~Property();

  Property *getParent() const {
    return m_parent;
  }

  const PropertyRepository &getRepository() const {
    return m_propertyRepository;
  }

  const std::string &getName() const {
    return m_name;
  }

  const AnyValue *getAnyValue() const {
    return m_value;
  }

  Property &setAnyValue(AnyValue *any_value) {
    if (m_value)
      delete m_value;

    m_value = any_value;
    return *this;
  }

  const Taxon *getTaxon() const {
    return m_taxon;
  }

  Property &setTaxon(const Taxon *taxon) {
    m_taxon = taxon;
    return *this;
  }

  template <class T>
  Property &addAttribute(AttributeInterface<T> *attribute) {
    m_attributeMap.addAttribute(attribute);
    return *this;
  }

  template <class T>
  const T *findAttribute() const {
    return m_attributeMap.findAttribute<T>();
  }

  const AttributeMap &getAttributeMap() const {
    return m_attributeMap;
  }

  template <class T>
  const T &getAttribute() const {
    return *m_attributeMap.findAttribute<T>();
  }

  iterator begin() {
    return m_children.begin();
  }

  const_iterator begin() const  {
    return m_children.begin();
  }

  iterator end() {
    return m_children.end();
  }

  const_iterator end() const {
    return m_children.end();
  }

  Children::size_type size() const {
    return m_children.size();
  }

  template <typename T>
  bool isType() const {
    return m_value && typeid(T) == m_value->type();
  }

  const std::type_info &type() const {
    return m_value ? m_value->type() : typeid(void);
  }

  bool exists(const std::string &name) const {
    return find(name) != 0;
  }

  Property &create(const std::string &name) {
    return add(name, 0);
  }

  const Property *find(const std::string &name) const;

  Property *find(const std::string &name);

  const Property &get(const std::string &name) const;

  Property &get(const std::string &name);

  bool valueExists() const {
    return getAnyValue() != 0;
  }

  template <typename T>
  const T &value() const {
    if (m_value && m_value->type() == typeid(T)) {
      return m_value->value<T>();
    }
    else
      throw cast_error(*this, m_value, DataTypeTraits<T>::name());
  }

  template <typename T>
  const T &value(const std::string &name) const {
    return get(name).value<T>();
  }

  template <class T>
  Property &create(const std::string &name, const T &t) {
    return add(name, new Value<T>(t));
  }

  Property &create(const std::string &name, const char *s) {
    return create(name, std::string(s));
  }

  Property &create(const std::string &name, AnyValue *any_value) {
    return add(name, any_value);
  }

  void reportValidationError(const std::string &message) const;

  std::string path() const;

  virtual const Taxon *match(const std::string &name) const;

private:
  Property *                    m_parent;
  std::string                   m_name;
  AnyValue *                    m_value;
  const Taxon *                 m_taxon;
  const PropertyRepository &    m_propertyRepository;
  AttributeMap                  m_attributeMap;
  Children                      m_children;
};

} // namespace sddm
} // namespace stk

#endif // STK_SDDM_PROPERTY_HPP
