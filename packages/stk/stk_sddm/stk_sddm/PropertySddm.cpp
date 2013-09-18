#include <iomanip>
#include <sstream>
#include <iostream>
#include <algorithm>

#include <stk_sddm/ci_traits.hpp>
#include <stk_sddm/Property.hpp>
#include <stk_sddm/Taxonomy.hpp>
#include <stk_sddm/PropertyRepository.hpp>
#include <stk_sddm/SourceAttribute.hpp>

namespace stk {
namespace sddm {

namespace {

struct find_pred {
  find_pred(const std::string &name)
    : m_name(name)
  {}

  bool operator()(const Property *property) {
    const size_t property_name_size = property->getName().size();
    const size_t name_size = m_name.size();
    const size_t len = std::min(property_name_size, name_size);
    
    int r = ignorecase_traits::compare(property->getName().data(), m_name.data(), len);
    if (!r)
      r =  property_name_size - name_size;
    return r == 0;
  }
    
  std::string         m_name;
};
  

std::runtime_error
not_found_error(
  const Property &              property,
  const std::string &           name) 
{
  std::ostringstream oss;
  oss << "Property " << name << " is not in " << property.path();
  
  return std::runtime_error(oss.str());
}


#if 0
std::runtime_error
exists_error(
  const Property &              property,
  const std::string &           name) 
{
  std::ostringstream oss;
  oss << "Property " << name << " already exists in " << property.path();
  
  return std::runtime_error(oss.str());
}
#endif

} // namespace <unnamed>

std::runtime_error
state_error(
  const std::string &   state)
{
  std::ostringstream oss;

  oss << "State " << state << " not found";
  return std::runtime_error(oss.str());
}
  

std::runtime_error
cast_error(
  const Property &		property,
  const AnyValue *       	value,
  const char *                  to_type_name)
{
  std::ostringstream oss;
  oss << "Cannot cast property '" << property.path();

  if (value)
    oss << "' of type " << value->type().name() << " to type " << to_type_name;
  else
    oss << "' to type " << to_type_name << " because it has no value assigned";

  const SourceAttribute *source_attribute = property.findAttribute<SourceAttribute>();

  if (source_attribute)
    source_attribute->describe(oss);
  
  return std::runtime_error(oss.str());
}


/**
 * @brief Member function <b>text</b> attempts to cast the data to the
 * specified type.  If the data is not of the specified type, an exception is thrown.
 *
 * @return		a <b>T</b> reference to the data.
 *
 * @throws		a <b>std::runtime_error</b> exception with description.
 */


Property::Property(
  Property &                    parent,
  const std::string &           name,
  AnyValue *                    any_value)
  : m_parent(&parent),
    m_name(name),
    m_value(any_value),
    m_propertyRepository(parent.m_propertyRepository),
    m_attributeMap(),
    m_children()
{}


Property::Property(
  const std::string &           name,
  AnyValue *                    any_value,
  const PropertyRepository &    property_repository)
  : m_parent(0),
    m_name(name),
    m_value(any_value),
    m_taxon(0),
    m_propertyRepository(property_repository),
    m_attributeMap(),
    m_children()
{}


Property::~Property()
{
  delete m_value;
    
  for (Children::const_iterator it = m_children.begin(); it != m_children.end(); ++it)
    delete (*it);
}


const Taxon *
Property::match(
  const std::string &           name) const
{
  return m_taxon ? m_taxon->match(name) : 0;
}


Property &
Property::add(
  const std::string &           name,
  AnyValue *                    value)
{
  const Taxon *taxon = match(name);
  
  Property *new_property =  new Property(*this, name, value);
  new_property->setTaxon(taxon);

  m_children.push_back(new_property);
  
  if (taxon)
    taxon->validate(*new_property);
  
  return *m_children.back();
}


const Property *
Property::find(
  const std::string &   name) const
{
  Children::const_iterator it = std::find_if(m_children.begin(), m_children.end(), find_pred(name));
  if (it != m_children.end())
    return *it;
  else
    return 0;
}


Property *
Property::find(
  const std::string &   name)
{
  Children::iterator it = std::find_if(m_children.begin(), m_children.end(), find_pred(name));
  if (it != m_children.end())
    return *it;
  else
    return 0;
}


const Property &
Property::get(
  const std::string &   name) const
{
  const Property *p = find(name);
  if (!p)
    throw not_found_error(*this, name);

  return *p;
}


Property &
Property::get(
  const std::string &   name)
{
  Property *p = find(name);
  if (!p)
    throw not_found_error(*this, name);

  return *p;
}


void
Property::reportValidationError(
  const std::string &   message) const
{
  return m_propertyRepository.reportValidationError(*this, message);
}


std::string
Property::path() const
{
  std::string s;
  const Property *p = this;
  while (p) {
    std::string t;
    if (p->m_name.find('/') == std::string::npos && p->m_name.find(' ') == std::string::npos)
      t = std::string(p->m_name);
    else
      t = std::string("\"").append(p->m_name).append("\"");

    if (s.empty())
      s = t;
    else
      s = t.append("/").append(s);
      
    p = p->m_parent;
  }
  return s;
}

} // namespace sddm
} // namespace stk
