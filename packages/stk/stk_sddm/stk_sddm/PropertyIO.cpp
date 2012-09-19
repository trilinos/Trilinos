#include <iomanip>
#include <sstream>

#include <stk_sddm/Taxonomy.hpp>
#include <stk_sddm/AttributeMap.hpp>
#include <stk_sddm/PropertyIO.hpp>


namespace stk {
namespace sddm {

namespace {

std::string
xml_clean(
  const std::string &s)
{
  std::string t;

  for (std::string::const_iterator c = s.begin(); c != s.end(); ++c)
    if (*c == '<')
      t += "&lt;";
    else if (*c == '>')
      t += "&gt;";
    else if (*c == '"')
      t += "&quot;";
    else
      t += *c;
  return t;
}


std::ostream &
print_property(
  std::ostream &                os,
  int                           depth,
  bool                          feature_coverage, 
  std::vector<const Property *> active_properties,
  const Property &              property)
{
  const AttributeMap &attribute_map = property.getAttributeMap();
  const AnyValue *any_value = property.getAnyValue();

  active_properties.push_back(&property);

  if (depth != 0)
    os << std::setw(depth*2) << "";

  if (feature_coverage) {
    for (std::vector<const Property *>::const_iterator it = active_properties.begin(); it != active_properties.end(); ++it)
      if ((*it)->getTaxon())
        os << (*it)->getTaxon()->getId() << " ";
    os << ", ";
  
    for (std::vector<const Property *>::const_iterator it = active_properties.begin(); it != active_properties.end(); ++it)
      if ((*it)->getTaxon())
        os << (*it)->getTaxon()->getName() << " ";
    os << ", ";
  }
  
  if (any_value)
    os << property.getName() << "(" << any_value->typeName() << "): " << *any_value;
  else
    os << property.getName();

  for (AttributeMap::const_iterator it = attribute_map.begin(); it != attribute_map.end(); ++it) {
    os << " [";
    (*it).second->dump(os);
    os << "]";
  }
  
  os << std::endl;
  
  for (Property::Children::const_iterator it = property.begin(); it != property.end(); ++it)
    print_property(os, depth + 1, feature_coverage, active_properties, *(*it));

  active_properties.pop_back();
  
  return os;
}

} // namespace <unnamed>

std::ostream &
print_property(
  std::ostream &        os,
  bool                  feature_coverage, 
  const Property &      property)
{  
  std::vector<const Property *> active_properties;
  
  print_property(os, 0, feature_coverage, active_properties, property);

  return os;
}


std::ostream &
print_xml(
  std::ostream &        os,
  const Property &      property,
  size_t                indent)
{
  os << std::setw(indent*2) << "" << "<Property>" << std::endl;
  ++indent;

  os << std::setw(indent*2) << "" << "<Name>" << xml_clean(property.getName()) << "</Name>" << std::endl;
  if (property.getTaxon()) {
    os << std::setw(indent*2) << "" << "<TaxonId>" << xml_clean(property.getTaxon()->getId()) << "</TaxonId>" << std::endl;
    os << std::setw(indent*2) << "" << "<Taxon>" << xml_clean(property.getTaxon()->getName()) << "</Taxon>" << std::endl;
  }

  for (AttributeMap::const_iterator it = property.getAttributeMap().begin(); it != property.getAttributeMap().end(); ++it)
    (*it).second->xml(os, indent);

  if (property.getAnyValue()) {
//    std::ostringstream oss;
    property.getAnyValue()->xml(os, indent);
      
//    os << std::setw(indent*2) << "" << "<Type>" << xml_clean(property.getAnyValue()->typeName()) << "</Type>" << std::endl;
//    os << std::setw(indent*2) << "" << "<Value>" << xml_clean(oss.str()) << "</Value>" << std::endl;
//    os << std::setw(indent*2) << "" << oss.str() << std::endl;
  }
    
  for (Property::Children::const_iterator it = property.begin(); it != property.end(); ++it)
    print_xml(os, *(*it), indent + 1);

  --indent;
  os << std::setw(indent*2) << "" << "</Property>" << std::endl;
  
  return os;
}


std::ostream &
operator<<(
  std::ostream &        os,
  const Property &      property)
{
  print_property(os, false, property);
  
  return os;
}

} // namespace sddm
} // namespace stk
