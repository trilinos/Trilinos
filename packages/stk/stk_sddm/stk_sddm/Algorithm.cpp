#include <stk_sddm/Algorithm.hpp>

#include <algorithm>

namespace stk {
namespace sddm {

void
remove(
  PropertyVector &              property_vector,
  const std::string &           name)
{
  PropertyVector::iterator it = std::remove_if(property_vector.begin(), property_vector.end(), find_pred(name));
  property_vector.erase(it, property_vector.end());
}


void
remove(
  ConstPropertyVector &         property_vector,
  const std::string &           name)
{
  ConstPropertyVector::iterator it = std::remove_if(property_vector.begin(), property_vector.end(), find_pred(name));
  property_vector.erase(it, property_vector.end());
}


bool
isDefaultProperty(
  const Property &      parent_property,
  const std::string &   name) 
{
  const Property *property = parent_property.find(name);

  return property && property->findAttribute<DefaultAttribute>();
}

} // namespace sddm
} // namespace stk
