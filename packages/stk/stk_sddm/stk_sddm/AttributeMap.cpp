#include <stk_sddm/AttributeMap.hpp>

namespace stk {
namespace sddm {

AttributeMap::~AttributeMap()
{
  for(Map::iterator it=m_attributeMap.begin(); it!=m_attributeMap.end(); ++it) {
    delete it->second;
  }
}

void
AttributeMap::addAttribute(
  const std::type_info *        type,
  AttributeInterfaceBase *      attribute)
{
  delete m_attributeMap[type];

  m_attributeMap[type] = attribute;
}

} // namespace sddm
} // namespace stk
