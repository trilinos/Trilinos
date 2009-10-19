/*    ------------------------------------------------------------
 *    Copyright 2002-2009 Sandia Corporation.
 *    Under the terms of Contract DE-AC04-94AL85000, there is a
 *    non-exclusive license for use of this work by or on behalf
 *    of the U.S. Government.  Export of this program may require
 *    a license from the United States Government.
 *    ------------------------------------------------------------
 */

#include <stk_util/environment/product_registry.h>
#include <stk_util/environment/ProductRegistry.hpp>
#include <stk_util/environment/sierra_version.hpp>

namespace stk {

const std::string
ProductRegistry::NAME = "Name";

const std::string
ProductRegistry::TITLE = "Title";

const std::string
ProductRegistry::VERSION = "Version";

const std::string
ProductRegistry::QUALIFIER = "Qualifier";

const std::string
ProductRegistry::BUILD_TIME = "Build Time";

const std::string
ProductRegistry::EXECUTABLE = "Executable";

const std::string
ProductRegistry::CONTACT = "Contact";

const std::string
ProductRegistry::ERROR = "Error";

const std::string
ProductRegistry::PRODUCT_TYPE = "Type";

const std::string
ProductRegistry::REGION_TITLE = "Region Title";

const std::string
ProductRegistry::BANNER_DETAIL = "Banner Detail";

const std::string
ProductRegistry::COPYRIGHT = "Copyright";

const std::string
ProductRegistry::PRODUCT_TYPE_REGION = "Region";


ProductRegistry &
ProductRegistry::instance()
{
  static ProductRegistry s_productRegistry;

  return s_productRegistry;
}


const char *
ProductRegistry::version()
{
  // SIERRA_VERSION should be a build-time define (i.e. -D flag) passed on
  // the compilation command line
  static const char *s_version = SIERRA_VERSION;
  
  return s_version;
}


ProductRegistry::AttributeMap &
ProductRegistry::addTPL(
  const std::string &   name,
  const std::string &   version,
  const std::string &   qualifier)
{
  std::pair<ProductMap::iterator, bool> iit = m_productMap.insert(std::make_pair(name, AttributeMap()));
  ProductMap::iterator it = iit.first;
  if (iit.second) {
    (*it).second[NAME] = name.c_str();
    (*it).second[VERSION] = version;
    (*it).second[QUALIFIER] = qualifier;
  }
  else {
    std::string &current_version = (*it).second[VERSION];
    std::string &current_qualifer = (*it).second[QUALIFIER];
    if (current_version.empty())
      current_version = version;
    if (current_qualifer.empty())
      current_qualifer = qualifier;
    if (current_version != version || current_qualifer != qualifier) {
      (*it).second[ERROR] = std::string("Product registration of ") + (*it).first + " version/qualifier conflict, "
        + " initially " + (*it).second[VERSION] + "/" + (*it).second[QUALIFIER]
        + " tried to change to " + version + "/" + qualifier;
      setRegistryInvalid();
    }
  }

  return (*it).second;
}


ProductRegistry::AttributeMap &
ProductRegistry::addProduct(const std::string &	name)
{
  std::pair<ProductMap::iterator, bool> iit = m_productMap.insert(std::make_pair(name, AttributeMap()));
  ProductMap::iterator it = iit.first;
  if (iit.second) {
    (*it).second[NAME] = name.c_str();
  }

  return (*it).second;
}


ProductRegistry::AttributeMap &
ProductRegistry::addRegion(
  const std::string &	name)
{
  AttributeMap &attribute_map = addProduct(name);
  attribute_map[ProductRegistry::PRODUCT_TYPE] = ProductRegistry::PRODUCT_TYPE_REGION;
  attribute_map[ProductRegistry::VERSION] = ProductRegistry::version();

  return attribute_map;
}


ProductRegistry::AttributeMap &
ProductRegistry::getProductAttributeMap(
  const std::string &	name)
{
  return m_productMap[name];
}


const std::string &
ProductRegistry::getProductAttribute(
  const std::string &	name,
  const std::string &	attribute) const
{
  return m_productMap[name][attribute];
}


std::string &
ProductRegistry::getProductAttribute(
  const std::string &	name,
  const std::string &	attribute)
{
  return m_productMap[name][attribute];
}


void
ProductRegistry::setProductAttribute(
  const std::string &	name,
  const std::string &	attribute,
  const std::string &	value)
{
  m_productMap[name][attribute] = value;
}

} // namespace stk

extern "C" {

void
product_registry_add(
  const char *		name )
{
  stk::ProductRegistry::instance().addProduct(name ? name : "<unknown>");
}


void
product_registry_add_tpl(
  const char *		name,
  const char *		version,
  const char *		qualifier )
{
  stk::ProductRegistry::instance().addTPL(name ? name : "<unknown>", version ? version : "", qualifier ? qualifier : "");
}


size_t
product_registry_size()
{
  return stk::ProductRegistry::instance().getProductMap().size();
}

} // extern "C"
