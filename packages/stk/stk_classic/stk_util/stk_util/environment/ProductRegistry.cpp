// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <stk_util/environment/product_registry.h>
#include <stk_util/environment/ProductRegistry.hpp>
#ifdef STK_BUILT_IN_SIERRA
#include <stk_util/environment/stk_version.hpp>
#else
#define STK_VERSION "0.1a"
#endif

namespace stk_classic {

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
  // STK_VERSION should be a build-time define (i.e. -D flag) passed on
  // the compilation command line
  static const char *s_version = STK_VERSION;
  
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

} // namespace stk_classic

extern "C" {

void
product_registry_add(
  const char *		name )
{
  stk_classic::ProductRegistry::instance().addProduct(name ? name : "<unknown>");
}

void
product_registry_add_tpl(
  const char *		name,
  const char *		version,
  const char *		qualifier )
{
  stk_classic::ProductRegistry::instance().addTPL(name ? name : "<unknown>", version ? version : "", qualifier ? qualifier : "");
}

} // extern "C"
