// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#ifndef STK_UTIL_REGISTRY_PRODUCTREGISTRY_HPP
#define STK_UTIL_REGISTRY_PRODUCTREGISTRY_HPP

#include <map>     // for map, map<>::const_iterator, map<>::value_compare
#include <string>  // for string

namespace stk {

/**
 * @brief Class <b>ProductRegistry</b> maps product names and attributes to
 * brief descriptive values.  Each added product has at a minimum the <b>NAME</b>
 * value initialized.  A region type product also has the <b>PRODUCT_TYPE</b> set to
 * <b>PRODUCT_TYPE_REGION</b>.
 *
 * A product may have any other attribute keyword and value assocaited with it.  In
 * particular, the main execution product often has the <b>/EXECUTABLE_FILE</b> and
 * <b>BUILD_TIME</b> attributes set.
 *
 * You may add any other attribute keyword and values to the product registry.
 */
class ProductRegistry
{
public:
  typedef std::map<std::string, std::string> AttributeMap;	///< Map of attribute keyword to value
  typedef std::map<std::string, AttributeMap> ProductMap;	///< Map of product name to attribute map

  /**
   * @brief Member function <b>instance</b> returns a reference to the registry singleton.
   *
   * @return			a <b>ProductRegistry</b> reference to the registry singleton
   */
  static ProductRegistry &instance();

private:
  /**
   * Creates a new <b>ProductRegistry</b> instance.
   */
  ProductRegistry()
    : m_productMap(),
      m_productName(),
      m_registryOK(true)
  {}

  ProductRegistry(const ProductRegistry &);

  ProductRegistry &operator=(const ProductRegistry &);

  ~ProductRegistry()
  {}

public:
  /**
   * @brief Member function <b>version</b> returns the version number of the combined product.
   *
   * @return			a <b>ProductRegistry</b> reference to the registry singleton
   */
  static const char *version();

  /**
   * @brief Member function <b>setRegistryInvalid</b> marks th registry as contain a
   * conflict of some sort.  The <b>ERROR</b> attribute of each product describes
   * each error.
   */
  void setRegistryInvalid() {
    m_registryOK = false;
  }

  /**
   * @ingroup EnvRuntimeInformationDetail
   * @brief Function <b>executable_date</b> returns the build date of the executable
   * file as a string in the form Mmm dd yyyy hh:mm::ss.
   *
   * @return			a <b>std::string</b> const reference to the executable
   *				file's build date and time.
   */
  static const std::string &executable_date();

  /**
   * @brief Member function <b>isRegistryOK</b> returns true if the registry has not
   * been flagged as having an error via <b>setRegistryInvalid</b>.
   *
   * @return			a <b>bool</b> value of true of the registry has no
   *				conflicts.
   */
  bool isRegistryOK() {
    return m_registryOK;
  }

  /**
   * @brief Member function <b>getProductName</b> returns the product name.
   *
   * @return			a <b>string</b> const reference to the product name.
   */
  const std::string &getProductName() const {
    return m_productName;
  }

  /**
   * @brief Member function <b>setProductName</b> sets the product name.
   */
  void setProductName(const std::string &product_name) {
    m_productName = product_name;
  }

  const ProductMap& productMap() const { return m_productMap; }
        ProductMap& productMap()       { return m_productMap; }

  ProductMap::const_iterator productMap_begin() const {
      return m_productMap.begin();
  }
  ProductMap::iterator productMap_begin() {
      return m_productMap.begin();
  }

  ProductMap::const_iterator productMap_end() const {
      return m_productMap.end();
  }
  ProductMap::iterator productMap_end() {
      return m_productMap.end();
  }

  /**
   * @brief Member function <b>addProduct</b> adds a product to the registry.
   *
   * @param name		a <b>std::string</b> const reference to the product's
   *				name.
   *
   * @return			an <b>AttributeMap</b> reference to the attribute
   *				map of the newly added product.
   */
  AttributeMap &addProduct(const std::string &name);

  /**
   * @brief Member function <b>addTPL</b> adds a product to the registry.  A
   * product always has the <b>VERSION</b> and <b>QUALIFIER</b> attributes set.
   *
   * @param name                a <b>std::string</b> const reference to the product's
   *                            name.
   *
   * @param version             a <b>std::string</b> const reference to the product's
   *                            version string.
   *
   * @param qualifier           a <b>std::string</b> const reference to the product's
   *                            qualifier string.
   *
   * @return                    an <b>AttributeMap</b> reference to the attribute
   *                            map of the newly added product.
   */
  AttributeMap &addTPL(const std::string &name, const std::string &version, const std::string &qualifier = "");

  /**
   * @brief Member function <b>addRegion</b> add a region as a product to the
   * registry.  A region product is a product with the <b>PRODUCT_TYPE</b> set to
   * <b>PRODUCT_TYPE_REGION</b>.
   *
   * @param name		a <b>std::string</b> const reference to the product's
   *				name.
   *
   * @return			an <b>AttributeMap</b> reference to the attribute
   *				map of the newly added region product.
   */
  AttributeMap &addRegion(const std::string &name);

  /**
   * @brief Member function <b>getProduct</b> returns a reference to the product
   * attribute map.
   *
   * @param name		a <b>std::string</b> const reference to the product's
   *				name.
   *
   * @return			an <b>AttributeMap</b> reference to the attribute
   *				map of the newly added region product.
   */
  AttributeMap &getProductAttributeMap(const std::string &name);

  /**
   * @brief Member function <b>getAttribute</b> returns the
   * attribute for the named product.
   *
   * @param name		a <b>std::string</b> const reference to the product's
   *				name.
   *
   * @param attribute		a <b>std::string</b> const reference to the product's
   *				attribute keyword.
   *
   * @return			a <b>std::string</b> const reference to the value of
   *				product's attribute.
   */
  std::string &getProductAttribute(const std::string &name, const std::string &attribute);

public:
  static const std::string	NAME;
  static const std::string	TITLE;
  static const std::string  VERSION;
  static const std::string  QUALIFIER;
  static const std::string	CONTACT;
  static const std::string	ERROR;
  static const std::string	PRODUCT_TYPE;

  static const std::string	EXECUTABLE;
  static const std::string	BUILD_TIME;

  static const std::string	BANNER_DETAIL;
  static const std::string	COPYRIGHT;

  static const std::string	REGION_TITLE;

  static const std::string	PRODUCT_TYPE_REGION;

private:
  mutable ProductMap m_productMap;
  std::string        m_productName;
  bool				 m_registryOK;
};

std::string get_version(const std::string& executableName);

} // namespace stk

#ifdef STK_BUILT_FOR_SIERRA
#undef VERSION // Nice, Trilinos leaves us this gem

namespace sierra {

typedef stk::ProductRegistry ProductRegistry;

} // namespace sierra
#endif // STK_BUILT_FOR_SIERRA


#endif // STK_UTIL_REGISTRY_PRODUCTREGISTRY_HPP
