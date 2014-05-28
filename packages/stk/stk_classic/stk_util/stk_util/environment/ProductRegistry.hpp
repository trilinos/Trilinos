/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_UTIL_ENVIRONMENT_PRODUCTREGISTRY_HPP
#define STK_UTIL_ENVIRONMENT_PRODUCTREGISTRY_HPP

#include <map>
#include <string>

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
 *
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
   *
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
   *
   */
  void setRegistryInvalid() {
    m_registryOK = false;
  }

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
   *
   */
  void setProductName(const std::string &product_name) {
    m_productName = product_name;
  }

  /**
   * @brief Member function <b>getProductMap</b> returns a reference to the map of
   * all products.
   *
   * @return			a <b>ProductMap</b> reference to the product map.
   */
  ProductMap &getProductMap() const {
    return m_productMap;
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
   *
   */
  const std::string &getProductAttribute(const std::string &name, const std::string &attribute) const;

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
   *
   */
  std::string &getProductAttribute(const std::string &name, const std::string &attribute);

  /**
   * @brief Member function <b>setAttribute</b> sets the
   * attribute for the named poduct to the specified value.
   *
   * @param name		a <b>std::string</b> const reference to the product's
   *				name.
   *
   * @param attribute		a <b>std::string</b> const reference to the product's
   *				attribute keyword.
   *
   * @param value		a <b>std::string</b> const reference to the value of
   *				product's attribute.
   *
   */
  void setProductAttribute(const std::string &name, const std::string &attribute, const std::string &value);

public:
  static const std::string	NAME;				///< Product's name attribute
  static const std::string	TITLE;				///< Product's title attribute
  static const std::string      VERSION;                        ///< TPL's version attribute
  static const std::string      QUALIFIER;                      ///< TPL's qualifier attribute
  static const std::string	CONTACT;			///< Product's contact attribute
  static const std::string	ERROR;				///< Product's error attribute
  static const std::string	PRODUCT_TYPE;			///< Product's product_type attribute

  static const std::string	EXECUTABLE;			///< Product's executable attribute
  static const std::string	BUILD_TIME;			///< Product's build_time attribute

  static const std::string	BANNER_DETAIL;			///< Product's additional banner info attribute
  static const std::string	COPYRIGHT;			///< Product's copyright information

  static const std::string	REGION_TITLE;			///< Product's region_title attribute

  static const std::string	PRODUCT_TYPE_REGION;		///< Region product_type value

private:
  mutable ProductMap		m_productMap;			///< Product map
  std::string                   m_productName;                  ///< Name of main product
  bool				m_registryOK;			///< Registry is OK
};

} // namespace stk

#ifdef STK_BUILT_IN_SIERRA
#undef VERSION // Nice, Trilinos leaves us this gem

namespace sierra {

typedef stk::ProductRegistry ProductRegistry;

} // namespace sierra
#endif // STK_BUILT_IN_SIERRA


#endif // STK_UTIL_ENVIRONMENT_PRODUCTREGISTRY_HPP
