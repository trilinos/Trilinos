/**   ------------------------------------------------------------
 *    Copyright 2002-2009 Sandia Corporation.
 *    Under the terms of Contract DE-AC04-94AL85000, there is a
 *    non-exclusive license for use of this work by or on behalf
 *    of the U.S. Government.  Export of this program may require
 *    a license from the United States Government.
 *    ------------------------------------------------------------
 */

#include <stk_util/environment/ProductRegistry.hpp>

#include <stk_util/environment/RegisterProduct.hpp>

namespace sierra {

const char *
get_product_name()
{
  return "UtilityLib";
}

void
register_product()
{
  // Register utility
  ProductRegistry::AttributeMap &attr_map = ProductRegistry::instance().addProduct(get_product_name());
  attr_map[ProductRegistry::VERSION]      = ProductRegistry::version();
  attr_map[ProductRegistry::TITLE] = "Utility library routines";
  attr_map[ProductRegistry::CONTACT] = "framework-developers@sourceforge.sandia.gov";

  // Register TPL's and other things which may not be properly registered but used directly.

}

} // namespace sierra
