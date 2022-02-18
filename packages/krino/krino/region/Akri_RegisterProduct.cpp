// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_RegisterProduct.hpp>

#include <stk_util/registry/ProductRegistry.hpp>

namespace krino {

const char *
get_product_name()
{ /* %TRACE% */  /* %TRACE% */
  return "Krino";
}


void
register_product()
{ /* %TRACE% */  /* %TRACE% */
  // Register krino
  stk::ProductRegistry::instance().addRegion(get_product_name());
}

} // namespace krino

