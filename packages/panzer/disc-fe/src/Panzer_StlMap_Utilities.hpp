// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_STLMAP_UTILITIES_HPP
#define PANZER_STLMAP_UTILITIES_HPP

#include "Teuchos_Assert.hpp"
#include <map>
#include <iostream>

namespace panzer {

  template<typename MapT>
  const typename MapT::mapped_type &
  getEntry(const MapT& in_map, const typename MapT::key_type& in_key) {
    
    typename MapT::const_iterator it = in_map.find(in_key);

    TEUCHOS_TEST_FOR_EXCEPTION(it == in_map.end(), std::runtime_error,
		       "Failed to find the key " << in_key << " in the map."
		       << std::endl);

    return it->second;
  }

}

#endif
