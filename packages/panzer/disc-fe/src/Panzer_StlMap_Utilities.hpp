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

  /** \brief Looks up a key in an STL-like associative container,
    * throwing std::runtime_error if the key is not found (instead of
    * the undefined behavior of operator[] on a const map, or the
    * silent insertion of operator[] on a non-const map).
    *
    * \tparam MapT an STL-like associative container type (e.g. std::map).
    * \param in_map the container to search.
    * \param in_key the key to look up.
    * \return the value associated with in_key.
    */
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
