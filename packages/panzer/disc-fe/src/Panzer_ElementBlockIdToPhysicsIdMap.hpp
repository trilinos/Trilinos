// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_ELEMENT_BLOCK_ID_TO_PHYSICS_ID_MAP_HPP
#define PANZER_ELEMENT_BLOCK_ID_TO_PHYSICS_ID_MAP_HPP

#include <map>
#include <string>

namespace Teuchos {
  class ParameterList;
}

namespace panzer {

  /** \brief Populates a map from element block ID to physics ID,
    * one entry per parameter from the given flat ParameterList,
    * where each parameter's name is the element block ID and its
    * string value is the physics ID.
    *
    * \param[out] b_to_p the map to populate; existing entries for block IDs also present in p are overwritten.
    * \param[in] p a flat ParameterList of "<element block ID>" = "<physics ID>" entries.
    */
  void buildBlockIdToPhysicsIdMap(std::map<std::string,std::string>& b_to_p,
				                          const Teuchos::ParameterList& p);

}

#endif
