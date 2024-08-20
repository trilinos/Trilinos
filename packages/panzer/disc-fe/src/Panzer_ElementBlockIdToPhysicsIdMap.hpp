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

  void buildBlockIdToPhysicsIdMap(std::map<std::string,std::string>& b_to_p, 
				  const Teuchos::ParameterList& p);

}

#endif
