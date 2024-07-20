// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_ElementBlockIdToPhysicsIdMap.hpp"
#include "Teuchos_ParameterList.hpp"

namespace panzer {

  void buildBlockIdToPhysicsIdMap(std::map<std::string,std::string>& b_to_p,
				  const Teuchos::ParameterList& p)
  {
    for (Teuchos::ParameterList::ConstIterator entry = p.begin();
	 entry != p.end(); ++entry) {
      std::string dummy_type;
      b_to_p[entry->first] = entry->second.getValue(&dummy_type);
    }
  }

}
