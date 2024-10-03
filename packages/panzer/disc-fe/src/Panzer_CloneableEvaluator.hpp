// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_CloneableEvaluator_H
#define PANZER_CloneableEvaluator_H

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace panzer {

  //! Non-templated empty base class for template managers
  class CloneableEvaluator {
    
  public:
    
    CloneableEvaluator() {}
    
    virtual ~CloneableEvaluator() {}
    
    virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const = 0;
  };
  
}

#endif
