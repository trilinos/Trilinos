// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_BCSTRATEGY_BASE_HPP
#define PANZER_BCSTRATEGY_BASE_HPP

#include "Panzer_Traits.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "Panzer_EvaluatorsRegistrar.hpp"

namespace PHX {
  template<typename T> class FieldManager;
}

namespace panzer {

  class PhysicsBlock;

  //! Non-templated empty base class for BCStrategy objects
class BCStrategyBase : public EvaluatorsRegistrar {
    
  public:
    
    BCStrategyBase() {}
    
    virtual ~BCStrategyBase() {}
    
    virtual void 
    setup(const panzer::PhysicsBlock& side_pb,
	  const Teuchos::ParameterList& user_data) = 0;

    virtual void 
      buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				 const panzer::PhysicsBlock& side_pb,
				 const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
				 const Teuchos::ParameterList& models,
				 const Teuchos::ParameterList& user_data) const = 0;

    virtual void 
    buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				      const panzer::PhysicsBlock& side_pb,
				      const LinearObjFactory<panzer::Traits> & lof,
				      const Teuchos::ParameterList& user_data) const = 0;

    virtual void
    buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					           const panzer::PhysicsBlock& side_pb,
						   const LinearObjFactory<panzer::Traits> & lof,
						   const Teuchos::ParameterList& user_data) const = 0;

  };
  
}

#endif
