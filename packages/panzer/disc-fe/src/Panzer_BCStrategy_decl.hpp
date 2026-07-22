// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_BCSTRATEGY_DECL_HPP
#define PANZER_BCSTRATEGY_DECL_HPP

#include "Panzer_BCStrategy_Base.hpp"
#include "Panzer_BC.hpp"

namespace panzer {

  template <typename EvalT>
  class BCStrategy : public panzer::BCStrategyBase {
    
  public:    
    
    BCStrategy(const panzer::BC& bc);
    
    virtual ~BCStrategy();
    
    /** Must be called before this->buildAndRegisterEvaluators() and
	this->buildAndRegisterGatherScatterEvaluators().
    */
    virtual void setup(const panzer::PhysicsBlock& side_pb,
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
    
  protected:
    
    const panzer::BC m_bc;

  };
  
}

#endif
