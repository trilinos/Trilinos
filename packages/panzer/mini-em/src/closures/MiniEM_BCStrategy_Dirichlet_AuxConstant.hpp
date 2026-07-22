// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MiniEM_BCSTRATEGY_DIRICHLET_AUXCONSTANT_DECL_HPP
#define MiniEM_BCSTRATEGY_DIRICHLET_AUXCONSTANT_DECL_HPP

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_BCStrategy.hpp"
#include "Panzer_Traits.hpp"
#include "Phalanx_FieldManager.hpp"

namespace mini_em {

  template <typename EvalT>
  class BCStrategy_Dirichlet_AuxConstant : public panzer::BCStrategy<EvalT> {
    
  public:    
    
    BCStrategy_Dirichlet_AuxConstant(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& /* global_data */);
    
    void setup(const panzer::PhysicsBlock& side_pb,
	       const Teuchos::ParameterList& /* user_data */);
    
    void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				    const panzer::PhysicsBlock& /* pb */,
				    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& /* factory */,
				    const Teuchos::ParameterList& /* models */,
				    const Teuchos::ParameterList& /* user_data */) const;

    void buildAndRegisterGatherScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                                 const panzer::PhysicsBlock& side_pb,
                                                 const panzer::LinearObjFactory<panzer::Traits> & lof,
                                                 const Teuchos::ParameterList& user_data) const;

    virtual void 
    buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				      const panzer::PhysicsBlock& side_pb,
				      const panzer::LinearObjFactory<panzer::Traits> & lof,
				      const Teuchos::ParameterList& user_data) const;

    virtual void
    buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					           const panzer::PhysicsBlock& side_pb,
						   const panzer::LinearObjFactory<panzer::Traits> & lof,
						   const Teuchos::ParameterList& user_data) const;

    std::string operatorName_;
    double value_;
    std::string fieldName_;
    Teuchos::RCP<panzer::PureBasis> basis_;

  };

template < >
void BCStrategy_Dirichlet_AuxConstant<panzer::Traits::Jacobian>::
buildAndRegisterGatherScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                        const panzer::PhysicsBlock& side_pb,
                                        const panzer::LinearObjFactory<panzer::Traits> & lof,
                                        const Teuchos::ParameterList& user_data) const;

template < >
void BCStrategy_Dirichlet_AuxConstant<panzer::Traits::Jacobian>::
buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                        const panzer::PhysicsBlock& side_pb,
                                        const panzer::LinearObjFactory<panzer::Traits> & lof,
                                        const Teuchos::ParameterList& user_data) const;
template < >
void BCStrategy_Dirichlet_AuxConstant<panzer::Traits::Jacobian>::
buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                        const panzer::PhysicsBlock& side_pb,
                                        const panzer::LinearObjFactory<panzer::Traits> & lof,
                                        const Teuchos::ParameterList& user_data) const;

}

#include "MiniEM_BCStrategy_Dirichlet_AuxConstant_impl.hpp"

#endif
