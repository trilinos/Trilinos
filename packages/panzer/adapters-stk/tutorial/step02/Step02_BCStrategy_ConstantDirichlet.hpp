// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Step02_BCStrategy_ConstantDirichlet_hpp__
#define __Step02_BCStrategy_ConstantDirichlet_hpp__

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_BCStrategy_Dirichlet_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_PureBasis.hpp"
#include "Phalanx_FieldManager.hpp"

namespace user_app {

  template <typename EvalT>
  class BCStrategy_ConstantDirichlet : public panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT> {
  public:    
    
    BCStrategy_ConstantDirichlet(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data);
    
    void setup(const panzer::PhysicsBlock& side_pb,
	       const Teuchos::ParameterList& user_data);
    
    void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				    const panzer::PhysicsBlock& pb,
				    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
				    const Teuchos::ParameterList& models,
				    const Teuchos::ParameterList& user_data) const;

    std::string residual_name;
    Teuchos::RCP<panzer::PureBasis> basis;
  };

}

#include "Step02_BCStrategy_ConstantDirichlet_impl.hpp"

#endif
