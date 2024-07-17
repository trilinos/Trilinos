// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_CLOSURE_MODEL_FACTORY_BASE_HPP
#define PANZER_CLOSURE_MODEL_FACTORY_BASE_HPP

#include "PanzerDiscFE_config.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Phalanx_Evaluator.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_Base.hpp"
#include "Panzer_EvaluatorsRegistrar.hpp"
#include <string>
#include <vector>

namespace panzer {

  class FieldLayoutLibrary;
  class IntegrationRule;
  struct GlobalData;

  class ClosureModelFactoryBase : public panzer::Base,
                                  public EvaluatorsRegistrar {

  public:

    ClosureModelFactoryBase() {}
    
    virtual ~ClosureModelFactoryBase() {}
    
    Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >
    virtual  buildClosureModels(const std::string& model_id,
                                const Teuchos::ParameterList& models,
                                const panzer::FieldLayoutLibrary& fl,
                                const Teuchos::RCP<panzer::IntegrationRule>& ir,
                                const Teuchos::ParameterList& equation_set_params,
                                const Teuchos::ParameterList& user_data,
                                const Teuchos::RCP<panzer::GlobalData>& global_data,
                                PHX::FieldManager<panzer::Traits>& fm) const = 0;

    /** This a convenience function for registering the evaluators. Essentially this
      * facilitates better usage of the ClosureModel TM and allows an easy registration
      * process externally without knowning the compile-time evaluation type.
      *
      * \param[in] evaluators Evaluators to register
      * \param[in] fm Field manager where the evaluators will be registered on completion.
      */
    virtual void registerEvaluators(const std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > & evaluators,
                                    PHX::FieldManager<panzer::Traits>& fm) const = 0;

  };
  
}

#endif
