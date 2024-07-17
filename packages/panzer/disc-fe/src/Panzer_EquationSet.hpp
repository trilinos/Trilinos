// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EQUATION_SET_HPP
#define PANZER_EQUATION_SET_HPP

#include "Panzer_EquationSet_Base.hpp"

namespace PHX {
  template<typename Traits> class FieldManager;
}

namespace panzer {

  template <typename EvalT>
  class EquationSet : public panzer::EquationSetBase {
    
  public:    
    
    EquationSet() {}

    virtual ~EquationSet() {}

    /// \name Initialization (derived from panzer::EquationSetBase)
    ///@{ 

    virtual void setElementBlockId(const std::string & blockId) = 0;

    ///@}

    /// \name Evaluator Construction and Registration Methods (derived from panzer::EquationSetBase)
    ///@{

    virtual void buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                                                const panzer::FieldLibrary& field_library,
                                                                const LinearObjFactory<panzer::Traits> & lof,
                                                                const Teuchos::ParameterList& user_data) const = 0;
    
    virtual void buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                                   const panzer::FieldLibrary& field_library,
                                                   const LinearObjFactory<panzer::Traits> & lof,
                                                   const Teuchos::ParameterList& user_data) const = 0;

    virtual void buildAndRegisterDOFProjectionsToIPEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                                              const panzer::FieldLayoutLibrary& field_library,
                                                              const Teuchos::RCP<panzer::IntegrationRule>& ir,
                                                              const Teuchos::Ptr<const panzer::LinearObjFactory<panzer::Traits> > & lof,
                                                              const Teuchos::ParameterList& user_data) const = 0;
    
    virtual void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                                       const panzer::FieldLibrary& field_library,
                                                       const Teuchos::ParameterList& user_data) const = 0;

    virtual void buildAndRegisterClosureModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                                              const panzer::FieldLayoutLibrary& field_library,
                                                        const Teuchos::RCP<panzer::IntegrationRule>& ir,
                                                        const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                                        const Teuchos::ParameterList& models,
                                                        const Teuchos::ParameterList& user_data) const = 0;

    virtual void buildAndRegisterClosureModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                                        const panzer::FieldLayoutLibrary& field_library,
                                                        const Teuchos::RCP<panzer::IntegrationRule>& ir,
                                                        const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                                        const std::string& model_name,
                                                        const Teuchos::ParameterList& models,
                                                        const Teuchos::ParameterList& user_data) const = 0;
    
    virtual void buildAndRegisterInitialConditionEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                                            const panzer::FieldLibrary& field_library,
                                                            const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                                            const std::string& model_name,
                                                            const Teuchos::ParameterList& models,
                                                            const LinearObjFactory<panzer::Traits> & lof,
                                                            const Teuchos::ParameterList& user_data) const = 0;

    ///@}

    /// \name Query Methods (derived from panzer::EquationSetBase)
    ///@{

    virtual const Teuchos::RCP<Teuchos::ParameterList> getEvaluatorParameterList() const = 0;
    
    virtual const std::vector<std::pair<std::string,Teuchos::RCP<panzer::PureBasis> > > & getProvidedDOFs() const = 0;

    virtual const std::vector<std::vector<std::string> > & getCoordinateDOFs() const = 0;

    virtual const std::map<int,Teuchos::RCP<panzer::IntegrationRule> > & getIntegrationRules() const = 0;

    virtual std::string getElementBlockId() const = 0;

    virtual std::string getType() const = 0;

    ///@}

  };
  
}

#endif
