// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EQUATION_SET_BASE_HPP
#define PANZER_EQUATION_SET_BASE_HPP

#include <string>
#include <vector>
#include <map>
#include <utility>
#include "Teuchos_RCP.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_FieldLibrary.hpp"
#include "Panzer_EvaluatorsRegistrar.hpp"

namespace Teuchos {
  class ParameterList;
}

namespace panzer {

  class PureBasis;
  class IntegrationRule;
  class FieldLibrary;
  class FieldLayoutLibrary;

  //! Non-templated empty base class for EquationSet objects
  class EquationSetBase : public EvaluatorsRegistrar {
    
  public:
    
    EquationSetBase() {}
    
    virtual ~EquationSetBase() {}
    
    /// \name Initialization
    ///@{ 

    virtual void setElementBlockId(const std::string & blockId) = 0;

    ///@}

    /// \name Evaluator Construction and Registration Methods
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
    
    //! Register closure model evaluators with the model name internally specified by the equation set
    virtual void buildAndRegisterClosureModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
							const panzer::FieldLayoutLibrary& field_library,
							const Teuchos::RCP<panzer::IntegrationRule>& ir,
							const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
							const Teuchos::ParameterList& models,
							const Teuchos::ParameterList& user_data) const = 0;

    //! Register closure model evaluators with the model name specified by an argument
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
							    const panzer::LinearObjFactory<panzer::Traits> & lof,
							    const Teuchos::ParameterList& user_data) const = 0;

    ///@}

    /// \name Query Methods
    ///@{

    //! Returns the parameter list that will be passed off from the equaiton set to the closure model evaluator factory.  This allows users to pass parameters from a particular equaiton set to its associated closure models.
    virtual const Teuchos::RCP<Teuchos::ParameterList> getEvaluatorParameterList() const = 0;
    
    //! Return the Basis for the equation set, key is the DOF name (note coordinate DOFs are NOT included)
    virtual const std::vector<std::pair<std::string,Teuchos::RCP<panzer::PureBasis> > > & getProvidedDOFs() const = 0;

    //! Return a vector of vectors that correspond to DOFs set as coordinate fields
    virtual const std::vector<std::vector<std::string> > & getCoordinateDOFs() const = 0;

    //! Return a map of unique integration rules for the equation set, key is the integration order
    virtual const std::map<int,Teuchos::RCP<panzer::IntegrationRule> > & getIntegrationRules() const = 0;

    virtual std::string getElementBlockId() const = 0;

    //! Returns the type of the equation set object.  Corresponds to the keyword used by the equation set factory to build a particular concrete equation set.
    virtual std::string getType() const = 0;

    ///@}

    //! Set the list of tangent parameter names
    virtual void setTangentParamNames(const std::vector<std::string>& tangent_param_names) = 0;

  };
  
}

#endif
