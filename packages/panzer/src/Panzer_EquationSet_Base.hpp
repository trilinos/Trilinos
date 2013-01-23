// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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

namespace Teuchos {
  class ParameterList;
}

namespace panzer {

  class PureBasis;
  class IntegrationRule;
  class BasisIRLayout;
  class FieldLibrary;
  class FieldLayoutLibrary;

  //! Non-templated empty base class for EquationSet objects
  class EquationSetBase {
    
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

    virtual const Teuchos::RCP<Teuchos::ParameterList> getEvaluatorParameterList() const = 0;
    
    //! Return the Basis for the equation set, key is the DOF name
    virtual const std::vector<std::pair<std::string,Teuchos::RCP<panzer::PureBasis> > > & getProvidedDOFs() const = 0;

    //! Return a map of unique integration rules for the equation set, key is the integration order
    virtual const std::map<int,Teuchos::RCP<panzer::IntegrationRule> > & getIntegrationRules() const = 0;

    virtual std::string getElementBlockId() const = 0;

    //! Returns the type of the equation set object.  Corresponds to the keyword used by the equation set factory to build a particular concrete equation set.
    virtual std::string getType() const = 0;

    ///@}

  };
  
}

#endif
