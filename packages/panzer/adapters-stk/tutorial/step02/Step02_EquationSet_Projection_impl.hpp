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

#ifndef __Step02_EquationSet_Projection_impl_hpp__
#define __Step02_EquationSet_Projection_impl_hpp__

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

// include evaluators here
#include "Panzer_Integrator_BasisTimesScalar.hpp"

// ***********************************************************************
template <typename EvalT>
user_app::EquationSet_Projection<EvalT>::
EquationSet_Projection(const Teuchos::RCP<Teuchos::ParameterList>& params,
		   const int& default_integration_order,
		   const panzer::CellData& cell_data,
		   const Teuchos::RCP<panzer::GlobalData>& global_data,
		   const bool build_transient_support) :
  panzer::EquationSet_DefaultImpl<EvalT>(params,default_integration_order,cell_data,global_data,build_transient_support )
{
  // ********************
  // Validate and parse parameter list
  // ********************
  {    
    Teuchos::ParameterList valid_parameters;
    this->setDefaultValidParameters(valid_parameters);

    valid_parameters.set("Model ID","","Closure model id associated with this equaiton set");
    valid_parameters.set("Prefix","","Prefix for using multiple instatiations of thei equation set");
    valid_parameters.set("Basis Type","HGrad","Type of Basis to use");
    valid_parameters.set("Basis Order",1,"Order of the basis");
    valid_parameters.set("Integration Order",-1,"Order of the integration rule");

    params->validateParametersAndSetDefaults(valid_parameters);
  }

  std::string model_id   = params->get<std::string>("Model ID");
  std::string prefix     = params->get<std::string>("Prefix");
  std::string basis_type = params->get<std::string>("Basis Type");
  int basis_order        = params->get<int>("Basis Order");
  int integration_order  = params->get<int>("Integration Order");

  // ********************
  // Setup DOFs and closure models
  // ********************
  {
    dof_name_ = prefix+"U";

    this->addDOF(dof_name_,basis_type,basis_order,integration_order);
  }

  this->addClosureModel(model_id);

  this->setupDOFs();
}

// ***********************************************************************
template <typename EvalT>
void user_app::EquationSet_Projection<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				      const panzer::FieldLibrary& fl,
				      const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // define some special strings to use
  const std::string residual_projection_term     = "RESIDUAL_"+dof_name_+"_PROJECTION";
  const std::string residual_projection_src_term = "RESIDUAL_"+dof_name_+"_PROJECTION_SOURCE";
  
  const std::string projection_src_name = dof_name_+"_SOURCE"; 
    // this must be satisfied by the closure model

  // ********************
  // Projection Equation
  // ********************

  RCP<panzer::IntegrationRule> ir  = this->getIntRuleForDOF(dof_name_); 
  RCP<panzer::BasisIRLayout> basis = this->getBasisIRLayoutForDOF(dof_name_); 

  // Projection operator (U,phi)
  {
    ParameterList p;
    p.set("Residual Name", residual_projection_term);
    p.set("Value Name",    dof_name_);
    p.set("Basis",         basis);
    p.set("IR",            ir);
    p.set("Multiplier",    1.0);

    RCP<PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));
    
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // Source operator -(u_source,phi)
  {
    ParameterList p;
    p.set("Residual Name", residual_projection_src_term);
    p.set("Value Name",    projection_src_name);
    p.set("Basis",         basis);
    p.set("IR",            ir);
    p.set("Multiplier",    -1.0);

    RCP<PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));
    
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // Use a sum operator to form the overall residual for the equation
  // - this way we avoid loading each operator separately into the
  // global residual and Jacobian
  {
    std::vector<std::string> residual_operator_names;

    residual_operator_names.push_back(residual_projection_term);
    residual_operator_names.push_back(residual_projection_src_term);

    // build a sum evaluator
    this->buildAndRegisterResidualSummationEvaluator(fm,dof_name_,residual_operator_names);
  }

}

// ***********************************************************************

#endif
