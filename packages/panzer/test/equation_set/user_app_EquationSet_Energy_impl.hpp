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

#ifndef USER_APP_EQUATIONSET_ENERGY_T_HPP
#define USER_APP_EQUATIONSET_ENERGY_T_HPP

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
#include "Panzer_Integrator_TransientBasisTimesScalar.hpp"
#include "Panzer_Integrator_GradBasisDotVector.hpp"
#include "Panzer_ScalarToVector.hpp"
#include "Panzer_Sum.hpp"
#include "Panzer_Constant.hpp"
#include "user_app_Convection.hpp"

// ***********************************************************************
template <typename EvalT>
user_app::EquationSet_Energy<EvalT>::
EquationSet_Energy(const Teuchos::RCP<Teuchos::ParameterList>& params,
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

    Teuchos::setStringToIntegralParameter<int>(
      "CONVECTION",
      "OFF",
      "Enables or disables convection term in the energy equation",
      Teuchos::tuple<std::string>("ON","OFF"),
      &valid_parameters
      );    

    params->validateParametersAndSetDefaults(valid_parameters);
  }

  m_do_convection = params->get<std::string>("CONVECTION");
  m_prefix = params->get<std::string>("Prefix");
  std::string basis_type = params->get<std::string>("Basis Type");
  int basis_order = params->get<int>("Basis Order");
  std::string model_id = params->get<std::string>("Model ID");
  int integration_order = params->get<int>("Integration Order");

  // ********************
  // Setup DOFs and closure models
  // ********************
  {
    m_dof_name = m_prefix+"TEMPERATURE";

    this->addDOF(m_dof_name,basis_type,basis_order,integration_order);
    this->addDOFGrad(m_dof_name);
    if (this->buildTransientSupport())
      this->addDOFTimeDerivative(m_dof_name);
  }

  this->addClosureModel(model_id);

  this->setupDOFs();
}

// ***********************************************************************
template <typename EvalT>
void user_app::EquationSet_Energy<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				      const panzer::FieldLibrary& fl,
				      const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  // ********************
  // Energy Equation
  // ********************

  Teuchos::RCP<const panzer::PureBasis> pb = fl.lookupBasis(m_dof_name);
  Teuchos::RCP<panzer::IntegrationRule> ir = this->getIntRuleForDOF(m_dof_name);
  Teuchos::RCP<panzer::BasisIRLayout> basis = panzer::basisIRLayout(pb,*ir); 

  // Transient Operator
  if (this->buildTransientSupport()) {
    ParameterList p("Transient Residual");
    p.set("Residual Name", "RESIDUAL_"+m_prefix+"TEMPERATURE_TRANSIENT_OP");
    p.set("Value Name", "DXDT_"+m_prefix+"TEMPERATURE");
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", 1.0);
    Teuchos::RCP<std::vector<std::string> > fms = 
      Teuchos::rcp(new std::vector<std::string>);
    fms->push_back(m_prefix+"DENSITY");
    fms->push_back(m_prefix+"HEAT_CAPACITY");
    p.set< Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers",fms);

    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_TransientBasisTimesScalar<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
  }

  // Diffusion Operator
  {
    double thermal_conductivity = 1.0;

    ParameterList p("Diffusion Residual");
    p.set("Residual Name", "RESIDUAL_"+m_prefix+"TEMPERATURE_DIFFUSION_OP");
    p.set("Flux Name", "GRAD_"+m_prefix+"TEMPERATURE");
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", thermal_conductivity);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_GradBasisDotVector<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }
  
  // Convection Operator
  if (m_do_convection == "ON") {

    // Combine scalar velocities into a velocity vector
    {
      ParameterList p("Velocity: ScalarToVector");
      RCP<std::vector<std::string> > scalar_names = rcp(new std::vector<std::string>);
      scalar_names->push_back(m_prefix+"UX");
      scalar_names->push_back(m_prefix+"UY");
      p.set<RCP<const std::vector<std::string> > >("Scalar Names", scalar_names);
      p.set("Vector Name", m_prefix+"U");
      p.set("Data Layout Scalar",ir->dl_scalar);
      p.set("Data Layout Vector",ir->dl_vector);

      RCP< PHX::Evaluator<panzer::Traits> > op = 
	rcp(new panzer::ScalarToVector<EvalT,panzer::Traits>(p));
      
      fm.template registerEvaluator<EvalT>(op);
    }

    // Evaluator to assemble convection term
    {
      ParameterList p("Convection Operator");
      p.set("IR", ir);
      p.set("Operator Name", m_prefix+"TEMPERATURE_CONVECTION_OP");
      p.set("A Name", m_prefix+"U");
      p.set("Gradient Name", "GRAD_"+m_prefix+"TEMPERATURE");
      p.set("Multiplier", 1.0);

      RCP< PHX::Evaluator<panzer::Traits> > op = 
	rcp(new user_app::Convection<EvalT,panzer::Traits>(p));
      
      fm.template registerEvaluator<EvalT>(op);
    }

    // Integration operator (could sum this into source for efficiency)
    {
      ParameterList p("Convection Residual");
      p.set("Residual Name","RESIDUAL_"+m_prefix+"TEMPERATURE_CONVECTION_OP");
      p.set("Value Name", m_prefix+"TEMPERATURE_CONVECTION_OP");
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);
      Teuchos::RCP<std::vector<std::string> > fms = 
	Teuchos::rcp(new std::vector<std::string>);
      fms->push_back(m_prefix+"DENSITY");
      fms->push_back(m_prefix+"HEAT_CAPACITY");
      p.set< Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers",fms);
      
      RCP< PHX::Evaluator<panzer::Traits> > op = 
	rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));
      
      fm.template registerEvaluator<EvalT>(op);
    }
  }

  // Source Operator
  {   
    ParameterList p("Source Residual");
    p.set("Residual Name", "RESIDUAL_"+m_prefix+"TEMPERATURE_SOURCE_OP");
    p.set("Value Name", "SOURCE_"+m_prefix+"TEMPERATURE");
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", -1.0);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
  }

  // Use a sum operator to form the overall residual for the equation
  // - this way we avoid loading each operator separately into the
  // global residual and Jacobian
  {
    std::vector<std::string> residual_operator_names;

    residual_operator_names.push_back("RESIDUAL_"+m_prefix+"TEMPERATURE_DIFFUSION_OP");
    residual_operator_names.push_back("RESIDUAL_"+m_prefix+"TEMPERATURE_SOURCE_OP");
    if (m_do_convection == "ON")
      residual_operator_names.push_back("RESIDUAL_"+m_prefix+"TEMPERATURE_CONVECTION_OP");
    if (this->buildTransientSupport())
      residual_operator_names.push_back("RESIDUAL_"+m_prefix+"TEMPERATURE_TRANSIENT_OP");

    this->buildAndRegisterResidualSummationEvalautor(fm,m_dof_name,residual_operator_names);
  }

}

// ***********************************************************************

#endif
