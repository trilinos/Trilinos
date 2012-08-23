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
EquationSet_Energy(const panzer::InputEquationSet& ies,
		   const panzer::CellData& cell_data,
		   const Teuchos::RCP<panzer::GlobalData>& global_data,
		   const bool build_transient_support) :
  panzer::EquationSet_DefaultImpl<EvalT>(ies, cell_data, global_data, build_transient_support )
{
  this->m_eqset_prefix = ies.prefix;

  // ********************
  // Assemble DOF names and Residual names
  // ********************
  this->m_dof_names->push_back(this->m_eqset_prefix+"TEMPERATURE");

  this->m_dof_gradient_names->push_back("GRAD_"+this->m_eqset_prefix+"TEMPERATURE");

  if (this->m_build_transient_support)
    this->m_dof_time_derivative_names->push_back("DOT_"+this->m_eqset_prefix+"TEMPERATURE");

  this->m_residual_names->push_back("RESIDUAL_"+this->m_eqset_prefix+"TEMPERATURE");

  this->m_scatter_name = "Scatter_RESIDUAL_"+this->m_eqset_prefix+"TEMPERATURE";

  // ********************
  // Build Basis Functions and Integration Rules
  // ********************
  
  this->setupDOFs(cell_data.baseCellDimension());

  // ********************
  // Parse valid options
  // ********************
  {    
    Teuchos::ParameterList valid_parameters;
    Teuchos::setStringToIntegralParameter<int>(
      "CONVECTION",
      "OFF",
      "Enables or disables convection term in the energy equation",
      Teuchos::tuple<std::string>("ON","OFF"),
      &valid_parameters
      );

    // Don't corrupt original input, create a copy
    Teuchos::ParameterList tmp_ies_params = this->getInputEquationSet().params;
    tmp_ies_params.setName(ies.prefix+"ENERGY");
    tmp_ies_params.validateParametersAndSetDefaults(valid_parameters);
    
    m_do_convection = tmp_ies_params.get<std::string>("CONVECTION");
  }

}

// ***********************************************************************
template <typename EvalT>
void user_app::EquationSet_Energy<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				      const std::vector<std::pair<std::string,Teuchos::RCP<panzer::BasisIRLayout> > > & dofs,
				      const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  
  const std::string& prefix = this->m_eqset_prefix;

  // ********************
  // Energy Equation
  // ********************

  // Transient Operator
  if (this->m_build_transient_support) {
    ParameterList p("Transient Residual");
    p.set("Residual Name", "RESIDUAL_"+prefix+"TEMPERATURE_TRANSIENT_OP");
    p.set("Value Name", "DOT_"+prefix+"TEMPERATURE");
    p.set("Basis", this->m_basis);
    p.set("IR", this->m_int_rule);
    p.set("Multiplier", 1.0);
    Teuchos::RCP<std::vector<std::string> > fms = 
      Teuchos::rcp(new std::vector<std::string>);
    fms->push_back(prefix+"DENSITY");
    fms->push_back(prefix+"HEAT_CAPACITY");
    p.set< Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers",fms);

    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_TransientBasisTimesScalar<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
  }

  // Diffusion Operator
  {
    double thermal_conductivity = 1.0;

    ParameterList p("Diffusion Residual");
    p.set("Residual Name", "RESIDUAL_"+prefix+"TEMPERATURE_DIFFUSION_OP");
    p.set("Flux Name", "GRAD_"+prefix+"TEMPERATURE");
    p.set("Basis", this->m_basis);
    p.set("IR", this->m_int_rule);
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
      scalar_names->push_back(prefix+"UX");
      scalar_names->push_back(prefix+"UY");
      p.set<RCP<const std::vector<std::string> > >("Scalar Names", scalar_names);
      p.set("Vector Name", prefix+"U");
      p.set("Data Layout Scalar",this->m_int_rule->dl_scalar);
      p.set("Data Layout Vector",this->m_int_rule->dl_vector);

      RCP< PHX::Evaluator<panzer::Traits> > op = 
	rcp(new panzer::ScalarToVector<EvalT,panzer::Traits>(p));
      
      fm.template registerEvaluator<EvalT>(op);
    }

    // Evaluator to assemble convection term
    {
      ParameterList p("Convection Operator");
      p.set("IR", this->m_int_rule);
      p.set("Operator Name", prefix+"TEMPERATURE_CONVECTION_OP");
      p.set("A Name", prefix+"U");
      p.set("Gradient Name", "GRAD_"+prefix+"TEMPERATURE");
      p.set("Multiplier", 1.0);

      RCP< PHX::Evaluator<panzer::Traits> > op = 
	rcp(new user_app::Convection<EvalT,panzer::Traits>(p));
      
      fm.template registerEvaluator<EvalT>(op);
    }

    // Integration operator (could sum this into source for efficiency)
    {
      ParameterList p("Convection Residual");
      p.set("Residual Name","RESIDUAL_"+prefix+"TEMPERATURE_CONVECTION_OP");
      p.set("Value Name", prefix+"TEMPERATURE_CONVECTION_OP");
      p.set("Basis", this->m_basis);
      p.set("IR", this->m_int_rule);
      p.set("Multiplier", 1.0);
      Teuchos::RCP<std::vector<std::string> > fms = 
	Teuchos::rcp(new std::vector<std::string>);
      fms->push_back(prefix+"DENSITY");
      fms->push_back(prefix+"HEAT_CAPACITY");
      p.set< Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers",fms);
      
      RCP< PHX::Evaluator<panzer::Traits> > op = 
	rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));
      
      fm.template registerEvaluator<EvalT>(op);
    }
  }

  // Source Operator
  {   
    ParameterList p("Source Residual");
    p.set("Residual Name", "RESIDUAL_"+prefix+"TEMPERATURE_SOURCE_OP");
    p.set("Value Name", "SOURCE_"+prefix+"TEMPERATURE");
    p.set("Basis", this->m_basis);
    p.set("IR", this->m_int_rule);
    p.set("Multiplier", -1.0);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
  }

  // Use a sum operator to form the overall residual for the equation
  // - this way we avoid loading each operator separately into the
  // global residual and Jacobian
  {
    ParameterList p;
    p.set("Sum Name", "RESIDUAL_"+prefix+"TEMPERATURE");

    RCP<std::vector<std::string> > sum_names = 
      rcp(new std::vector<std::string>);

    sum_names->push_back("RESIDUAL_"+prefix+"TEMPERATURE_DIFFUSION_OP");
    sum_names->push_back("RESIDUAL_"+prefix+"TEMPERATURE_SOURCE_OP");
    if (m_do_convection == "ON")
      sum_names->push_back("RESIDUAL_"+prefix+"TEMPERATURE_CONVECTION_OP");
    if (this->m_build_transient_support)
      sum_names->push_back("RESIDUAL_"+prefix+"TEMPERATURE_TRANSIENT_OP");

    p.set("Values Names", sum_names);
    p.set("Data Layout", this->m_basis->functional);

    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

}

// ***********************************************************************

#endif
