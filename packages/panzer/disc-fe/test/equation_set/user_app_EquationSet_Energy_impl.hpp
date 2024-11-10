// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    valid_parameters.set("CONVECTION", "OFF",
      "Enables or disables convection term in the energy equation",
      rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("ON", "OFF"))));

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
                                      const panzer::FieldLibrary& /* fl */,
                                      const Teuchos::ParameterList& /* user_data */) const
{
  using panzer::BasisIRLayout;
  using panzer::EvaluatorStyle;
  using panzer::IntegrationRule;
  using panzer::Integrator_BasisTimesScalar;
  using panzer::Integrator_GradBasisDotVector;
  using panzer::ScalarToVector;
  using panzer::Traits;
  using PHX::Evaluator;
  using std::string;
  using std::vector;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using user_app::Convection;
  // ********************
  // Energy Equation
  // ********************

  RCP<IntegrationRule> ir = this->getIntRuleForDOF(m_dof_name); 
  RCP<BasisIRLayout> basis = this->getBasisIRLayoutForDOF(m_dof_name); 

  // Transient Operator
  if (this->buildTransientSupport())
  {
    string resName("RESIDUAL_" + m_dof_name),
           valName("DXDT_" + m_prefix + "TEMPERATURE");
    double multiplier(1);
    vector<string> fieldMultipliers{m_prefix + "DENSITY",
      m_prefix + "HEAT_CAPACITY"};
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
      resName, valName, *basis, *ir, multiplier, fieldMultipliers));
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // Diffusion Operator
  {
    double thermal_conductivity = 1.0;

    ParameterList p("Diffusion Residual");
    p.set("Residual Name", "RESIDUAL_"+m_dof_name);
    p.set("Flux Name", "GRAD_"+m_prefix+"TEMPERATURE");
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", thermal_conductivity);
    
    RCP< Evaluator<Traits> > op = 
      rcp(new Integrator_GradBasisDotVector<EvalT,Traits>(p));

    this->template registerEvaluator<EvalT>(fm, op);
  }
  
  // Convection Operator
  if (m_do_convection == "ON") {

    // Combine scalar velocities into a velocity vector
    {
      ParameterList p("Velocity: ScalarToVector");
      RCP<vector<string> > scalar_names = rcp(new vector<string>);
      scalar_names->push_back(m_prefix+"UX");
      scalar_names->push_back(m_prefix+"UY");
      p.set<RCP<const vector<string> > >("Scalar Names", scalar_names);
      p.set("Vector Name", m_prefix+"U");
      p.set("Data Layout Scalar",ir->dl_scalar);
      p.set("Data Layout Vector",ir->dl_vector);

      RCP< Evaluator<Traits> > op = 
        rcp(new ScalarToVector<EvalT,Traits>(p));
      
      this->template registerEvaluator<EvalT>(fm, op);
    }

    // Evaluator to assemble convection term
    {
      ParameterList p("Convection Operator");
      p.set("IR", ir);
      p.set("Operator Name", m_prefix+"TEMPERATURE_CONVECTION_OP");
      p.set("A Name", m_prefix+"U");
      p.set("Gradient Name", "GRAD_"+m_prefix+"TEMPERATURE");
      p.set("Multiplier", 1.0);

      RCP< Evaluator<Traits> > op = 
        rcp(new Convection<EvalT,Traits>(p));
      
      this->template registerEvaluator<EvalT>(fm, op);
    }

    // Integration operator (could sum this into source for efficiency)
    {
      string resName("RESIDUAL_" + m_dof_name),
             valName(m_prefix + "TEMPERATURE_CONVECTION_OP");
      double multiplier(1);
      vector<string> fieldMultipliers{m_prefix + "DENSITY",
        m_prefix + "HEAT_CAPACITY"};
      RCP<Evaluator<Traits>> op = rcp(new
        Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
        resName, valName, *basis, *ir, multiplier, fieldMultipliers));
      this->template registerEvaluator<EvalT>(fm, op);
    }
  }

  // Source Operator
  {   
    string resName("RESIDUAL_" + m_dof_name),
           valName("SOURCE_" + m_prefix + "TEMPERATURE");
    double multiplier(-1);
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
      resName, valName, *basis, *ir, multiplier));
    this->template registerEvaluator<EvalT>(fm, op);
  }
}

// ***********************************************************************

#endif
