// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _MiniEM_EquationSet_Darcy_impl_hpp_
#define _MiniEM_EquationSet_Darcy_impl_hpp_


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
#include "Panzer_Integrator_BasisTimesVector.hpp"
#include "Panzer_Integrator_BasisTimesTensorTimesVector.hpp"
#include "Panzer_Integrator_DivBasisTimesScalar.hpp"
#include "Panzer_ScalarToVector.hpp"
#include "Panzer_Sum.hpp"
#include "Panzer_ScalarToVector.hpp"
#include "Panzer_Product.hpp"
#include "Panzer_Constant.hpp"

// ***********************************************************************
template <typename EvalT>
mini_em::EquationSet_Darcy<EvalT>::
EquationSet_Darcy(const Teuchos::RCP<Teuchos::ParameterList>& params,
    const int& default_integration_order,
    const panzer::CellData& cell_data,
    const Teuchos::RCP<panzer::GlobalData>& global_data,
    const bool build_transient_support) :
    panzer::EquationSet_DefaultImpl<EvalT>(params,default_integration_order,cell_data,global_data,build_transient_support ) {
  // ********************
  // Validate and parse parameter list
  // ********************
  {
    Teuchos::ParameterList valid_parameters;
    this->setDefaultValidParameters(valid_parameters);
    /*  Equations are
     *   du/dt     - \nabla \dot sigma = f
     *  - \nabla u + 1/\kappa sigma    = 0
     *
     *   in weak form
     */


    valid_parameters.set("Model ID","","Closure model id associated with this equation set");
    valid_parameters.set("Basis Order",1,"Order of the basis");
    valid_parameters.set("Integration Order",2,"Order of the integration");
    valid_parameters.set("Diffusivity","kappa","Diffusivity");
    valid_parameters.set("Inverse Diffusivity","1/kappa","Inverse Diffusivity");
    valid_parameters.set("Forcing","forcing","Forcing");

    params->validateParametersAndSetDefaults(valid_parameters);
  }

  int basis_order = params->get<int>("Basis Order");
  std::string model_id = params->get<std::string>("Model ID");
  int int_order = params->get<int>("Integration Order");
  dimension = cell_data.baseCellDimension();

  // ********************
  // Setup DOFs and closure models
  // ********************
  {
    m_u_field_dof_name = "u";
    std::string basis_type = "HVol";
    this->addDOF(m_u_field_dof_name,basis_type,basis_order-1,int_order);
    this->addDOFTimeDerivative(m_u_field_dof_name);
  }

  {
    m_sigma_field_dof_name = "sigma";
    std::string basis_type = "HDiv";
    this->addDOF(m_sigma_field_dof_name,basis_type,basis_order,int_order);
    this->addDOFDiv(m_sigma_field_dof_name);
  }

  inverse_diffusivity_ = params->get<std::string>("Inverse Diffusivity");
  forcing_ = params->get<std::string>("Forcing");
  this->addClosureModel(model_id);

  this->setupDOFs();
}

// ***********************************************************************
template <typename EvalT>
void mini_em::EquationSet_Darcy<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
    const panzer::FieldLibrary& /* fl */,
    const Teuchos::ParameterList& /* user_data */) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  // ********************
  //   du/dt - \nabla \dot sigma = f
  // ********************
  {
    Teuchos::RCP<panzer::IntegrationRule> ir = this->getIntRuleForDOF(m_u_field_dof_name);
    Teuchos::RCP<panzer::BasisIRLayout> basis = this->getBasisIRLayoutForDOF(m_u_field_dof_name);

    std::vector<std::string> residual_operator_names;
    {
      std::string resid="RESIDUAL_"+m_u_field_dof_name+"_TIME_OP";
      ParameterList p("Time Derivative"+m_u_field_dof_name);
      p.set("Residual Name", resid);
      p.set("Value Name", "DXDT_"+m_u_field_dof_name);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);

      RCP< PHX::Evaluator<panzer::Traits> > op =
          rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
      residual_operator_names.push_back(resid);
    }
    {
      std::string resName("RESIDUAL_" + m_u_field_dof_name + "_DIV_SIGMA_OP");
      ParameterList p("Div sigma" + m_u_field_dof_name);
      p.set("Residual Name", resName);
      p.set("Value Name", "DIV_" + m_sigma_field_dof_name);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", -1.0);
      RCP< PHX::Evaluator<panzer::Traits> > op = rcp(new panzer::Integrator_BasisTimesScalar<EvalT, panzer::Traits>(p));
      this->template registerEvaluator<EvalT>(fm, op);
      residual_operator_names.push_back(resName);
    }
    {
      std::string resid="RESIDUAL_"+m_u_field_dof_name+"_FORCING";
      ParameterList p("Forcing"+m_u_field_dof_name);
      p.set("Residual Name", resid);
      p.set("Value Name", forcing_);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", -1.0);

      RCP< PHX::Evaluator<panzer::Traits> > op =
          rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
      residual_operator_names.push_back(resid);
    }

    this->buildAndRegisterResidualSummationEvaluator(fm,m_u_field_dof_name,residual_operator_names);
  }
  {
    using panzer::BasisIRLayout;
    using panzer::EvaluatorStyle;
    using panzer::IntegrationRule;
    using panzer::Integrator_BasisTimesScalar;
    using panzer::Integrator_BasisTimesVector;
    using panzer::Traits;
    using PHX::Evaluator;
    using std::string;
    using std::vector;
    using Teuchos::RCP;
    // ********************
    // - \nabla u + 1/\kappa sigma = 0
    // ********************
    RCP<IntegrationRule> ir = this->getIntRuleForDOF(m_sigma_field_dof_name);
    RCP<BasisIRLayout> basis = this->getBasisIRLayoutForDOF(m_sigma_field_dof_name);

    vector<string> residual_operator_names;
    {
      RCP<Evaluator<Traits>> op;
      string resName("RESIDUAL_" + m_sigma_field_dof_name + "_GRAD_U_OP");
      ParameterList p("Grad u" + m_sigma_field_dof_name);
      p.set("Residual Name", resName);
      p.set("Value Name", m_u_field_dof_name);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);
      op = rcp(new panzer::Integrator_DivBasisTimesScalar<EvalT, Traits>(p));
      residual_operator_names.push_back(resName);

      this->template registerEvaluator<EvalT>(fm, op);
    }
    {
      RCP<Evaluator<Traits> > op;
      string resName("RESIDUAL_" + m_sigma_field_dof_name + "_MASS_OP");
      ParameterList p("Mass" + m_sigma_field_dof_name);
      p.set("Residual Name", resName);
      p.set("Value Name", m_sigma_field_dof_name);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);
      const std::vector<std::string> fieldMultiplier = {inverse_diffusivity_};
      p.set("Field Multipliers",Teuchos::rcpFromRef(fieldMultiplier));
      op = rcp(new Integrator_BasisTimesVector<EvalT, Traits>(p));
      residual_operator_names.push_back(resName);

      this->template registerEvaluator<EvalT>(fm, op);
    }

    this->buildAndRegisterResidualSummationEvaluator(fm, m_sigma_field_dof_name, residual_operator_names);
  }



}

// ***********************************************************************



#endif /* _MiniEM_EquationSet_Darcy_impl_hpp_ */
