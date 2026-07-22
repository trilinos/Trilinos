// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _MiniEM_EquationSet_Maxwell_impl_hpp_
#define _MiniEM_EquationSet_Maxwell_impl_hpp_


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
#include "Panzer_Integrator_CurlBasisDotVector.hpp"
#include "Panzer_ScalarToVector.hpp"
#include "Panzer_Sum.hpp"
#include "Panzer_ScalarToVector.hpp"
#include "Panzer_Product.hpp"
#include "Panzer_Constant.hpp"

// ***********************************************************************
template <typename EvalT>
mini_em::EquationSet_Maxwell<EvalT>::
EquationSet_Maxwell(const Teuchos::RCP<Teuchos::ParameterList>& params,
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
     *   epsilon dE/dt = 1/mu \nabla \times B - J - sigma E
     *   dB/dt = - \nabla \times E
     *
     *   in weak form
     */


    valid_parameters.set("Model ID","","Closure model id associated with this equation set");
    valid_parameters.set("Basis Order",1,"Order of the basis");
    valid_parameters.set("Integration Order",2,"Order of the integration");
    valid_parameters.set("Permittivity","epsilon","Permittivity");
    valid_parameters.set("Conductivity","sigma","Conductivity");
    valid_parameters.set("Permeability","mu","Permeability");
    valid_parameters.set("Inverse Permeability","1/mu","Inverse Permeability");
    valid_parameters.set("Current","J","Current source");

    params->validateParametersAndSetDefaults(valid_parameters);
  }

  int basis_order = params->get<int>("Basis Order");
  std::string model_id = params->get<std::string>("Model ID");
  int int_order = params->get<int>("Integration Order");
  dimension = cell_data.baseCellDimension();
  TEUCHOS_ASSERT(dimension == 2 || dimension == 3);

  // ********************
  // Setup DOFs and closure models
  // ********************
  {
    m_Efield_dof_name = "E_edge";
    std::string basis_type = "HCurl";
    this->addDOF(m_Efield_dof_name,basis_type,basis_order,int_order);
    this->addDOFCurl(m_Efield_dof_name);
    this->addDOFTimeDerivative(m_Efield_dof_name);
  }

  {
    m_Bfield_dof_name = "B_face";
    std::string basis_type = "HDiv";
    if ( dimension == 2)
      basis_type="Const";
    this->addDOF(m_Bfield_dof_name,basis_type, basis_order,int_order);
    this->addDOFTimeDerivative(m_Bfield_dof_name);
  }

  permittivity_ = params->get<std::string>("Permittivity");
  conductivity_ = params->get<std::string>("Conductivity");
  inverse_permeability_ = params->get<std::string>("Inverse Permeability");
  current_ = params->get<std::string>("Current");
  this->addClosureModel(model_id);

  this->setupDOFs();
}

// ***********************************************************************
template <typename EvalT>
void mini_em::EquationSet_Maxwell<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
    const panzer::FieldLibrary& /* fl */,
    const Teuchos::ParameterList& /* user_data */) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  // ********************
  // EField Equation
  //  eps dE/dt = 1/mu \nabla \times B - J - sigma E
  // ********************
  { 
    Teuchos::RCP<panzer::IntegrationRule> ir = this->getIntRuleForDOF(m_Efield_dof_name);
    Teuchos::RCP<panzer::BasisIRLayout> basis = this->getBasisIRLayoutForDOF(m_Efield_dof_name);

    std::vector<std::string> residual_operator_names;
    {
      std::string resid="RESIDUAL_"+m_Efield_dof_name+"_TIME_OP";
      ParameterList p("Time Derivative"+m_Efield_dof_name);
      p.set("Residual Name", resid);
      p.set("Value Name", "DXDT_"+m_Efield_dof_name);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);
      const std::vector<std::string> fieldMultiplier = {permittivity_};
      p.set("Field Multipliers",Teuchos::rcpFromRef(fieldMultiplier));

      RCP< PHX::Evaluator<panzer::Traits> > op =
          rcp(new panzer::Integrator_BasisTimesVector<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
      residual_operator_names.push_back(resid);
    }
    {
      std::string resid="RESIDUAL_"+m_Efield_dof_name+"_CONDUCTIVITY";
      ParameterList p(m_Efield_dof_name);
      p.set("Residual Name", resid);
      p.set("Value Name", m_Efield_dof_name);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Tensor Name", conductivity_);
      RCP<PHX::Evaluator<panzer::Traits> > op =
        rcp(new panzer::Integrator_BasisTimesTensorTimesVector<EvalT, panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
      residual_operator_names.push_back(resid);
    }
    {
      using panzer::EvaluatorStyle;
      using panzer::Integrator_CurlBasisDotVector;
      using panzer::Traits;
      using PHX::Evaluator;
      using std::string;
      string resName("RESIDUAL_" + m_Efield_dof_name),
             valName(m_Bfield_dof_name);
      double multiplier(-1.0);
      std::vector<std::string> fieldMultiplier = {inverse_permeability_};
      RCP<Evaluator<Traits>> op = rcp(new
        Integrator_CurlBasisDotVector<EvalT, Traits>(
        EvaluatorStyle::CONTRIBUTES, resName, valName, *basis, *ir,
        multiplier,fieldMultiplier));
      this->template registerEvaluator<EvalT>(fm, op);
    }
    {
      std::string resid="RESIDUAL_"+m_Efield_dof_name+"_CURRENT_SOURCE";
      ParameterList p("Curl B"+m_Efield_dof_name);
      p.set("Residual Name", resid);
      p.set("Value Name", current_);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);
     
      RCP< PHX::Evaluator<panzer::Traits> > op =
          rcp(new panzer::Integrator_BasisTimesVector<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
      residual_operator_names.push_back(resid);
    }

    this->buildAndRegisterResidualSummationEvaluator(fm,m_Efield_dof_name,residual_operator_names);
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
    //Build the B Field equation
    //dB/dt = - \nabla \times E
    RCP<IntegrationRule> ir = this->getIntRuleForDOF(m_Bfield_dof_name);
    RCP<BasisIRLayout> basis = this->getBasisIRLayoutForDOF(m_Bfield_dof_name);

    vector<string> residual_operator_names;
    {
      string valName("DXDT_" + m_Bfield_dof_name);
      double multiplier(1);
      RCP<Evaluator<Traits>> op;
      if (dimension == 3)
      {
        string resName("RESIDUAL_" + m_Bfield_dof_name + "_TIME_OP");
        ParameterList p("Time Derivative" + m_Bfield_dof_name);
        p.set("Residual Name", resName);
        p.set("Value Name", valName);
        p.set("Basis", basis);
        p.set("IR", ir);
        p.set("Multiplier", multiplier);
        op = rcp(new Integrator_BasisTimesVector<EvalT, Traits>(p));
        residual_operator_names.push_back(resName);
      }
      else
      {
        string resName("RESIDUAL_" + m_Bfield_dof_name);
        op = rcp(new Integrator_BasisTimesScalar<EvalT, Traits>(
          EvaluatorStyle::EVALUATES, resName, valName, *basis, *ir,
          multiplier));
      }
      this->template registerEvaluator<EvalT>(fm, op);
    }
    {
      string valName("CURL_" + m_Efield_dof_name);
      double multiplier(1);
      RCP<Evaluator<Traits>> op;
      if (dimension == 3)
      {
        string resName("RESIDUAL_" + m_Bfield_dof_name + "_CURLE_OP");
        ParameterList p("Curl B" + m_Bfield_dof_name);
        p.set("Residual Name", resName);
        p.set("Value Name", "CURL_" + m_Efield_dof_name);
        p.set("Basis", basis);
        p.set("IR", ir);
        p.set("Multiplier", multiplier);
        op = rcp(new Integrator_BasisTimesVector<EvalT, Traits>(p));
        residual_operator_names.push_back(resName);
      }
      else
      {
        string resName("RESIDUAL_" + m_Bfield_dof_name);
        op = rcp(new Integrator_BasisTimesScalar<EvalT, Traits>(
          EvaluatorStyle::CONTRIBUTES, resName, valName, *basis, *ir,
          multiplier));
      }
      this->template registerEvaluator<EvalT>(fm, op);
    }
    if (dimension == 3)
      this->buildAndRegisterResidualSummationEvaluator(fm, m_Bfield_dof_name, residual_operator_names);
  }



}

// ***********************************************************************



#endif /* SRC_SIMPLEEM_EQUATIONSET_ELECTROMAGNETIC_IMPL_HPP_ */
