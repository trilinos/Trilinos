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
#include "Panzer_Integrator_CurlBasisDotVector.hpp"
#include "Panzer_ScalarToVector.hpp"
#include "Panzer_Sum.hpp"
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
     *   1/eps dE/dt = 1/mu \nabla \times B -J
     *   dB/dt = - \nabla \times E
     *
     *   in weak form
     */


    valid_parameters.set("Model ID","","Closure model id associated with this equation set");
    valid_parameters.set("Basis Order",1,"Order of the basis");
    valid_parameters.set("Integration Order",2,"Order of the integration");
    valid_parameters.set("Epsilon",8.854187817e-12, "Permittivity of free space");
    valid_parameters.set("Mu",1.2566370614e-6, "Permeability of free space");

    params->validateParametersAndSetDefaults(valid_parameters);
  }

  int basis_order = params->get<int>("Basis Order");
  std::string model_id = params->get<std::string>("Model ID");
  int int_order = params->get<int>("Integration Order");
  eps = params->get<double>("Epsilon");
  mu = params->get<double>("Mu");
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

  this->addClosureModel(model_id);

  this->setupDOFs();
}

// ***********************************************************************
template <typename EvalT>
void mini_em::EquationSet_Maxwell<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
    const panzer::FieldLibrary& fl,
    const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  // ********************
  // EField Equation
  //  eps dE/dt = 1/mu \nabla \times B -J
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
      p.set("Multiplier", eps);
      RCP< PHX::Evaluator<panzer::Traits> > op =
          rcp(new panzer::Integrator_BasisTimesVector<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
      residual_operator_names.push_back(resid);
    }
    {
      std::string resid="RESIDUAL_"+m_Efield_dof_name+"_CURLB_OP";
      ParameterList p("Curl B"+m_Efield_dof_name);
      p.set("Residual Name", resid);
      p.set("Value Name", m_Bfield_dof_name);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", -1.0/mu);
      RCP< PHX::Evaluator<panzer::Traits> > op =
          rcp(new panzer::Integrator_CurlBasisDotVector<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
      residual_operator_names.push_back(resid);
    }
    {
      std::string resid="RESIDUAL_"+m_Efield_dof_name+"_CURRENT_SOURCE";
      ParameterList p("Curl B"+m_Efield_dof_name);
      p.set("Residual Name", resid);
      p.set("Value Name", "CURRENT");
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);
     
      RCP< PHX::Evaluator<panzer::Traits> > op =
          rcp(new panzer::Integrator_BasisTimesVector<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
      residual_operator_names.push_back(resid);
    }

    this->buildAndRegisterResidualSummationEvalautor(fm,m_Efield_dof_name,residual_operator_names);
  }
  {
    //Build the B Field equation
    //dB/dt = - \nabla \times E
    Teuchos::RCP<panzer::IntegrationRule> ir = this->getIntRuleForDOF(m_Bfield_dof_name);
    Teuchos::RCP<panzer::BasisIRLayout> basis = this->getBasisIRLayoutForDOF(m_Bfield_dof_name);

    std::vector<std::string> residual_operator_names;
    {
      std::string resid="RESIDUAL_"+m_Bfield_dof_name+"_TIME_OP";
      ParameterList p("Time Derivative"+m_Bfield_dof_name);
      p.set("Residual Name", resid);
      p.set("Value Name", "DXDT_"+m_Bfield_dof_name);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);
      RCP< PHX::Evaluator<panzer::Traits> > op;
      if ( dimension == 3 )
        op = rcp(new panzer::Integrator_BasisTimesVector<EvalT,panzer::Traits>(p));
      else
        op = rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
      residual_operator_names.push_back(resid);
    }
    {
      std::string resid="RESIDUAL_"+m_Bfield_dof_name+"_CURLE_OP";
      ParameterList p("Curl B"+m_Bfield_dof_name);
      p.set("Residual Name", resid);
      p.set("Value Name", "CURL_"+m_Efield_dof_name);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);
      RCP< PHX::Evaluator<panzer::Traits> > op;
      if ( dimension == 3 )
        op = rcp(new panzer::Integrator_BasisTimesVector<EvalT,panzer::Traits>(p));
      else
        op = rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
      residual_operator_names.push_back(resid);
    }
    this->buildAndRegisterResidualSummationEvalautor(fm,m_Bfield_dof_name,residual_operator_names);
  }



}

// ***********************************************************************



#endif /* SRC_SIMPLEEM_EQUATIONSET_ELECTROMAGNETIC_IMPL_HPP_ */
