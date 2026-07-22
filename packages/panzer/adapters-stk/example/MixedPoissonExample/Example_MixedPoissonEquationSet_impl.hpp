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
#include "Panzer_Integrator_GradBasisDotVector.hpp"
#include "Panzer_Integrator_DivBasisTimesScalar.hpp"
#include "Panzer_ScalarToVector.hpp"
#include "Panzer_Sum.hpp"
#include "Panzer_Constant.hpp"

// ***********************************************************************
template <typename EvalT>
Example::MixedPoissonEquationSet<EvalT>::
MixedPoissonEquationSet(const Teuchos::RCP<Teuchos::ParameterList>& params,
			 const int& default_integration_order,
			 const panzer::CellData& cell_data,
			 const Teuchos::RCP<panzer::GlobalData>& global_data,
			 const bool build_transient_support) :
  panzer::EquationSet_DefaultImpl<EvalT>(params, default_integration_order, cell_data, global_data, build_transient_support )
{
  // ********************
  // Validate and parse parameter list
  // ********************
  {    
    Teuchos::ParameterList valid_parameters;
    this->setDefaultValidParameters(valid_parameters);
    
    valid_parameters.set("Model ID","","Closure model id associated with this equaiton set");
    valid_parameters.set("Integration Order",-1,"Order of the integration rule");
    valid_parameters.set("HGrad Basis Order",-1,"Polynomial order of hgrad basis");
    valid_parameters.set("HDiv Basis Order",-1,"Polynomial order of hdiv basis");
    
    params->validateParametersAndSetDefaults(valid_parameters);
  }
  
  std::string basis_type = "HGrad";
  std::string grad_basis_type = "HDiv";
  int basis_order = params->get<int>("HGrad Basis Order");
  int grad_basis_order = params->get<int>("HDiv Basis Order");
  int integration_order = params->get<int>("Integration Order");
  std::string model_id = params->get<std::string>("Model ID");

  // ********************
  // Panzer uses strings to match fields. In this section we define the
  // name of the fields provided by this equation set. This is a bit strange
  // in that this is not the fields necessarily required by this equation set.
  // For instance for the momentum equations in Navier-Stokes only the velocity
  // fields are added, the pressure field is added by continuity.
  //
  //
  // After all the equation set evaluators are added to a given field manager, the
  // panzer code adds in appropriate scatter evaluators to distribute the
  // entries into the residual and the Jacobian operator. These operators will be
  // "required" by the field manager and will serve as roots of evaluation tree.
  // The leaves of this tree will include the gather evaluators whose job it is to
  // gather the solution from a vector.
  // ********************

  // ********************
  // Assemble DOF names and Residual names
  // ********************

  this->addDOF("PHI",basis_type,basis_order,integration_order);
  this->addDOFGrad("PHI");

  this->addDOF("GRADPHI_FIELD",grad_basis_type,grad_basis_order,integration_order);
  this->addDOFDiv("GRADPHI_FIELD");

   // ********************
   // Build Basis Functions and Integration Rules
   // ********************
   
   this->addClosureModel(model_id);

   this->setupDOFs();
}

// ***********************************************************************
template <typename EvalT>
void Example::MixedPoissonEquationSet<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				      const panzer::FieldLibrary& /* fl */,
				      const Teuchos::ParameterList& /* user_data */) const
{
  using panzer::EvaluatorStyle;
  using panzer::Integrator_DivBasisTimesScalar;
  using panzer::Traits;
  using PHX::Evaluator;
  using std::string;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
   
  Teuchos::RCP<panzer::IntegrationRule> ir = this->getIntRuleForDOF("GRADPHI_FIELD");
  Teuchos::RCP<panzer::BasisIRLayout> basis_v = this->getBasisIRLayoutForDOF("GRADPHI_FIELD"); 
  Teuchos::RCP<panzer::BasisIRLayout> basis_phi = this->getBasisIRLayoutForDOF("PHI"); 

  // This implements the Least-Squares approach for a mixed formulation of
  //
  //  -\nabla\cdot\nabla \phi = f, \phi(D) = 0
  //
  // Specifically the weak form is
  //
  //    (-\nabla\cdot v-f,-\nabla\cdot w) + (\nabla \phi - v,\nabla q) 
  //
  // where w \in HDIV and q \in HGRAD

  // "diffusion" operator (-\nabla\cdot v,-\nabla\cdot w)
  {   
    string resName("RESIDUAL_GRADPHI_FIELD"),
           valName("DIV_GRADPHI_FIELD");
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_DivBasisTimesScalar<EvalT, Traits>(EvaluatorStyle::EVALUATES,
      resName, valName, *basis_v, *ir));
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // Source operator (-f,-\nabla\cdot w)
  {   
    string resName("RESIDUAL_GRADPHI_FIELD"), valName("SOURCE");
    double multiplier(-1);
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_DivBasisTimesScalar<EvalT, Traits>(
      EvaluatorStyle::CONTRIBUTES, resName, valName, *basis_v, *ir,
      multiplier));
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // "diffusion" operator (\nabla \phi,\nabla q)
  {   
    ParameterList p("Source Residual");
    p.set("Residual Name", "RESIDUAL_PHI_DIFFUSION_OP");
    p.set("Flux Name", "GRAD_PHI"); 
    p.set("Basis", basis_phi);
    p.set("IR", ir);
    p.set("Multiplier", 1.0);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_GradBasisDotVector<EvalT,panzer::Traits>(p));
    
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // "diffusion" operator (-v,\nabla q)
  {   
    ParameterList p("Source Residual");
    p.set("Residual Name", "RESIDUAL_PHI_MASS_OP");
    p.set("Flux Name", "GRADPHI_FIELD"); 
    p.set("Basis", basis_phi);
    p.set("IR", ir);
    p.set("Multiplier", -1.0);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_GradBasisDotVector<EvalT,panzer::Traits>(p));
    
    this->template registerEvaluator<EvalT>(fm, op);
  }
  
  // Use a sum operator to form the overall residual for the equation
  {
    std::vector<std::string> sum_names;
    
    // these are the names of the residual values to sum together
    sum_names.push_back("RESIDUAL_PHI_DIFFUSION_OP");
    sum_names.push_back("RESIDUAL_PHI_MASS_OP");

    this->buildAndRegisterResidualSummationEvaluator(fm,"PHI",sum_names);
  }

}

// ***********************************************************************

#endif
