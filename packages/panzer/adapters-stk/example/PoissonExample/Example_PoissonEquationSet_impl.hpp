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
#include "Panzer_Integrator_GradBasisDotVector.hpp"
#include "Panzer_ScalarToVector.hpp"
#include "Panzer_Sum.hpp"
#include "Panzer_Constant.hpp"

// ***********************************************************************
template <typename EvalT>
Example::PoissonEquationSet<EvalT>::
PoissonEquationSet(const Teuchos::RCP<Teuchos::ParameterList>& params,
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
    valid_parameters.set("Basis Type","HGrad","Type of Basis to use");
    valid_parameters.set("Basis Order",1,"Order of the basis");
    valid_parameters.set("Integration Order",-1,"Order of the integration rule");
    
    params->validateParametersAndSetDefaults(valid_parameters);
  }
  
  std::string basis_type = params->get<std::string>("Basis Type");
  int basis_order = params->get<int>("Basis Order");
  int integration_order = params->get<int>("Integration Order");
  std::string model_id = params->get<std::string>("Model ID");

   // ********************
   // Panzer uses strings to match fields. In this section we define the
   // name of the fields provided by this equation set. This is a bit strange
   // in that this is not the fields necessarily required by this equation set.
   // For instance for the momentum equations in Navier-Stokes only the velocity
   // fields are added, the pressure field is added by continuity.
   //
   // In this case "TEMPERATURE" is the lone field.  We also name the gradient
   // for this field. These names automatically generate evaluators for "TEMPERATURE"
   // and "GRAD_TEMPERATURE" gathering the basis coefficients of "TEMPERATURE" and
   // the values of the TEMPERATURE and GRAD_TEMPERATURE fields at quadrature points.
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

   this->addDOF("TEMPERATURE",basis_type,basis_order,integration_order);
   this->addDOFGrad("TEMPERATURE");
   if (this->buildTransientSupport())
     this->addDOFTimeDerivative("TEMPERATURE");

   // ********************
   // Build Basis Functions and Integration Rules
   // ********************
   
   this->addClosureModel(model_id);

   this->setupDOFs();
}

// ***********************************************************************
template <typename EvalT>
void Example::PoissonEquationSet<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				      const panzer::FieldLibrary& /* fl */,
				      const Teuchos::ParameterList& /* user_data */) const
{
  using panzer::BasisIRLayout;
  using panzer::EvaluatorStyle;
  using panzer::IntegrationRule;
  using panzer::Integrator_BasisTimesScalar;
  using panzer::Integrator_GradBasisDotVector;
  using panzer::Traits;
  using PHX::Evaluator;
  using std::string;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  
  RCP<IntegrationRule> ir = this->getIntRuleForDOF("TEMPERATURE");
  RCP<BasisIRLayout> basis = this->getBasisIRLayoutForDOF("TEMPERATURE"); 

  // ********************
  // Energy Equation
  // ********************

  // Transient Operator: Assembles \int \dot{T} v
  if (this->buildTransientSupport())
  {
    string resName("RESIDUAL_TEMPERATURE"), valName("DXDT_TEMPERATURE");
    double multiplier(1);
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
      resName, valName, *basis, *ir, multiplier));
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // Diffusion Operator: Assembles \int \nabla T \cdot \nabla v
  {
    double thermal_conductivity = 1.0;

    ParameterList p("Diffusion Residual");
    p.set("Residual Name", "RESIDUAL_TEMPERATURE");
    p.set("Flux Name", "GRAD_TEMPERATURE"); // this field is constructed by the panzer library
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", thermal_conductivity);
    
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_GradBasisDotVector<EvalT, Traits>(p));

    this->template registerEvaluator<EvalT>(fm, op);
  }
  
  // Source Operator
  {   
    string resName("RESIDUAL_TEMPERATURE"),
           valName("SOURCE_TEMPERATURE");
    double multiplier(-1);
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
      resName, valName, *basis, *ir, multiplier));
    this->template registerEvaluator<EvalT>(fm, op);
  }
}

// ***********************************************************************

#endif
