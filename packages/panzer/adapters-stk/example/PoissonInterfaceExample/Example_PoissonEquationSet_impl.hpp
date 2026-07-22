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
    
    valid_parameters.set("Model ID","","Closure model id associated with this equation set");
    valid_parameters.set("Basis Type","HGrad","Type of Basis to use");
    valid_parameters.set("Basis Order",1,"Order of the basis");
    valid_parameters.set("Integration Order",-1,"Order of the integration rule");
    valid_parameters.set("DOF Name","DOF Name");
    
    params->validateParametersAndSetDefaults(valid_parameters);
  }
  
  std::string basis_type = params->get<std::string>("Basis Type");
  int basis_order = params->get<int>("Basis Order");
  int integration_order = params->get<int>("Integration Order");
  std::string model_id = params->get<std::string>("Model ID");
  dof_name = params->get<std::string>("DOF Name");

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

   this->addDOF(dof_name,basis_type,basis_order,integration_order);
   this->addDOFGrad(dof_name);
   if (this->buildTransientSupport())
     this->addDOFTimeDerivative(dof_name);

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
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  
  Teuchos::RCP<panzer::IntegrationRule> ir = this->getIntRuleForDOF(dof_name);
  Teuchos::RCP<panzer::BasisIRLayout> basis = this->getBasisIRLayoutForDOF(dof_name); 

  // ********************
  // Energy Equation
  // ********************

  // Transient Operator: Assembles \int \dot{T} v
  if (this->buildTransientSupport()) {
    ParameterList p("Transient Residual");
    p.set("Residual Name", "RESIDUAL_" + dof_name + "_TRANSIENT_OP"); // we are defining the name of this operator
    p.set("Value Name", "DXDT_" + dof_name); // this field is constructed by the panzer library
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", 1.0);

    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));
    
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // Diffusion Operator: Assembles \int \nabla T \cdot \nabla v
  {
    double thermal_conductivity = 1.0;

    ParameterList p("Diffusion Residual");
    p.set("Residual Name", "RESIDUAL_" + dof_name + "_DIFFUSION_OP");
    p.set("Flux Name", "GRAD_" + dof_name); // this field is constructed by the panzer library
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", thermal_conductivity);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_GradBasisDotVector<EvalT,panzer::Traits>(p));

    this->template registerEvaluator<EvalT>(fm, op);
  }
  
  // Source Operator
  {   
    ParameterList p("Source Residual");
    p.set("Residual Name", "RESIDUAL_" + dof_name + "_SOURCE_OP");
    p.set("Value Name", "SOURCE_" + dof_name); // this field must be provided by the closure model factory
                                               // and specified by the user
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", -1.0);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));
    
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // Use a sum operator to form the overall residual for the equation
  {
    std::vector<std::string> sum_names;

    // these are the names of the residual values to sum together
    sum_names.push_back("RESIDUAL_" + dof_name + "_DIFFUSION_OP");
    sum_names.push_back("RESIDUAL_" + dof_name + "_SOURCE_OP");
    if (this->buildTransientSupport())
      sum_names.push_back("RESIDUAL_" + dof_name + "_TRANSIENT_OP");

    this->buildAndRegisterResidualSummationEvaluator(fm,dof_name,sum_names);
  }

}

// ***********************************************************************

#endif
