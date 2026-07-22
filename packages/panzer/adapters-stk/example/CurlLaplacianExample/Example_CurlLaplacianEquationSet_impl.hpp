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
#include "Panzer_Integrator_BasisTimesVector.hpp"
#include "Panzer_Integrator_CurlBasisDotVector.hpp"
#include "Panzer_ScalarToVector.hpp"
#include "Panzer_Sum.hpp"
#include "Panzer_Constant.hpp"

// ***********************************************************************
template <typename EvalT>
Example::CurlLaplacianEquationSet<EvalT>::
CurlLaplacianEquationSet(const Teuchos::RCP<Teuchos::ParameterList>& params,
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
   // In this case "EFIELD" is the lone field.  We also name the curl
   // for this field. These names automatically generate evaluators for "EFIELD"
   // and "CURL_EFIELD" gathering the basis coefficients of "EFIELD" and
   // the values of the EFIELD and CURL_EFIELD fields at quadrature points.
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

   this->addDOF("EFIELD",basis_type,basis_order,integration_order);
   this->addDOFCurl("EFIELD");
   if (this->buildTransientSupport())
     this->addDOFTimeDerivative("EFIELD");

   // ********************
   // Build Basis Functions and Integration Rules
   // ********************
   
   this->addClosureModel(model_id);

   this->setupDOFs();
}

// ***********************************************************************
template <typename EvalT>
void Example::CurlLaplacianEquationSet<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				      const panzer::FieldLibrary& /* fl */,
				      const Teuchos::ParameterList& /* user_data */) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
   
  Teuchos::RCP<panzer::IntegrationRule> ir = this->getIntRuleForDOF("EFIELD");
  Teuchos::RCP<panzer::BasisIRLayout> basis = this->getBasisIRLayoutForDOF("EFIELD"); 

  // ********************
  // Energy Equation
  // ********************

  // Transient Operator: Assembles \int \dot{T} v
  if (this->buildTransientSupport()) {
    ParameterList p("Transient Residual");
    p.set("Residual Name", "RESIDUAL_EFIELD_TRANSIENT_OP"); // we are defining the name of this operator
    p.set("Value Name", "DXDT_EFIELD"); // this field is constructed by the panzer library
//    p.set("Test Field Name", "EFIELD"); 
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", 1.0);

    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_BasisTimesVector<EvalT,panzer::Traits>(p));
    
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // Diffusion Operator: Assembles \int \nabla T \cdot \nabla v
  {
    using panzer::EvaluatorStyle;
    using panzer::Integrator_CurlBasisDotVector;
    using panzer::Traits;
    using PHX::Evaluator;
    using std::string;
    string resName("RESIDUAL_EFIELD"), valName("CURL_EFIELD");
    double thermalConductivity(1.0), multiplier(-thermalConductivity);
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_CurlBasisDotVector<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
      resName, valName, *basis, *ir, multiplier));
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // Mass operator
  {   
    ParameterList p("Source Residual");
    p.set("Residual Name", "RESIDUAL_EFIELD_MASS_OP");
    p.set("Value Name", "EFIELD"); 
//    p.set("Test Field Name", "EFIELD"); 
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", -1.0);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_BasisTimesVector<EvalT,panzer::Traits>(p));
    
    this->template registerEvaluator<EvalT>(fm, op);
  }
  
  // Source Operator
  {   
    ParameterList p("Source Residual");
    p.set("Residual Name", "RESIDUAL_EFIELD_SOURCE_OP");
    p.set("Value Name", "SOURCE_EFIELD"); // this field must be provided by the closure model factory
                                               // and specified by the user
//    p.set("Test Field Name", "EFIELD"); 
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", 1.0);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_BasisTimesVector<EvalT,panzer::Traits>(p));
    
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // Use a sum operator to form the overall residual for the equation
  {
    std::vector<std::string> sum_names;
    
    // these are the names of the residual values to sum together
    sum_names.push_back("RESIDUAL_EFIELD_MASS_OP");
    sum_names.push_back("RESIDUAL_EFIELD_SOURCE_OP");
    if (this->buildTransientSupport())
      sum_names.push_back("RESIDUAL_EFIELD_TRANSIENT_OP");

    this->buildAndRegisterResidualSummationEvaluator(fm,"EFIELD",sum_names);
  }

}

// ***********************************************************************

#endif
