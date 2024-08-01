// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef USER_APP_EQUATIONSET_MESHCOORDS_T_HPP
#define USER_APP_EQUATIONSET_MESHCOORDS_T_HPP

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
#include "Panzer_ZeroContributedField.hpp"
#include "user_app_Convection.hpp"

// ***********************************************************************
template <typename EvalT>
user_app::EquationSet_MeshCoords<EvalT>::
EquationSet_MeshCoords(const Teuchos::RCP<Teuchos::ParameterList>& params,
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
    valid_parameters.set("Integration Order",-1,"Order of the integration rule");

    params->validateParametersAndSetDefaults(valid_parameters);
  }

  std::string basis_type = "HGrad"; // use nodal linears
  int basis_order = 1;
  std::string model_id = params->get<std::string>("Model ID");
  int integration_order = params->get<int>("Integration Order");

  dimension_ = cell_data.baseCellDimension();

  // ********************
  // Setup DOFs and closure models
  // ********************
  std::string dof_names[3] = {"COORDX","COORDY","COORDZ"};
  std::vector<std::string> coord;
  for(int i=0;i<dimension_;i++) {
    std::string dof_name = dof_names[i];
    coord.push_back(dof_name);

    this->addDOF(dof_name,basis_type,basis_order,integration_order);
    if (this->buildTransientSupport())
      this->addDOFTimeDerivative(dof_name);
    else {
      TEUCHOS_ASSERT(false); // what does this even mean?
    }
  }

  this->setCoordinateDOFs(coord);

  this->addClosureModel(model_id);

  this->setupDOFs();
}

// ***********************************************************************
template <typename EvalT>
void user_app::EquationSet_MeshCoords<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				      const panzer::FieldLibrary& /* fl */,
				      const Teuchos::ParameterList& /* user_data */) const
{
  using panzer::BasisIRLayout;
  using panzer::EvaluatorStyle;
  using panzer::IntegrationRule;
  using panzer::Integrator_BasisTimesScalar;
  using panzer::Traits;
  using panzer::ZeroContributedField;
  using PHX::Evaluator;
  using std::string;
  using std::vector;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  // ********************
  // MeshCoords Equation
  // ********************

  string coord_names[3] = {"X","Y","Z"};
  string dof_names[3] = {"COORDX","COORDY","COORDZ"};

  RCP<IntegrationRule> ir = this->getIntRuleForDOF(dof_names[0]); 
  RCP<BasisIRLayout> basis = this->getBasisIRLayoutForDOF(dof_names[0]); 

  // Initialize the residual to zero.
  for (int i(0); i < dimension_; ++i)
  {
    string fieldName("RESIDUAL_" + dof_names[i]);
    auto layout = ir->dl_scalar;
    RCP<Evaluator<Traits>> op =
      rcp(new ZeroContributedField<EvalT, Traits>(fieldName, *layout));
    this->template registerEvaluator<EvalT>(fm, op);
  } // end loop over dimension_

  // Transient Operator
  if (this->buildTransientSupport())
  {
    for (int i(0); i < dimension_; ++i)
    {
      string resName("RESIDUAL_" + dof_names[i]),
             valName("DXDT_" + dof_names[i]);
      double multiplier(1);
      RCP<Evaluator<Traits>> op = rcp(new
        Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
        resName, valName, *basis, *ir, multiplier));
      this->template registerEvaluator<EvalT>(fm, op);
    }
  }
  else
  {
    TEUCHOS_ASSERT(false); // what does this even mean?
  }

  for (int i(0); i < dimension_; ++i)
  {
    string resName("RESIDUAL_" + dof_names[i]),
           valName("NODE_VELOCITY_" + coord_names[i]);
    double multiplier(-1);
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
      resName, valName, *basis, *ir, multiplier));
    this->template registerEvaluator<EvalT>(fm, op);
  }
}

// ***********************************************************************

#endif
