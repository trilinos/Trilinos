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
				      const panzer::FieldLibrary& fl,
				      const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  // ********************
  // MeshCoords Equation
  // ********************

  std::string coord_names[3] = {"X","Y","Z"};
  std::string dof_names[3] = {"COORDX","COORDY","COORDZ"};

  Teuchos::RCP<panzer::IntegrationRule> ir = this->getIntRuleForDOF(dof_names[0]); 
  Teuchos::RCP<panzer::BasisIRLayout> basis = this->getBasisIRLayoutForDOF(dof_names[0]); 

  // Transient Operator
  if (this->buildTransientSupport()) {
    for(int i=0;i<dimension_;i++) {
      ParameterList p("Transient Residual");
      p.set("Residual Name", "RESIDUAL_"+dof_names[i]+"_TRANSIENT_OP");
      p.set("Value Name", "DXDT_"+dof_names[i]);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);
  
      RCP< PHX::Evaluator<panzer::Traits> > op = 
        rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));
      
      fm.template registerEvaluator<EvalT>(op);
    }
  }
  else {
    TEUCHOS_ASSERT(false); // what does this even mean?
  }

  for(int i=0;i<dimension_;i++) {
    ParameterList p("Momentum Source Residual");
    p.set("Residual Name", "RESIDUAL_"+dof_names[i]+"_SOURCE_OP");
    p.set("Value Name", "NODE_VELOCITY_"+coord_names[i]);
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", -1.0);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

  // Use a sum operator to form the overall residual for the equation
  // - this way we avoid loading each operator separately into the
  // global residual and Jacobian
  {
    for(int i=0;i<dimension_;i++) {
      std::vector<std::string> residual_operator_names;

      residual_operator_names.push_back("RESIDUAL_"+dof_names[i]+"_TRANSIENT_OP");
      residual_operator_names.push_back("RESIDUAL_"+dof_names[i]+"_SOURCE_OP");

      this->buildAndRegisterResidualSummationEvalautor(fm,dof_names[i],residual_operator_names);
    }
  }
}

// ***********************************************************************

#endif
