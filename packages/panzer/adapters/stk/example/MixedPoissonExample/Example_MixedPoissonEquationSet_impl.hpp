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
    
    params->validateParametersAndSetDefaults(valid_parameters);
  }
  
  std::string basis_type = "HGrad";
  std::string grad_basis_type = "HDiv";
  int basis_order = 1;
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

  this->addDOF("GRADPHI_FIELD",grad_basis_type,basis_order,integration_order);
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
				      const panzer::FieldLibrary& fl,
				      const Teuchos::ParameterList& user_data) const
{
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
    ParameterList p("Source Residual");
    p.set("Residual Name", "RESIDUAL_GRADPHI_FIELD_DIFFUSION_OP");
    p.set("Value Name", "DIV_GRADPHI_FIELD"); 
    p.set("Test Field Name", "GRADPHI_FIELD"); 
    p.set("Basis", basis_v);
    p.set("IR", ir);
    p.set("Multiplier", 1.0);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_DivBasisTimesScalar<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
  }

  // Source operator (-f,-\nabla\cdot w)
  {   
    ParameterList p("Source Residual");
    p.set("Residual Name", "RESIDUAL_GRADPHI_FIELD_MASS_OP");
    p.set("Value Name", "SOURCE"); 
    p.set("Test Field Name", "GRADPHI_FIELD"); 
    p.set("Basis", basis_v);
    p.set("IR", ir);
    p.set("Multiplier", -1.0); // scale by the right hand side 
                                            // when phi = sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_DivBasisTimesScalar<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
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
    
    fm.template registerEvaluator<EvalT>(op);
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
    
    fm.template registerEvaluator<EvalT>(op);
  }
  
  // Use a sum operator to form the overall residual for the equation
  {
    std::vector<std::string> sum_names;
    
    // these are the names of the residual values to sum together
    sum_names.push_back("RESIDUAL_GRADPHI_FIELD_DIFFUSION_OP");
    sum_names.push_back("RESIDUAL_GRADPHI_FIELD_MASS_OP");

    this->buildAndRegisterResidualSummationEvalautor(fm,"GRADPHI_FIELD",sum_names);
  }

  // Use a sum operator to form the overall residual for the equation
  {
    std::vector<std::string> sum_names;
    
    // these are the names of the residual values to sum together
    sum_names.push_back("RESIDUAL_PHI_DIFFUSION_OP");
    sum_names.push_back("RESIDUAL_PHI_MASS_OP");

    this->buildAndRegisterResidualSummationEvalautor(fm,"PHI",sum_names);
  }

}

// ***********************************************************************

#endif
