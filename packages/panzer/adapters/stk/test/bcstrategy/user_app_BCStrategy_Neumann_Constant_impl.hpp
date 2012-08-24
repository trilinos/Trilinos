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

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Panzer_PhysicsBlock.hpp"

#include "Panzer_BasisIRLayout.hpp"

// Evaluators
#include "Panzer_ConstantFlux.hpp"

#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

// ***********************************************************************
template <typename EvalT>
user_app::BCStrategy_Neumann_Constant<EvalT>::
BCStrategy_Neumann_Constant(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) :
  panzer::BCStrategy_Neumann_DefaultImpl<EvalT>(bc,global_data)
{

}

// ***********************************************************************
template <typename EvalT>
void user_app::BCStrategy_Neumann_Constant<EvalT>::
setup(const panzer::PhysicsBlock& side_pb,
      const Teuchos::ParameterList& user_data)
{
  // need the dof value to form the residual
  this->requireDOFGather(this->m_bc.equationSetName());

  const std::string residual_name = "Residual_" + this->m_bc.identifier();
  const std::string dof_name = this->m_bc.equationSetName();
  const std::string flux_name = "Constant_" + this->m_bc.equationSetName();
  const int integration_order = this->m_bc.params()->template get<int>("Integration Order");

  this->addResidualContribution(residual_name,dof_name,flux_name,integration_order,side_pb);
}

// ***********************************************************************
template <typename EvalT>
void user_app::BCStrategy_Neumann_Constant<EvalT>::
buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
			   const panzer::PhysicsBlock& pb,
			   const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
			   const Teuchos::ParameterList& models,
			   const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  const std::vector<boost::tuples::tuple<std::string,std::string,std::string,int,Teuchos::RCP<panzer::PureBasis>,Teuchos::RCP<panzer::IntegrationRule> > > data = this->getResidualContributionData();

  std::string flux_name = data[0].get<2>();
  Teuchos::RCP<panzer::IntegrationRule> ir = data[0].get<5>();

  // provide a constant flux target value to map into residual
  {
    ParameterList p("BC Constant Neumann");
    p.set("Flux Field Name", flux_name);
    p.set("Data Layout", ir->dl_vector);
    TEUCHOS_ASSERT(this->m_bc.params()->isSublist("Flux Values"));
    p.sublist("Flux Values") = this->m_bc.params()->sublist("Flux Values");
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::ConstantFlux<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
  }

}

// ***********************************************************************
template <typename EvalT>
void user_app::BCStrategy_Neumann_Constant<EvalT>::
postRegistrationSetup(typename panzer::Traits::SetupData d,
		      PHX::FieldManager<panzer::Traits>& vm)
{
  
}


// ***********************************************************************
template <typename EvalT>
void user_app::BCStrategy_Neumann_Constant<EvalT>::
evaluateFields(typename panzer::Traits::EvalData d)
{
  
}
