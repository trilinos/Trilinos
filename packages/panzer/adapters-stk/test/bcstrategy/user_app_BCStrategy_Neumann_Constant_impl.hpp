// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
      const Teuchos::ParameterList& /* user_data */)
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
			   const panzer::PhysicsBlock& /* pb */,
			   const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& /* factory */,
			   const Teuchos::ParameterList& /* models */,
			   const Teuchos::ParameterList& /* user_data */) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  const std::vector<std::tuple<std::string,std::string,std::string,int,Teuchos::RCP<panzer::PureBasis>,Teuchos::RCP<panzer::IntegrationRule> > > data = this->getResidualContributionData();

  std::string flux_name = std::get<2>(data[0]);
  Teuchos::RCP<panzer::IntegrationRule> ir = std::get<5>(data[0]);

  // provide a constant flux target value to map into residual
  {
    ParameterList p("BC Constant Neumann");
    p.set("Flux Field Name", flux_name);
    p.set("Data Layout", ir->dl_vector);
    TEUCHOS_ASSERT(this->m_bc.params()->isSublist("Flux Values"));
    p.sublist("Flux Values") = this->m_bc.params()->sublist("Flux Values");
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::ConstantFlux<EvalT,panzer::Traits>(p));
    
    this->template registerEvaluator<EvalT>(fm, op);
  }

}

// ***********************************************************************
template <typename EvalT>
void user_app::BCStrategy_Neumann_Constant<EvalT>::
postRegistrationSetup(typename panzer::Traits::SetupData /* d */,
		      PHX::FieldManager<panzer::Traits>& /* vm */)
{
  
}


// ***********************************************************************
template <typename EvalT>
void user_app::BCStrategy_Neumann_Constant<EvalT>::
evaluateFields(typename panzer::Traits::EvalData /* d */)
{
  
}
