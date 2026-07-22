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

#include "Panzer_PureBasis.hpp"

// Evaluators
#include "Panzer_Integrator_BasisTimesScalar.hpp"
#include "Panzer_Constant.hpp"

#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

// ***********************************************************************
template <typename EvalT>
Example::BCStrategy_Neumann_Constant<EvalT>::
BCStrategy_Neumann_Constant(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) :
  panzer::BCStrategy_Neumann_DefaultImpl<EvalT>(bc,global_data)
{
  TEUCHOS_ASSERT(this->m_bc.strategy() == "Neumann Constant");
}

// ***********************************************************************
template <typename EvalT>
void Example::BCStrategy_Neumann_Constant<EvalT>::
setup(const panzer::PhysicsBlock& side_pb,
      const Teuchos::ParameterList& /* user_data */)
{
  using Teuchos::RCP;
  using std::vector;
  using std::string;
  using std::pair;

  // obtain the dof name
  const string dof_name = this->m_bc.equationSetName();

  // need the dof value to form the residual
  this->requireDOFGather(dof_name);

  // unique residual name
  const string residual_name = "Residual_" + dof_name;
  const string flux_name = "Constant_Flux";

  const std::map<int,RCP< panzer::IntegrationRule > >& ir = side_pb.getIntegrationRules();
  TEUCHOS_ASSERT(ir.size() == 1); 
  
  const int integration_order = ir.begin()->second->order();

  this->addResidualContribution(residual_name,dof_name,flux_name,integration_order,side_pb);

}

// ***********************************************************************
template <typename EvalT>
void Example::BCStrategy_Neumann_Constant<EvalT>::
buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
			   const panzer::PhysicsBlock& pb,
			   const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& /* factory */,
			   const Teuchos::ParameterList& /* models */,
			   const Teuchos::ParameterList& /* user_data */) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
using std::string; 

  const std::vector<std::tuple<string,string,string,int,Teuchos::RCP<panzer::PureBasis>,
        Teuchos::RCP<panzer::IntegrationRule> > > data = this->getResidualContributionData();

  string residual_name = std::get<0>(data[0]);
  string dof_name = std::get<1>(data[0]);
  string flux_name = std::get<2>(data[0]);

  RCP<panzer::IntegrationRule> ir = std::get<5>(data[0]);
  RCP<const panzer::FieldLayoutLibrary> fll = pb.getFieldLibrary()->buildFieldLayoutLibrary(*ir);
  RCP<panzer::BasisIRLayout> basis = fll->lookupLayout(dof_name);

  // provide a constant flux target value to map into residual
  {
    ParameterList p("Constant Neumann BC");
    p.set("Data Layout", ir->dl_scalar);
    p.set("Name", flux_name);
    p.set("Value",  this->m_bc.params()->template get<double>("Value"));
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Constant<EvalT,panzer::Traits>(p));
    
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // add contribution to the residual 
  {
    using panzer::EvaluatorStyle;
    using panzer::Integrator_BasisTimesScalar;
    using panzer::Traits;
    using PHX::Evaluator;
    double multiplier(1);
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::EVALUATES,
      residual_name, flux_name, *basis, *ir, multiplier));
    this->template registerEvaluator<EvalT>(fm, op);
  }
}

// ***********************************************************************
template <typename EvalT>
void Example::BCStrategy_Neumann_Constant<EvalT>::
buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
			                       const panzer::PhysicsBlock& pb,
				               const panzer::LinearObjFactory<panzer::Traits> & lof,
				               const Teuchos::ParameterList& user_data) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::vector;
  using std::map;
  using std::string;
  using std::pair;

  // Gather
  pb.buildAndRegisterGatherAndOrientationEvaluators(fm,lof,user_data);

}

// ***********************************************************************
template <typename EvalT>
void Example::BCStrategy_Neumann_Constant<EvalT>::
postRegistrationSetup(typename panzer::Traits::SetupData /* d */,
		      PHX::FieldManager<panzer::Traits>& /* vm */)
{
  
}


// ***********************************************************************
template <typename EvalT>
void Example::BCStrategy_Neumann_Constant<EvalT>::
evaluateFields(typename panzer::Traits::EvalData /* d */)
{
  
}
