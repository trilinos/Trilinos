// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_BCSTRATEGY_WEAKDIRICHLET_DEFAULT_IMPL_IMPL_HPP
#define PANZER_BCSTRATEGY_WEAKDIRICHLET_DEFAULT_IMPL_IMPL_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Panzer_PhysicsBlock.hpp"

#include "Panzer_PureBasis.hpp"

#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include <sstream>

// Evaluators
#include "Panzer_WeakDirichlet_Residual.hpp"

#ifdef PANZER_HAVE_EPETRA_STACK
#include "Panzer_GatherSolution_Epetra.hpp"
#include "Panzer_ScatterResidual_Epetra.hpp"
#endif

#include "Panzer_Normals.hpp"

// ***********************************************************************
template <typename EvalT>
panzer::BCStrategy_WeakDirichlet_DefaultImpl<EvalT>::
BCStrategy_WeakDirichlet_DefaultImpl(const panzer::BC& bc,
			       const Teuchos::RCP<panzer::GlobalData>& global_data) :
  panzer::BCStrategy<EvalT>(bc),
  panzer::GlobalDataAcceptorDefaultImpl(global_data)
{

}

// ***********************************************************************
template <typename EvalT>
panzer::BCStrategy_WeakDirichlet_DefaultImpl<EvalT>::
~BCStrategy_WeakDirichlet_DefaultImpl()
{

}

// ***********************************************************************
template <typename EvalT>
void panzer::BCStrategy_WeakDirichlet_DefaultImpl<EvalT>::
buildAndRegisterGatherScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
			                const panzer::PhysicsBlock& pb,
				        const panzer::LinearObjFactory<panzer::Traits> & lof,
					const Teuchos::ParameterList& user_data) const
{
  buildAndRegisterGatherAndOrientationEvaluators(fm,pb,lof,user_data);
  buildAndRegisterScatterEvaluators(fm,pb,lof,user_data);
}

// ***********************************************************************
template <typename EvalT>
void panzer::BCStrategy_WeakDirichlet_DefaultImpl<EvalT>::
buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
			                       const panzer::PhysicsBlock& pb,
				               const panzer::LinearObjFactory<panzer::Traits> & lof,
				               const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::vector;
  using std::map;
  using std::string;
  using std::pair;

  // Gather
  pb.buildAndRegisterGatherAndOrientationEvaluators(fm,lof,user_data);

  // Iterate over each residual contribution
  for (vector<std::tuple<std::string,std::string,std::string,int,Teuchos::RCP<panzer::PureBasis>,Teuchos::RCP<panzer::IntegrationRule> > >::const_iterator eq =
	 m_residual_contributions.begin(); eq != m_residual_contributions.end(); ++eq) {

    const string& residual_name = std::get<0>(*eq);
    const string& dof_name = std::get<1>(*eq);
    const string& flux_name = std::get<2>(*eq);
    //const int& integration_order = std::get<3>(*eq);
    const RCP<const panzer::PureBasis> basis = std::get<4>(*eq);
    const RCP<const panzer::IntegrationRule> ir = std::get<5>(*eq);

    // Normals evaluator
    {
      std::stringstream s;
      s << "Side Normal:" << pb.cellData().side();
      ParameterList p(s.str());
      p.set<std::string>("Name","Side Normal");
      p.set<int>("Side ID",pb.cellData().side());
      p.set< Teuchos::RCP<panzer::IntegrationRule> >("IR", Teuchos::rcp_const_cast<panzer::IntegrationRule>(ir));
      p.set<bool>("Normalize",true);

      RCP< PHX::Evaluator<panzer::Traits> > op = rcp(new panzer::Normals<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
    }

    // Neumann Residual evaluator: residual += phi n dot flux
    {
      ParameterList p("Neumann Residual: " + residual_name + " to DOF: " + dof_name);
      p.set("Residual Name", residual_name);
      p.set("DOF Name",dof_name);
      p.set("Flux Name", flux_name);
      p.set("Normal Name", "Side Normal");
      p.set("Basis", basis);
      p.set("IR", ir);

      RCP< PHX::Evaluator<panzer::Traits> > op =
	rcp(new panzer::WeakDirichletResidual<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
    }

  }
}

// ***********************************************************************
template <typename EvalT>
void panzer::BCStrategy_WeakDirichlet_DefaultImpl<EvalT>::
buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
  		                  const panzer::PhysicsBlock& /* pb */,
			          const panzer::LinearObjFactory<panzer::Traits> & lof,
			  	  const Teuchos::ParameterList& /* user_data */) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::vector;
  using std::map;
  using std::string;
  using std::pair;

  // Iterate over each residual contribution
  for (vector<std::tuple<std::string,std::string,std::string,int,Teuchos::RCP<panzer::PureBasis>,Teuchos::RCP<panzer::IntegrationRule> > >::const_iterator eq =
	 m_residual_contributions.begin(); eq != m_residual_contributions.end(); ++eq) {

    const string& residual_name = std::get<0>(*eq);
    const string& dof_name = std::get<1>(*eq);
    const RCP<const panzer::PureBasis> basis = std::get<4>(*eq);
    const RCP<const panzer::IntegrationRule> ir = std::get<5>(*eq);

    // Scatter evaluator
    {
      ParameterList p("Scatter: "+ residual_name + " to " + dof_name);

      // Set name
      string scatter_field_name = "Dummy Scatter: " + this->m_bc.identifier() + residual_name;
      p.set("Scatter Name", scatter_field_name);
      p.set("Basis", basis);

      RCP<vector<string> > residual_names = rcp(new vector<string>);
      residual_names->push_back(residual_name);
      p.set("Dependent Names", residual_names);

      RCP<map<string,string> > names_map = rcp(new map<string,string>);
      names_map->insert(std::pair<string,string>(residual_name,dof_name));
      p.set("Dependent Map", names_map);

      RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildScatter<EvalT>(p);

      this->template registerEvaluator<EvalT>(fm, op);

      // Require variables
      {
	using panzer::Dummy;
	PHX::Tag<typename EvalT::ScalarT> tag(scatter_field_name,
					      rcp(new PHX::MDALayout<Dummy>(0)));
	fm.template requireField<EvalT>(tag);
      }

    }  // end of Scatter

  }
}

// ***********************************************************************
template <typename EvalT>
void panzer::BCStrategy_WeakDirichlet_DefaultImpl<EvalT>::
requireDOFGather(const std::string required_dof_name)
{
  m_required_dof_names.push_back(required_dof_name);
}

// ***********************************************************************
template <typename EvalT>
void panzer::BCStrategy_WeakDirichlet_DefaultImpl<EvalT>::
addResidualContribution(const std::string residual_name,
			const std::string dof_name,
			const std::string flux_name,
			const int integration_order,
			const panzer::PhysicsBlock& side_pb)
{
  Teuchos::RCP<panzer::PureBasis> basis = this->getBasis(dof_name,side_pb);

  Teuchos::RCP<panzer::IntegrationRule> ir = buildIntegrationRule(integration_order,side_pb);

  m_residual_contributions.push_back(std::make_tuple(residual_name,
						     dof_name,
						     flux_name,
						     integration_order,
						     basis,
						     ir));
}

// ***********************************************************************
template <typename EvalT>
const std::vector<std::tuple<std::string,std::string,std::string,int,Teuchos::RCP<panzer::PureBasis>,Teuchos::RCP<panzer::IntegrationRule> > >
panzer::BCStrategy_WeakDirichlet_DefaultImpl<EvalT>::getResidualContributionData() const
{
  return m_residual_contributions;
}

// ***********************************************************************
template <typename EvalT>
Teuchos::RCP<panzer::PureBasis>
panzer::BCStrategy_WeakDirichlet_DefaultImpl<EvalT>::
getBasis(const std::string dof_name,const panzer::PhysicsBlock& side_pb) const
{
  const std::vector<std::pair<std::string,Teuchos::RCP<panzer::PureBasis> > >& dofBasisPair = side_pb.getProvidedDOFs();
  Teuchos::RCP<panzer::PureBasis> basis;
  for (std::vector<std::pair<std::string,Teuchos::RCP<panzer::PureBasis> > >::const_iterator it =
	 dofBasisPair.begin(); it != dofBasisPair.end(); ++it) {
    if (it->first == dof_name)
      basis = it->second;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(is_null(basis), std::runtime_error,
			     "Error the name \"" << dof_name
			     << "\" is not a valid DOF for the boundary condition:\n"
			     << this->m_bc << "\n");

  return basis;
}

// ***********************************************************************
template <typename EvalT>
Teuchos::RCP<panzer::IntegrationRule>
panzer::BCStrategy_WeakDirichlet_DefaultImpl<EvalT>::
buildIntegrationRule(const int integration_order,const panzer::PhysicsBlock& side_pb) const
{
  TEUCHOS_ASSERT(side_pb.cellData().isSide());
  Teuchos::RCP<panzer::IntegrationRule> ir = Teuchos::rcp(new panzer::IntegrationRule(integration_order,side_pb.cellData()));
  return ir;
}

// ***********************************************************************
template <typename EvalT>
const panzer::BC
panzer::BCStrategy_WeakDirichlet_DefaultImpl<EvalT>::bc() const
{
  return this->m_bc;
}

// ***********************************************************************

#endif
