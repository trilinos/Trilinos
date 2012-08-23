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

#ifndef PANZER_BCSTRATEGY_NEUMANN_DEFAULT_IMPL_IMPL_HPP
#define PANZER_BCSTRATEGY_NEUMANN_DEFAULT_IMPL_IMPL_HPP

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
#include "Panzer_Neumann_Residual.hpp"
#include "Panzer_GatherSolution_Epetra.hpp"
#include "Panzer_ScatterResidual_Epetra.hpp"
#include "Panzer_Normals.hpp"

// ***********************************************************************
template <typename EvalT>
panzer::BCStrategy_Neumann_DefaultImpl<EvalT>::
BCStrategy_Neumann_DefaultImpl(const panzer::BC& bc,
			       const Teuchos::RCP<panzer::GlobalData>& global_data) :
  panzer::BCStrategy<EvalT>(bc),
  panzer::GlobalDataAcceptorDefaultImpl(global_data)
{
  
}

// ***********************************************************************
template <typename EvalT>
panzer::BCStrategy_Neumann_DefaultImpl<EvalT>::
~BCStrategy_Neumann_DefaultImpl()
{
  
}

// ***********************************************************************
template <typename EvalT>
void panzer::BCStrategy_Neumann_DefaultImpl<EvalT>::
buildAndRegisterGatherScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
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
  for (vector<boost::tuples::tuple<std::string,std::string,std::string,int,Teuchos::RCP<panzer::PureBasis>,Teuchos::RCP<panzer::IntegrationRule> > >::const_iterator eq = 
	 m_residual_contributions.begin(); eq != m_residual_contributions.end(); ++eq) {
    
    const string& residual_name = eq->get<0>();
    const string& dof_name = eq->get<1>();
    const string& flux_name = eq->get<2>();
    //const int& integration_order = eq->get<3>();
    const RCP<const panzer::PureBasis> basis = eq->get<4>();
    const RCP<const panzer::IntegrationRule> ir = eq->get<5>();

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
      
      fm.template registerEvaluator<EvalT>(op);
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
	rcp(new panzer::NeumannResidual<EvalT,panzer::Traits>(p));
    
      fm.template registerEvaluator<EvalT>(op);
    }

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
    
      fm.template registerEvaluator<EvalT>(op);
      
      // Require variables
      {
	using panzer::Dummy;
	PHX::Tag<typename EvalT::ScalarT> tag(scatter_field_name, 
					      rcp(new PHX::MDALayout<Dummy>(0)));
	fm.template requireField<EvalT>(tag);
      }

    }  // end of Scatter

  }  // end of residual contribution 

}

// ***********************************************************************
template <typename EvalT>
void panzer::BCStrategy_Neumann_DefaultImpl<EvalT>::
requireDOFGather(const std::string required_dof_name)
{
  m_required_dof_names.push_back(required_dof_name);
}

// ***********************************************************************
template <typename EvalT>
void panzer::BCStrategy_Neumann_DefaultImpl<EvalT>::
addResidualContribution(const std::string residual_name,
			const std::string dof_name,
			const std::string flux_name,
			const int integration_order,
			const panzer::PhysicsBlock& side_pb)
{
  Teuchos::RCP<panzer::PureBasis> basis = this->getBasis(dof_name,side_pb);

  Teuchos::RCP<panzer::IntegrationRule> ir = buildIntegrationRule(integration_order,side_pb);
  
  m_residual_contributions.push_back(boost::tuples::make_tuple(residual_name,
							       dof_name,
							       flux_name,
							       integration_order,
							       basis,
							       ir));
}

// ***********************************************************************
template <typename EvalT>
const std::vector<boost::tuples::tuple<std::string,std::string,std::string,int,Teuchos::RCP<panzer::PureBasis>,Teuchos::RCP<panzer::IntegrationRule> > >
panzer::BCStrategy_Neumann_DefaultImpl<EvalT>::getResidualContributionData() const
{
  return m_residual_contributions;
}

// ***********************************************************************
template <typename EvalT>
Teuchos::RCP<panzer::PureBasis> 
panzer::BCStrategy_Neumann_DefaultImpl<EvalT>::
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
panzer::BCStrategy_Neumann_DefaultImpl<EvalT>::
buildIntegrationRule(const int integration_order,const panzer::PhysicsBlock& side_pb) const
{
  TEUCHOS_ASSERT(side_pb.cellData().isSide());
  Teuchos::RCP<panzer::IntegrationRule> ir = Teuchos::rcp(new panzer::IntegrationRule(integration_order,side_pb.cellData()));
  return ir;
}

// ***********************************************************************

#endif
