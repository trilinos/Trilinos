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

#ifndef PANZER_BCSTRATEGY_DIRICHLET_DEFAULT_IMPL_IMPL_HPP
#define PANZER_BCSTRATEGY_DIRICHLET_DEFAULT_IMPL_IMPL_HPP

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

// Evaluators
#include "Panzer_Dirichlet_Residual.hpp"
#include "Panzer_GatherSolution_Epetra.hpp"
#include "Panzer_GatherBasisCoordinates.hpp"
#include "Panzer_ScatterDirichletResidual_Epetra.hpp"

// ***********************************************************************
template <typename EvalT>
panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>::
BCStrategy_Dirichlet_DefaultImpl(const panzer::BC& bc,
				 const Teuchos::RCP<panzer::GlobalData>& global_data) :
  panzer::BCStrategy<EvalT>(bc),
  panzer::GlobalDataAcceptorDefaultImpl(global_data)
{

}

// ***********************************************************************
template <typename EvalT>
panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>::
~BCStrategy_Dirichlet_DefaultImpl()
{

}

// ***********************************************************************
template <typename EvalT>
void panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>::
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

  // **************************
  // Coordinates for basis functions (no integration points needed)
  // **************************
  {
    const std::map<std::string,Teuchos::RCP<panzer::PureBasis> > & bases = pb.getBases();
    for (std::map<std::string,Teuchos::RCP<panzer::PureBasis> >::const_iterator it=bases.begin();
         it!=bases.end();it++) {

       // add basis coordinates
       RCP< PHX::Evaluator<panzer::Traits> > basis_op
          = rcp(new panzer::GatherBasisCoordinates<EvalT,panzer::Traits>(*it->second));
       fm.template registerEvaluator<EvalT>(basis_op);
    }
  }

  // Gather
  for (vector<string>::const_iterator dof_name = required_dof_names.begin();
       dof_name != required_dof_names.end(); ++dof_name) {
    
    ParameterList p("BC Gather");
    
    RCP<vector<string> > gather_names_vec = rcp(new vector<string>);
    gather_names_vec->push_back(*dof_name);
    
    p.set("DOF Names", gather_names_vec);
    p.set("Indexer Names", gather_names_vec);
    
    const vector<pair<string,RCP<panzer::PureBasis> > >& dofBasisPair = pb.getProvidedDOFs();
    RCP<panzer::PureBasis> basis;
    for (vector<pair<string,RCP<panzer::PureBasis> > >::const_iterator it = 
	   dofBasisPair.begin(); it != dofBasisPair.end(); ++it) {
      if (it->first == *dof_name)
	basis = it->second;
    }
    
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(basis), std::runtime_error,
		       "Error the name \"" << *dof_name
		       << "\" is not a valid DOF for the boundary condition:\n"
		       << this->m_bc << "\n");
    
    p.set("Basis", basis);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildGather<EvalT>(p);
    
    fm.template registerEvaluator<EvalT>(op);
  }
  
  // Dirichlet Residual: residual = dof_value - target_value
  map<string,string>::const_iterator res_to_target = residual_to_target_field_map.begin();
  for (map<string,string>::const_iterator res_to_dof = residual_to_dof_names_map.begin();
       res_to_dof != residual_to_dof_names_map.end(); ++res_to_dof, ++res_to_target) {

    ParameterList p("Dirichlet Residual: "+res_to_dof->first + " to " + res_to_dof->second);
    p.set("Residual Name", res_to_dof->first);
    p.set("DOF Name", res_to_dof->second);
    p.set("Value Name", res_to_target->second);

    const vector<pair<string,RCP<panzer::PureBasis> > >& dofBasisPair = pb.getProvidedDOFs();
    RCP<panzer::PureBasis> basis;
    for (vector<pair<string,RCP<panzer::PureBasis> > >::const_iterator it = 
	   dofBasisPair.begin(); it != dofBasisPair.end(); ++it) {
      if (it->first == res_to_dof->second)
	basis = it->second;
    }
    
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(basis), std::runtime_error,
		       "Error the name \"" << res_to_dof->second
		       << "\" is not a valid DOF for the boundary condition:\n"
		       << this->m_bc << "\n");
    
    p.set("Data Layout", basis->functional);

    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::DirichletResidual<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);

  }

  // Scatter
 
  for (map<string,string>::const_iterator res_to_dof = residual_to_dof_names_map.begin();
       res_to_dof != residual_to_dof_names_map.end(); ++res_to_dof) {

    ParameterList p("Scatter: "+res_to_dof->first + " to " + res_to_dof->second);
    
    // Set name
    string scatter_field_name = "Dummy Scatter: " + this->m_bc.identifier() + res_to_dof->first; 
    p.set("Scatter Name", scatter_field_name);

    // Set basis
    const vector<pair<string,RCP<panzer::PureBasis> > >& dofBasisPair = pb.getProvidedDOFs();
    RCP<panzer::PureBasis> basis;
    for (vector<pair<string,RCP<panzer::PureBasis> > >::const_iterator it = 
	   dofBasisPair.begin(); it != dofBasisPair.end(); ++it) {
      if (it->first == res_to_dof->second)
	basis = it->second;
    }
    
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(basis), std::runtime_error,
		       "Error the name \"" << res_to_dof->second
		       << "\" is not a valid DOF for the boundary condition:\n"
		       << this->m_bc << "\n");
    
    p.set("Basis", basis);

    RCP<vector<string> > residual_names = rcp(new vector<string>);
    residual_names->push_back(res_to_dof->first);
    p.set("Dependent Names", residual_names);

    RCP<map<string,string> > names_map = rcp(new map<string,string>);
    names_map->insert(*res_to_dof);
    p.set("Dependent Map", names_map);
    
    TEUCHOS_TEST_FOR_EXCEPTION(!pb.cellData().isSide(), std::logic_error,
		       "Error - physics block is not a side set!");
    
    p.set<int>("Side Subcell Dimension", 
	       pb.getBaseCellTopology().getDimension() - 1);
    p.set<int>("Local Side ID", pb.cellData().side());

    RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildScatterDirichlet<EvalT>(p);
      // rcp(new panzer::ScatterDirichletResidual_Epetra<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
    
    // Require variables
    {
      using panzer::Dummy;
      PHX::Tag<typename EvalT::ScalarT> tag(scatter_field_name, 
					    rcp(new PHX::MDALayout<Dummy>(0)));
      fm.template requireField<EvalT>(tag);
    }
  
  }

}

// ***********************************************************************

#endif
