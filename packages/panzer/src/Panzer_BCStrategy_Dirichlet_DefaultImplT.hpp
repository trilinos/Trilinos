#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Panzer_PhysicsBlock.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Basis.hpp"

// Evaluators
#include "Panzer_DOF.hpp"
#include "Panzer_DOFGradient.hpp"
#include "Panzer_Dirichlet_Constant.hpp"

#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

// ***********************************************************************
template <typename EvalT>
panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>::
BCStrategy_Dirichlet_DefaultImpl(const panzer::BC& bc) :
  panzer::BCStrategy<EvalT>(bc)
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

  // Gather
  for (vector<string>::const_iterator dof_name = required_dof_names.begin();
       dof_name != required_dof_names.end(); ++dof_name) {
    
    ParameterList p("BC Gather");
    
    RCP<vector<string> > gather_names_vec = rcp(new vector<string>);
    gather_names_vec->push_back(*dof_name);
    
    p.set("DOF Names", gather_names_vec);
    p.set("Indexer Names", gather_names_vec);
    
    const vector<pair<string,RCP<panzer::Basis> > >& dofBasisPair = pb.getProvidedDOFs();
    RCP<panzer::Basis> basis;
    for (vector<pair<string,RCP<panzer::Basis> > >::const_iterator it = 
	   dofBasisPair.begin(); it != dofBasisPair.end(); ++it) {
      if (it->first == *dof_name)
	basis = it->second;
    }
    
    TEST_FOR_EXCEPTION(Teuchos::is_null(basis), std::runtime_error,
		       "Error the name \"" << *dof_name
		       << "\" is not a valid DOF for the boundary condition:\n"
		       << this->m_bc << "\n");
    
    p.set("Basis", basis);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildGather<EvalT>(p);
    
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
    const vector<pair<string,RCP<panzer::Basis> > >& dofBasisPair = pb.getProvidedDOFs();
    RCP<panzer::Basis> basis;
    for (vector<pair<string,RCP<panzer::Basis> > >::const_iterator it = 
	   dofBasisPair.begin(); it != dofBasisPair.end(); ++it) {
      if (it->first == res_to_dof->second)
	basis = it->second;
    }
    
    TEST_FOR_EXCEPTION(Teuchos::is_null(basis), std::runtime_error,
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
    
    TEST_FOR_EXCEPTION(!pb.cellData().isSide(), std::logic_error,
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
