#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Panzer_PhysicsBlock.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Basis.hpp"

// Evaluators
// #include "Panzer_DOF.hpp"
// #include "Panzer_DOFGradient.hpp"
// #include "Panzer_GatherSolution.hpp"
// #include "Panzer_ScatterDirichletResidual.hpp"
// #include "Panzer_Dirichlet_Constant.hpp"

// ***********************************************************************
template <typename EvalT>
user_app::BCStrategy_Dirichlet_Constant<EvalT>::
BCStrategy_Dirichlet_Constant(const panzer::BC& bc) :
  panzer::BCStrategy<EvalT>(bc)
{ }

// ***********************************************************************
template <typename EvalT>
void user_app::BCStrategy_Dirichlet_Constant<EvalT>::
buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
			   const panzer::PhysicsBlock& pb) const
{
  /*
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::vector;
  using std::map;
  using std::string;
  using std::pair;

  // ********************
  // DOFs (unknowns)
  // ********************

  RCP<vector<string> > field_names = rcp(new vector<string>);
  field_names->push_back(this->m_bc.equationSetName());

  const vector<pair<string,RCP<panzer::Basis> > >& dofs = pb.getProvidedDOFs();

  RCP<panzer::Basis> basis;
  for (vector<pair<string,RCP<panzer::Basis> > >::const_iterator dof_it = 
	 dofs.begin(); dof_it != dofs.end(); ++dof_it) {
    if (dof_it->first == this->m_bc.equationSetName())
      basis = dof_it->second;
  }

  TEST_FOR_EXCEPTION(Teuchos::is_null(basis), std::runtime_error,
		     "Error the name \"" << this->m_bc.equationSetName()
		     << "\" is not a valid DOF for the boundary condition:\n"
		     << this->m_bc << "\n");

  // Gather
  {
    ParameterList p("BC Gather");
    p.set("Basis", basis);
    p.set("DOF Names", field_names);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::rf::GatherSolution<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
  }
  
  string bc_name = this->m_bc.identifier();

  // Scatter
  RCP<map<string,string> > names_map = rcp(new map<string,string>);
  RCP<vector<string> > residual_names = rcp(new vector<string>);
  for (vector<string>::iterator i=field_names->begin();
       i != field_names->end(); ++i) {
    residual_names->push_back(bc_name + *i);
    names_map->insert(std::make_pair(bc_name + *i,*i));
  }

  string scatter_field_name = "Dummy Scatter: " + bc_name;
  
  {
    ParameterList p("Scatter");
    p.set("Scatter Name", scatter_field_name);
    p.set("Basis", basis);
    p.set("Evaluated Names", residual_names);
    p.set("Evaluated Map", names_map);
    
    TEST_FOR_EXCEPTION(!pb.cellData().isSide(), std::logic_error,
		       "Error - physics block is not a side set!");
    
    p.set<int>("Side Subcell Dimension", 
	       pb.getBaseCellTopology().getDimension() - 1);
    p.set<int>("Local Side ID", pb.cellData().side());



    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::ScatterDirichletResidual<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
  }

  // Evaluator for Constant dirichlet BCs
  {
    ParameterList p("BC Constant Dirichlet");
    p.set("Residual Name", (*residual_names)[0]);
    p.set("DOF Name", (*names_map)[(*residual_names)[0]]);
    p.set("Data Layout", basis->functional);
    p.set("Value", this->m_bc.constantValue());
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::DirichletConstant<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
  }

  // Require variables
  {

    PHX::Tag<typename EvalT::ScalarT> tag(scatter_field_name, 
					  rcp(new PHX::MDALayout<Dummy>(0)));
    fm.template requireField<EvalT>(tag);
  }
  */
}
