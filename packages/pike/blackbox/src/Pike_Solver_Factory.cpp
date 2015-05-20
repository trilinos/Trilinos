#include "Pike_BlackBox_config.hpp"
#include "Pike_Solver_Factory.hpp"
#include "Teuchos_Assert.hpp"

// Pike solvers
#include "Pike_Solver_BlockGaussSeidel.hpp"
#include "Pike_Solver_BlockJacobi.hpp"
#include "Pike_Solver_TransientStepper.hpp"

#include <algorithm>

namespace pike {

  SolverFactory::SolverFactory()
  {
    supportedTypes_.push_back("Block Gauss Seidel");
    supportedTypes_.push_back("Block Jacobi");
  }

  void SolverFactory::addFactory(const Teuchos::RCP<pike::SolverAbstractFactory>& f)
  {
    userFactories_.push_back(f);
  }
  
  bool SolverFactory::supportsType(const std::string& type) const
  {
    if (std::find(supportedTypes_.begin(),supportedTypes_.end(),type) != supportedTypes_.end())
      return true;

    for (std::vector<Teuchos::RCP<pike::SolverAbstractFactory> >::const_iterator i=userFactories_.begin();
	 i != userFactories_.end(); ++i) {
      if ( (*i)->supportsType(type) )
	return true;
    }

    return false;
  }

  Teuchos::RCP<pike::Solver> 
  SolverFactory::buildSolver(const Teuchos::RCP<Teuchos::ParameterList>& p,
			     const std::string& inSolverSublistName) const
  {
    std::string solverSublistName = inSolverSublistName;
    if (solverSublistName == "") {
      TEUCHOS_TEST_FOR_EXCEPTION(!p->isType<std::string>("Solver Sublist Name"), std::logic_error,
				 "The \"Solver Sublist Name\" must be specified but is missing!");
      
      solverSublistName = p->get<std::string>("Solver Sublist Name");
    }

    TEUCHOS_TEST_FOR_EXCEPTION(!p->isSublist(solverSublistName), std::logic_error,
			       "Error the requested \"Solver Sublist Name\" with value \"" 
			       << solverSublistName << "\" is not valid!");
    
    Teuchos::RCP<Teuchos::ParameterList> solverSublist = Teuchos::sublist(p,solverSublistName);
    
    TEUCHOS_TEST_FOR_EXCEPTION(!solverSublist->isType<std::string>("Type"),std::logic_error,
			       "The \"Type\" key must be specified for the solver sublist \"" 
			       << solverSublistName << "\"!");
    
    std::string type = solverSublist->get<std::string>("Type");
    
    Teuchos::RCP<pike::Solver> solver;

    if (type == "Block Gauss Seidel") {
      Teuchos::RCP<pike::BlockGaussSeidel> gs = Teuchos::rcp(new pike::BlockGaussSeidel);
      gs->setParameterList(solverSublist);
      solver = gs;
    }
    else if (type == "Block Jacobi") {
      Teuchos::RCP<pike::BlockJacobi> jacobi = Teuchos::rcp(new pike::BlockJacobi);
      jacobi->setParameterList(solverSublist);
      solver = jacobi;
    }
    else if (type == "Transient Stepper") {
      Teuchos::RCP<pike::TransientStepper> trans = Teuchos::rcp(new pike::TransientStepper);
      trans->setParameterList(solverSublist);
      
      std::string internalSolverSublistName = solverSublist->get<std::string>("Internal Solver Sublist");

      Teuchos::RCP<pike::Solver> internalSolver = this->buildSolver(p,internalSolverSublistName);
      trans->setSolver(internalSolver);

      solver = trans;
    }
    else {
      for (std::vector<Teuchos::RCP<pike::SolverAbstractFactory> >::const_iterator i=userFactories_.begin();
	   i != userFactories_.end(); ++i) {
	if ( (*i)->supportsType(type) )
	  solver = (*i)->buildSolver(p);
      }
    }

    TEUCHOS_TEST_FOR_EXCEPTION(is_null(solver), std::runtime_error,
			       "Error - the requested solver \"Type\" with value \""<< type 
			       << "\" is not supported by the Solver Factory!");

    return solver;
  }

}
