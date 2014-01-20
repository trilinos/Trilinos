#include "Pike_Mock_UserSolverFactory.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Pike_Solver_BlockJacobi.hpp"

namespace pike_test {

  UserSolverFactory::UserSolverFactory(const std::string& mySolverType) :
    mySolverType_(mySolverType)
  { }

  Teuchos::RCP<pike::Solver> 
  UserSolverFactory::buildSolver(const Teuchos::RCP<Teuchos::ParameterList>& p) const
  {
    // for real classes there should be safety checks on the
    // parameter list access.  Since this is a mock object for unit
    // testing, not going to bother writing the safety checks.
    
    std::string solverSublistName = p->get<std::string>("Solver Sublist Name");
    std::string type = p->sublist(solverSublistName).get<std::string>("Type");
    
    Teuchos::RCP<pike::Solver> solver;
    
    if (type == mySolverType_) {
      // instead of implementing a mock solver for testing we just
      // reuse a block Jacobi solver.  Need to change the "Type" to
      // "Block Jacobi" so that the parmeter list validation passes.
      p->set("Type","Block Jacobi");
      solver = Teuchos::rcp(new pike::BlockJacobi);
      solver->setParameterList(Teuchos::sublist(p,solverSublistName));
    }
    
    return solver;
  }

}
