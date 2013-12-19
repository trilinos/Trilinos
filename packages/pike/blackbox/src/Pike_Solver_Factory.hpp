#ifndef PIKE_SOLVER_FACTORY_HPP
#define PIKE_SOLVER_FACTORY_HPP

#include "Teuchos_RCP.hpp"

namespace Teuchos {
  class ParameterList;
}

namespace pike {

  class StatusTest;
  
  /** \brief Allocates a new solver based on parameter list
      
      \param[in] paramList (Required) Parameter list the determines the solver to build.
      \param[in] statusTests (Optional) Status tests for stopping criteria.  If this value is null, the status tests will be built from the parameter list. 
  */
  Teuchos::RCP<pike::Solver> 
  buildSolver(const Teuchos::RCP<Teuchos::ParameterList>& paramList,
	      const Teuchos::RCP<pike::StatusTest>& statusTests = Teuchos::Null);

}

#endif
