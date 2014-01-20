#ifndef PIKE_TEST_MOCK_USER_SOLVER_FACTORY
#define PIKE_TEST_MOCK_USER_SOLVER_FACTORY

#include "Pike_Solver_AbstractFactory.hpp"

namespace pike_test {
  
  /** \brief Mock class to test the user defined solver factory support. */
  class UserSolverFactory : public pike::SolverAbstractFactory {

  public:

    UserSolverFactory(const std::string& mySolverType);

    Teuchos::RCP<pike::Solver> 
    buildSolver(const Teuchos::RCP<Teuchos::ParameterList>& p) const;

  private:

    std::string mySolverType_;

  };

}

#endif
