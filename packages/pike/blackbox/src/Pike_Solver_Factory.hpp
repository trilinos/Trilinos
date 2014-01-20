#ifndef PIKE_SOLVER_FACTORY_HPP
#define PIKE_SOLVER_FACTORY_HPP

#include "Pike_Solver_AbstractFactory.hpp"
#include <vector>

namespace pike {

  class SolverFactory : public pike::SolverAbstractFactory {
    
  public:
    
    //! Regsiter a user defined solver factory.
    void addFactory(const Teuchos::RCP<pike::SolverAbstractFactory>& f);
    
    virtual
    Teuchos::RCP<pike::Solver> 
    buildSolver(const Teuchos::RCP<Teuchos::ParameterList>& p) const;

  private:

    void validateParameterList(const Teuchos::RCP<Teuchos::ParameterList>& p);

    std::vector<Teuchos::RCP<pike::SolverAbstractFactory> > userFactories_;
  };

}

#endif
