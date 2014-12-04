#ifndef PIKE_SOLVER_FACTORY_HPP
#define PIKE_SOLVER_FACTORY_HPP

#include "Pike_Solver_AbstractFactory.hpp"
#include <vector>

namespace pike {

  class SolverFactory : public pike::SolverAbstractFactory {
    
  public:
    
    SolverFactory();

    //! Regsiter a user defined solver factory.
    void addFactory(const Teuchos::RCP<pike::SolverAbstractFactory>& f);
    
    bool supportsType(const std::string& type) const;

    virtual
    Teuchos::RCP<pike::Solver> 
    buildSolver(const Teuchos::RCP<Teuchos::ParameterList>& p,
		const std::string& solverSublistName = "") const;

  private:

    void validateParameterList(const Teuchos::RCP<Teuchos::ParameterList>& p);

    std::vector<Teuchos::RCP<pike::SolverAbstractFactory> > userFactories_;

    std::vector<std::string> supportedTypes_;
  };

}

#endif
