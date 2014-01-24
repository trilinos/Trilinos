#ifndef PIKE_SOLVER_ABSTRACT_FACTORY_HPP
#define PIKE_SOLVER_ABSTRACT_FACTORY_HPP

#include "Teuchos_RCP.hpp"

namespace Teuchos { class ParameterList; }

namespace pike {

  class Solver;

  /** \brief Abstract factory design pattern for building solvers. */
  class SolverAbstractFactory {

  public:

    virtual ~SolverAbstractFactory() {}

    //! Returns true of the type can be built by this factory.
    virtual bool supportsType(const std::string& type) const = 0;

    /** \brief Allocates and returns a pike::Solver object.
	
	\param[in] p Parameter list the determines the solver to build.
	\returns A newly allocated pike::Solver 
    */
    virtual 
    Teuchos::RCP<pike::Solver> 
    buildSolver(const Teuchos::RCP<Teuchos::ParameterList>& p) const = 0;

  };

}

#endif
