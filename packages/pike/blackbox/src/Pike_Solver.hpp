#ifndef PIKE_SOLVER_HPP
#define PIKE_SOLVER_HPP

#include "Pike_BlackBox_config.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Pike_StatusTest.hpp"
#include <vector>

namespace pike {

  class BlackBoxModelEvaluator;
  class DataTransfer;
  class Observer;


  /** \brief Pure virtual base class for iterative solvers */
  class Solver : public Teuchos::ParameterListAcceptorDefaultBase,
		 public Teuchos::Describable,
		 public Teuchos::VerboseObject<pike::Solver> {

  public:

    virtual ~Solver() {}

    /** \brief Register an application/model evaluator with the solver. 

	Only allowed to be called before the call to completeRegistrion()
    */
    virtual void registerModelEvaluator(const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me) = 0;

    /** \brief Register a DataTransfer with the solver. 

	Only allowed to be called before the call to completeRegistrion()
    */
    virtual void registerDataTransfer(const Teuchos::RCP<pike::DataTransfer>& dt) = 0;

    /** \brief Finalizes the solver so that it can begin solving the problem.
	
	Once this is called, the methods registerModelEvaluator() and
	registerDataTransfer() can no longer be called.
     */
    virtual void completeRegistration() = 0;

    /** \brief Take one step of the solve iteration sequence. 
    
         @return Current SolveStatus.  
    */
    virtual pike::SolveStatus step() = 0;

    //! Solve the system (step() until either converged or failed). Returns the current SolveStatus.
    virtual pike::SolveStatus solve() = 0;

    //! Reset the solver to reuse for another solve.
    virtual void reset() = 0;

    //! Returns the current SolveStatus.
    virtual pike::SolveStatus getStatus() const = 0;

    //! Returns the current number of iterations.
    virtual int getNumberOfIterations() const = 0;

    virtual void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList) = 0;

    virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const = 0;

    //! Register an observer with the solver.
    virtual void addObserver(const Teuchos::RCP<pike::Observer>& observer) = 0;

    //! Returns the observers registered with the solver.
    virtual std::vector<Teuchos::RCP<pike::Observer> > getObservers() const = 0;

    /** \brief Sets the status test to use with the solver.

	Users can optionally set their own status tests via this
	method.  If no status test is set, then the solver will build
	the status tests based either on input parameters or on the
	default.

	@param[in] statusTests The status test used to determine when to stop the solver.
     */ 
    virtual void setStatusTests(const Teuchos::RCP<pike::StatusTest>& statusTests) = 0;

    //! Returns the status tests.
    virtual Teuchos::RCP<const pike::StatusTest> getStatusTests() const = 0;
    
  };

  std::ostream& operator<<(std::ostream& os, const pike::Solver& solver);
}

#endif
