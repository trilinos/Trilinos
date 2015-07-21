#ifndef PIKE_SOLVER_HPP
#define PIKE_SOLVER_HPP

#include "Pike_BlackBox_config.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Pike_StatusTest.hpp"
#include <vector>

namespace Teuchos {template<typename> class Comm;}

namespace pike {

  class BlackBoxModelEvaluator;
  class DataTransfer;
  class SolverObserver;


  /** \brief Pure virtual base class (strategy design pattern) for iterative solvers */
  class Solver : public Teuchos::Describable,
		 public Teuchos::VerboseObject<pike::Solver> {

  public:

    virtual ~Solver() {}

    /** \brief Optionally register a comm.

	Some solvers provide the option to globally barrier between solver actions.  For such solvers, set the comm with this method.   
    */
    virtual void registerComm(const Teuchos::RCP<const Teuchos::Comm<int> >& comm) = 0;

    /** \brief Register an application/model evaluator with the solver. 

	Only allowed to be called before the call to completeRegistrion().
    */
    virtual void registerModelEvaluator(const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me) = 0;

    /** \brief Register a DataTransfer with the solver. 

	Only allowed to be called before the call to completeRegistration().
    */
    virtual void registerDataTransfer(const Teuchos::RCP<pike::DataTransfer>& dt) = 0;

    /** \brief Finalizes the solver so that it can begin solving the problem.
	
	Once this is called, the methods registerModelEvaluator() and
	registerDataTransfer() can no longer be called.
     */
    virtual void completeRegistration() = 0;

    //! Returns the requested model evaluator.
    virtual Teuchos::RCP<const pike::BlackBoxModelEvaluator> getModelEvaluator(const std::string& name) const = 0;

    //! Returns all registered model evaluators.
    virtual const std::vector<Teuchos::RCP<const pike::BlackBoxModelEvaluator> > getModelEvaluators() const = 0;

    //! Return the requested data transfer.
    virtual Teuchos::RCP<const pike::DataTransfer> getDataTransfer(const std::string& name) const = 0;

    //! Return all registered data transfers.
    virtual const std::vector<Teuchos::RCP<const pike::DataTransfer> > getDataTransfers() const = 0;

    /** \brief Take one step of the solve iteration sequence.
    
         @return Current SolveStatus.  
    */
    virtual pike::SolveStatus step() = 0;

    /** \brief Initialize the system. This is typically used to
	trigger one time initialization requirements in the
	pike::SolverObserver::observeInitialization(). */
    virtual void initialize() = 0;

    //! Solve the system (step() until either converged or failed). Returns the current SolveStatus.
    virtual pike::SolveStatus solve() = 0;

    /** \brief Prepare to terminate the simulation. This is typically
	used to trigger one time stopping requirements in the
	pike::SolverObserver::observeFianlize(). */
    virtual void finalize() = 0;

    //! Reset the solver to reuse for another solve.
    virtual void reset() = 0;

    //! Returns the current SolveStatus.
    virtual pike::SolveStatus getStatus() const = 0;

    //! Returns the current number of iterations.
    virtual int getNumberOfIterations() const = 0;

    //! Register an observer with the solver.
    virtual void addObserver(const Teuchos::RCP<pike::SolverObserver>& observer) = 0;

    //! Returns the observers registered with the solver.
    virtual std::vector<Teuchos::RCP<pike::SolverObserver> > getObservers() const = 0;

    /** \brief Sets/overrides the status test to use with the solver.

	Users can optionally set their own status tests via this
	method.  If no status test is set, then the solver will build
	the status tests based on input parameters and/or defaults.

	@param[in] statusTests The status test used to determine when to stop the solver.
     */ 
    virtual void setStatusTests(const Teuchos::RCP<pike::StatusTest>& statusTests) = 0;

    //! Returns the status tests.
    virtual Teuchos::RCP<const pike::StatusTest> getStatusTests() const = 0;
    
    /** \brief Returns the name of the solver.  This is set by the
        input parameter list.  This can be used to debug hierarchical
        solves.
    */
    virtual std::string name() const = 0;
  };

}

#endif
