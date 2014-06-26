#ifndef PIKE_SOLVER_DEFAULT_BASE_HPP
#define PIKE_SOLVER_DEFAULT_BASE_HPP

#include "Pike_Solver.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"

namespace Teuchos {
  class ParameterList;
}

namespace pike {

  class SolverDefaultBase : public pike::Solver,
			    public Teuchos::ParameterListAcceptorDefaultBase {

  public:

    SolverDefaultBase();

    virtual ~SolverDefaultBase();

    virtual void registerComm(const Teuchos::RCP<const Teuchos::Comm<int> >& comm);

    virtual void registerModelEvaluator(const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me);

    virtual void registerDataTransfer(const Teuchos::RCP<pike::DataTransfer>& dt);

    virtual void completeRegistration();

    virtual Teuchos::RCP<const pike::BlackBoxModelEvaluator> 
    getModelEvaluator(const std::string& name) const;

    virtual const std::vector<Teuchos::RCP<const pike::BlackBoxModelEvaluator> > getModelEvaluators() const;

    virtual Teuchos::RCP<const pike::DataTransfer> 
    getDataTransfer(const std::string& name) const;

    virtual const std::vector<Teuchos::RCP<const pike::DataTransfer> > getDataTransfers() const;

    virtual void stepImplementation() = 0;

    virtual pike::SolveStatus step();

    /** \brief If not called manually by the user, this will be called
	internally on the first call to solve() */
    virtual void initialize();

    virtual pike::SolveStatus solve();

    /** \brief In not called manually by the user, this will be called
	in the destructor of this solver. */
    virtual void finalize();

    virtual void reset();

    virtual pike::SolveStatus getStatus() const; 

    virtual int getNumberOfIterations() const;

    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList);

    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

    Teuchos::RCP<Teuchos::ParameterList> getNonconstValidParameters();

    virtual void addObserver(const Teuchos::RCP<pike::SolverObserver>& observer);

    virtual std::vector<Teuchos::RCP<pike::SolverObserver> > getObservers() const;

    virtual void setStatusTests(const Teuchos::RCP<pike::StatusTest>& statusTests);

    virtual Teuchos::RCP<const pike::StatusTest> getStatusTests() const;

    virtual std::string name() const;

  protected:

    typedef std::vector<Teuchos::RCP<SolverObserver> >::iterator ObserverIterator;
    typedef std::vector<Teuchos::RCP<pike::BlackBoxModelEvaluator> >::iterator ModelIterator;
    typedef std::vector<Teuchos::RCP<pike::BlackBoxModelEvaluator> >::const_iterator ModelConstIterator;
    typedef std::vector<Teuchos::RCP<pike::DataTransfer> >::iterator TransferIterator;
    typedef std::vector<Teuchos::RCP<pike::DataTransfer> >::const_iterator TransferConstIterator;

    int numberOfIterations_;
    Teuchos::RCP<Teuchos::ParameterList> validParameters_;
    Teuchos::RCP<pike::StatusTest> statusTests_;
    std::vector<Teuchos::RCP<SolverObserver> > observers_;
    pike::SolveStatus status_;
    std::vector<Teuchos::RCP<pike::BlackBoxModelEvaluator> > models_;
    std::vector<Teuchos::RCP<pike::DataTransfer> > transfers_;
    bool registrationComplete_;
    std::string name_;

    bool isInitialized_;
    bool isFinalized_;

    // Output
    bool printBeginSolveStatus_;
    bool printStepStatus_;
    bool printEndSolveStatus_;
  };

}

#endif
