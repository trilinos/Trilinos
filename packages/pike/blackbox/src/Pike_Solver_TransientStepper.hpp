#ifndef PIKE_SOLVER_TRANSIENT_STEPPER_HPP
#define PIKE_SOLVER_TRANSIENT_STEPPER_HPP

#include "Pike_Solver.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"

namespace pike {

  class TransientStepper : public pike::Solver,
			   public Teuchos::ParameterListAcceptorDefaultBase {
    
  public:

    TransientStepper();
    ~TransientStepper();

    //! Must be called before any other calls to this object.
    void setSolver(const Teuchos::RCP<pike::Solver>& solver);

    // Derived from pike::Solver
    void registerComm(const Teuchos::RCP<const Teuchos::Comm<int> >& comm);
    void registerModelEvaluator(const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me);
    void registerDataTransfer(const Teuchos::RCP<pike::DataTransfer>& dt);
    void completeRegistration();
    Teuchos::RCP<const pike::BlackBoxModelEvaluator> getModelEvaluator(const std::string& name) const;
    const std::vector<Teuchos::RCP<const pike::BlackBoxModelEvaluator> > getModelEvaluators() const;
    Teuchos::RCP<const pike::DataTransfer> getDataTransfer(const std::string& name) const;
    const std::vector<Teuchos::RCP<const pike::DataTransfer> > getDataTransfers() const;
    pike::SolveStatus step();
    void initialize();
    pike::SolveStatus solve();
    void finalize();
    void reset();
    pike::SolveStatus getStatus() const;
    int getNumberOfIterations() const;
    void addObserver(const Teuchos::RCP<pike::SolverObserver>& observer);
    std::vector<Teuchos::RCP<pike::SolverObserver> > getObservers() const;
    void setStatusTests(const Teuchos::RCP<pike::StatusTest>& statusTests);
    Teuchos::RCP<const pike::StatusTest> getStatusTests() const;
    std::string name() const;

    // Derived from ParameterListAcceptorDefaultBase
    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList);
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
    Teuchos::RCP<Teuchos::ParameterList> getNonconstValidParameters();

  private:

    int currentTimeStep_;
    int maxTimeSteps_;
    double beginTime_;
    double endTime_;
    double currentTime_;
    double initialStepSize_;
    double currentStepSize_;
    double minStepSize_;
    double maxStepSize_;
    double stepGrowthFactor_;
    double stepDecreaseFactor_;
    Teuchos::Array<double> checkPointsVec_;
    std::list<double> checkPoints_;
    bool printTimeStepSummary_;
    bool printTimeStepDetails_;
    int numConvergedTimeStepsBeforeGrowth_;

    pike::SolveStatus overallStatus_;
    pike::SolveStatus timeStepStatus_;
    std::vector<Teuchos::RCP<pike::BlackBoxModelEvaluator> > transientModels_;

    int totalNumFailedSteps_;
    int numConsecutiveFailedTimeSteps_;
    int numConsecutiveConvergedTimeSteps_;

    Teuchos::RCP<pike::Solver> solver_;
    Teuchos::RCP<Teuchos::ParameterList> validParameters_;
    bool registrationComplete_;
  };

}

#endif
