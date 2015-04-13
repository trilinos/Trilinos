#ifndef PIKE_MOCK_MODEL_EVALUATOR_HPP
#define PIKE_MOCK_MODEL_EVALUATOR_HPP

#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"
#include <string>

namespace pike { class Solver; }

namespace pike_test {

  /** \brief Mock model evaluator for unit testing

      Users can set the solver iteration to either converge or fail on
      (and can choose whether it is a local or global convergence
      failure).
   */
  class MockModelEvaluator : public pike::BlackBoxModelEvaluator {

  public:

    enum Mode {
      //! Triggers failure of the local solve
      LOCAL_FAILURE,
      //! Triggers convergence of this model's global coupled problem
      GLOBAL_CONVERGENCE
    };

    MockModelEvaluator(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
		       const std::string& myName,
		       const Mode mode,
		       const int iterationTrigger,
		       const int responseFreezeIteration);
    
    //@{ BlackBoxModelEvaluator derived methods
    
    std::string name() const;

    void solve();

    bool isLocallyConverged() const;

    bool isGloballyConverged() const;

    Teuchos::ArrayView<const double> getResponse(const int i) const;
    int getResponseIndex(const std::string& rName) const;
    std::string getResponseName(const int i) const;
    bool supportsResponse(const std::string& rName) const;
    int getNumberOfResponses() const;

    bool supportsParameter(const std::string& pName) const;
    int getNumberOfParameters() const;
    std::string getParameterName(const int l) const;    
    int getParameterIndex(const std::string& pName) const;
    void setParameter(const int l, const Teuchos::ArrayView<const double>& p);

    bool isTransient() const;
    double getCurrentTime() const;
    double getTentativeTime() const;
    bool solvedTentativeStep() const;
    double getCurrentTimeStepSize() const;
    double getDesiredTimeStepSize() const;    
    double getMaxTimeStepSize() const;
    void setNextTimeStepSize(const double& dt);
    void acceptTimeStep();

    void setSolver(const Teuchos::RCP<pike::Solver>& solver);

    //@}

  private:
    Teuchos::RCP<const Teuchos::Comm<int> > comm_;
    std::string name_;
    Mode mode_;
    int iterationTrigger_;
    int responseFreezeIteration_;
    Teuchos::RCP<pike::Solver> solver_;
    std::map<std::string,int> responseMap_;
    std::vector<std::string> responseNames_;
    std::vector<std::vector<double> > responseValues_; 
    std::map<std::string,int> parameterMap_;
    std::vector<std::string> parameterNames_;
    std::vector<std::vector<double> > parameterValues_; 
  };

  /** \brief non-member ctor
      \relates MockModelEvaluator
  */
  Teuchos::RCP<pike_test::MockModelEvaluator> 
  mockModelEvaluator(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
		     const std::string& name,
		     const pike_test::MockModelEvaluator::Mode mode,
		     const int iterationTrigger,
		     const int responseFreezeIteration);
  

}

#endif
