#ifndef PIKE_RXN_MODEL_EVALUATOR_ALL_HPP
#define PIKE_RXN_MODEL_EVALUATOR_ALL_HPP

#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Pike_MultiphysicsDistributor.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"
#include <string>

namespace pike_test {

  class RxnAll : public pike::BlackBoxModelEvaluator {

  public:

    RxnAll(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd);

    //@{ BlackBoxModelEvaluator derived methods
    
    std::string name() const;
    void solve();
    bool isLocallyConverged() const;

    bool supportsParameter(const std::string& pName) const;
    int getNumberOfParameters() const;
    std::string getParameterName(const int l) const;
    int getParameterIndex(const std::string& pName) const;
    void setParameter(const int l, const Teuchos::ArrayView<const double>& p);

    Teuchos::ArrayView<const double> getResponse(const int i) const;
    int getResponseIndex(const std::string& rName) const;
    std::string getResponseName(const int i) const;
    bool supportsResponse(const std::string& rName) const;
    int getNumberOfResponses() const;

    bool isTransient() const;
    double getCurrentTime() const;
    double getTentativeTime() const;
    bool solvedTentativeStep() const;
    double getCurrentTimeStepSize() const;
    double getDesiredTimeStepSize() const;
    double getMaxTimeStepSize() const;
    void setNextTimeStepSize(const double& dt);
    void acceptTimeStep();
    //@}

    double evaluateError();
    void reset();

  private:

    void evaluateF(const double& t, 
		   const std::vector<double>& x, 
		   std::vector<double>& f);

    Teuchos::RCP<pike::MultiphysicsDistributor> mpd_;

    // Solution (and tentative solution if x not accepted yet)
    std::vector<double> x_;
    // Old solution at n-1
    std::vector<double> xOld_;
    // Analytic Solution
    std::vector<double> xAnalytic_;
    // RK tmp vectors
    std::vector<double> Y_;
    std::vector<double> f_;

    // Parameters
    double p_k1_;
    double p_k2_;
    double p_CA0_;
    double p_CB0_;
    double p_CC0_;

    // Responses

    // Transient members
    int timeStepNumber_;
    bool solvedTentativeStep_;
    double currentTime_;
    double tentativeTime_;
    double currentTimeStepSize_;
    bool isLocallyConverged_;
  };

  /** \brief non-member ctor
      \relates RxnAll
  */
  Teuchos::RCP<pike_test::RxnAll> 
  rxnAll(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd);

}

#endif
