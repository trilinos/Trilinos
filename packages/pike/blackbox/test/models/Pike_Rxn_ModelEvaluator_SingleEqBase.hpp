#ifndef PIKE_RXN_MODEL_EVALUATOR_SINGLE_PHYSICS_BASE_HPP
#define PIKE_RXN_MODEL_EVALUATOR_SINGLE_PHYSICS_BASE_HPP

#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Pike_MultiphysicsDistributor.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"
#include <string>

namespace pike_test {

  class RxnSingleEqBase : public pike::BlackBoxModelEvaluator {

  public:

    RxnSingleEqBase(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd);

    //@{ BlackBoxModelEvaluator derived methods
    
    virtual std::string name() const = 0;
    void solve();
    bool isLocallyConverged() const;

    virtual bool supportsParameter(const std::string& pName) const = 0;
    virtual int getNumberOfParameters() const = 0;
    virtual std::string getParameterName(const int l) const = 0;
    virtual int getParameterIndex(const std::string& pName) const = 0;
    virtual void setParameter(const int l, const Teuchos::ArrayView<const double>& p) = 0;

    virtual Teuchos::ArrayView<const double> getResponse(const int i) const = 0;
    virtual int getResponseIndex(const std::string& rName) const = 0;
    virtual std::string getResponseName(const int i) const = 0;
    virtual bool supportsResponse(const std::string& rName) const = 0;
    virtual int getNumberOfResponses() const = 0;

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

    virtual double evaluateError() = 0;
    virtual void reset() = 0;

  protected:

    virtual void evaluateF(const double& t, 
			   const std::vector<double>& x, 
			   std::vector<double>& f) = 0;

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

}

#endif
