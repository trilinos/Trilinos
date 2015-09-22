#ifndef PIKE_RXN_MODEL_EVALUATOR_SINGLE_PHYSICS_EQ2_HPP
#define PIKE_RXN_MODEL_EVALUATOR_SINGLE_PHYSICS_EQ2_HPP

#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Pike_MultiphysicsDistributor.hpp"
#include "Pike_Rxn_ModelEvaluator_SingleEqBase.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"
#include <string>

namespace pike_test {

  class RxnSingleEq2 : public pike_test::RxnSingleEqBase {

  public:

    RxnSingleEq2(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd);

    //@{ BlackBoxModelEvaluator derived methods
    
    std::string name() const override;

    bool supportsParameter(const std::string& pName) const override;
    int getNumberOfParameters() const override;
    std::string getParameterName(const int l) const override;
    int getParameterIndex(const std::string& pName) const override;
    void setParameter(const int l, const Teuchos::ArrayView<const double>& p) override;

    Teuchos::ArrayView<const double> getResponse(const int i) const override;
    int getResponseIndex(const std::string& rName) const override;
    std::string getResponseName(const int i) const override;
    bool supportsResponse(const std::string& rName) const override;
    int getNumberOfResponses() const override;

    //@}

    double evaluateError() override;
    void reset() override;

  private:

    void evaluateF(const double& t, 
		   const std::vector<double>& x, 
		   std::vector<double>& f) override;

    // Parameters
    double p_CA_;
    double p_CC_;
  };

  /** \brief non-member ctor
      \relates RxnSingleEq2
  */
  Teuchos::RCP<pike_test::RxnSingleEq2> 
  rxnSingleEq2(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd);

}

#endif
