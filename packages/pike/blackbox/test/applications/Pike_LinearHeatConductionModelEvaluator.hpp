#ifndef PIKE_LINEAR_HEAT_CONDUCTION_MODEL_EVALUATOR_HPP
#define PIKE_LINEAR_HEAT_CONDUCTION_MODEL_EVALUATOR_HPP

#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Pike_Response_Scalar.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"
#include <string>

namespace pike {

  /** \brief Simple example of model evaluator for unit testing

      This model simulates linear heat conduction through a 1D slab.
      The solution is analytic so we can use this as a simple test
      problem. 

      We can string this model together multiple times to simulate
      multiphysics coupling through a composite wall via interfacial
      coupling.  This object can operate in two modes:

      1. T_RIGHT_IS_RESPONSE: Given input parameters of thermal
      conductivity, \f$k\f$, heat flux, \f$q\f$, and temperature
      T_left, compute a response temperature T_right at the interface.

      2. Q_IS_RESPONSE: To close the model, right most wall will be
      given the thermal conductivity, temperatures and compute a
      response heat flux that it will pass to all the other
      applications.

      See any book on heat transfer.  For example, "Transport
      Phenomena" by Bird, Stewart and ightfoot, 2nd edition, section
      10.6 "Heat Conduction Through Composite Walls", pages 303-305.
   */
  class LinearHeatConductionModelEvaluator : public pike::BlackBoxModelEvaluator {

  public:

    enum Mode {
      T_RIGHT_IS_RESPONSE,
      Q_IS_RESPONSE
    };

    LinearHeatConductionModelEvaluator(Teuchos::RCP<Teuchos::Comm<int> > comm,
				       std::string name,
				       Mode mode);

    ~LinearHeatConductionModelEvaluator();

    //@{ BlackBoxModelEvaluator derived methods
    
    virtual const std::string name() const;

    void solve();

    bool isConverged() const;

    Teuchos::RCP<pike::Response> getResponse(const int i) const;

    int getResponseIndex(const std::string name) const;

    bool supportsResponse(const std::string name) const;

    //@}

    // Possible paramters to set
    void set_k(const double& k);
    void set_q(const double& q);
    void set_T_left(const double& T_left);
    void set_T_right(const double& T_right);
    
    // Possible responses
    double get_q() const;
    double get_K() const;
    double get_T_left() const;
    double get_T_right() const;
    
  private:
    Teuchos::RCP<Teuchos::Comm<int> > comm_;
    std::string name_;
    Mode mode_;

    double k_;
    double q_;
    double T_left_;
    double T_right_;
    
    std::map<std::string,int> responseMap_;
    std::vector<Teuchos::RCP<ScalarResponse<double> > > responseValue_;
    
  };

  /** \brief non-member ctor
      \relates LinearHeatConductionModelEvaluator
  */
  Teuchos::RCP<pike::LinearHeatConductionModelEvaluator> 
  linearHeatConductionModelEvaluator(Teuchos::RCP<Teuchos::Comm<int> > comm,
				     std::string name,
				     pike::LinearHeatConductionModelEvaluator::Mode mode);
  

}

#endif
