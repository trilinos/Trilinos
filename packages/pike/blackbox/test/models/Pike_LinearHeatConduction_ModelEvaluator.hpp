#ifndef PIKE_LINEAR_HEAT_CONDUCTION_MODEL_EVALUATOR_HPP
#define PIKE_LINEAR_HEAT_CONDUCTION_MODEL_EVALUATOR_HPP

#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"
#include <string>

namespace pike_test {

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
      Phenomena" by Bird, Stewart and Lightfoot, 2nd edition, section
      10.6 "Heat Conduction Through Composite Walls", pages 303-305.
   */
  class LinearHeatConductionModelEvaluator : public pike::BlackBoxModelEvaluator {

  public:

    enum Mode {
      T_RIGHT_IS_RESPONSE,
      Q_IS_RESPONSE
    };

    LinearHeatConductionModelEvaluator(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
				       const std::string& myName,
				       const Mode mode);

    //@{ BlackBoxModelEvaluator derived methods
    
    std::string name() const;
    void solve();
    bool isLocallyConverged() const;
    bool isGloballyConverged() const;

    virtual bool supportsParameter(const std::string& pName) const;
    virtual int getNumberOfParameters() const;
    virtual std::string getParameterName(const int l) const;
    virtual int getParameterIndex(const std::string& pName) const;
    virtual void setParameter(const int l, const Teuchos::ArrayView<const double>& p);

    Teuchos::ArrayView<const double> getResponse(const int i) const;
    int getResponseIndex(const std::string& rName) const;
    std::string getResponseName(const int i) const;
    bool supportsResponse(const std::string& rName) const;
    int getNumberOfResponses() const;

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
    Teuchos::RCP<const Teuchos::Comm<int> > comm_;
    std::string name_;
    Mode mode_;

    double k_;
    double q_;
    double T_left_;
    double T_right_;
    
    std::map<std::string,int> parameterMap_;
    std::vector<std::string> parameterNames_;

    std::map<std::string,int> responseMap_;
    std::vector<std::string> responseNames_;
    std::vector<std::vector<double> > responseValues_;
    
  };

  /** \brief non-member ctor
      \relates LinearHeatConductionModelEvaluator
  */
  Teuchos::RCP<pike_test::LinearHeatConductionModelEvaluator> 
  linearHeatConductionModelEvaluator(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
				     const std::string& name,
				     const pike_test::LinearHeatConductionModelEvaluator::Mode mode);
  

}

#endif
