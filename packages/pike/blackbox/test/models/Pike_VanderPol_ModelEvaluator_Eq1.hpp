#ifndef PIKE_VANDERPOL_MODEL_EVALUATOR_1_HPP
#define PIKE_VANDERPOL_MODEL_EVALUATOR_1_HPP

#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Pike_MultiphysicsDistributor.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"
#include <string>

namespace pike_test {

  /** \brief Simple example of Van der Pol model evaluator for unit testing

   * This is the Van der Pol Equation from 
   * "Solving Differential Equations I, Nonstiff Problems" by
   * E. Hairer, S.P. Norsett, and G. Wanner.  2000, Second edition, Springer.
   * 
   * \ddot{x_1} = x_2
   * \ddot{x_2} = \epsilon*(1-x_1^2)*x_2-x_1, \epsilon > 0
   *
   * Hairer, etal. show the solution x(t) for large \epsilon:
   *
   * log(x_1) - \frac{x_1^2}{2} = \frac{t-t_0}{\epsilon} + C
   *
   * There is one parameter in this model, \epsilon.
   *
   */
// Nonlinear ODE system with manufactured solution based on asymptotic solution
// for small epsilon
//
// f[0] = x_dot[0] - x[1];
// f[1] = x_dot[1] - (eps*(1.0-x[0]*x[0])*x[1]-x[0]) - forcing_term;
//
// forcing term is defined so that exact solution is given by
//
// x_0(t) =  2*cos(t)+eps*(0.75*sin(t)-0.25*sin(3.0*t))
// x_1(t) = -2*sin(t)+eps*(0.75*cos(t)-0.75*cos(3.0*t))
//
// initial conditions for exact solution
//
// x_0(0) = 2
// x_1(0) = 0
//
// exact time derivatives
//
// x_dot_0(t) = -2*sin(t)+eps*(0.75*cos(t)-0.75*cos(3.0*t))
// x_dot_1(t) = -2*cos(t)+eps*(-0.75*sin(t)+3.0*0.75*sin(3.0*t))
//

  class VanderPolME1 : public pike::BlackBoxModelEvaluator {

  public:

    VanderPolME1(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd);

    //@{ BlackBoxModelEvaluator derived methods
    
    std::string name() const;
    void solve();
    bool isLocallyConverged() const;

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

    virtual bool isTransient() const;
    virtual double getCurrentTime() const;
    virtual double getTentativeTime() const;
    virtual bool solvedTentativeStep() const;
    virtual double getCurrentTimeStepSize() const;
    virtual double getDesiredTimeStepSize() const;
    virtual double getMaxTimeStepSize() const;
    virtual void setNextTimeStepSize(const double& dt);
    virtual void acceptTimeStep();
    //@}

    void setUnitTestFailureForTimeStep(const int& timeStepNumber);

  private:
    Teuchos::RCP<pike::MultiphysicsDistributor> mpd_;
    pike::MultiphysicsDistributor::ApplicationIndex eq1_;

    // My internal solution (also the only response)
    double x1_;

    // Parameters
    double p_epsilon_;
    double p_x2_;

    // Transient members
    int timeStepNumber_;
    bool sovledTentativeStep_;
    double currentTime_;
    double tentativeTime_;
    double currentTimeStepSize_;
    bool isLocallyConverged_;
    
    // Model evaluator reports a unit test solve failure on the
    // specified iteration.  This is for UNIT TESTING ONLY!
    int unitTestSolveFailureOnTimeStepNumber_;
    bool unitTestFailureTriggered_;
  };

  /** \brief non-member ctor
      \relates VanderPolME1
  */
  Teuchos::RCP<pike_test::VanderPolME1> 
  vanderPolME1(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd);

}

#endif
