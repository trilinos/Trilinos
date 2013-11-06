//@HEADER
// ***********************************************************************
//
//                     Rapid Optimization Library
//
// Questions? Contact:    Drew Kouri (dpkouri@sandia.gov)
//                      Denis Ridzal (dridzal@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef ROL_LINESEARCHSTEP_H
#define ROL_LINESEARCHSTEP_H

#include "ROL_Step.hpp"
#include "ROL_Secant.hpp"
#include <sstream>
#include <iomanip>

/** \class ROL::LineSearchStep
    \brief Provides the interface to compute optimization steps
           with line search.
*/


namespace ROL {

template <class Real>
class LineSearchStep : public Step<Real> {
private:

  bool useSecant_;
  Teuchos::RCP<Secant<Real> > secant_;
  
  int it_;
  int maxit_;
  Real rho_;
  Real alpha0_;
  Real c_;

  std::string step_;

public:

  virtual ~LineSearchStep() {}

  LineSearchStep(int maxit = 20, Real rho = 0.5, Real alpha0 = 1.0, Real c = 1.e-4, 
                 bool useSecant = true, SecantType type = Secant_lBFGS, int L = 10, int BBtype = 1 ) 
    : useSecant_(useSecant), maxit_(maxit), rho_(rho), alpha0_(alpha0), c_(c) {
     
    if ( useSecant_ ) {
      if ( type == Secant_lBFGS ) {
        secant_ = Teuchos::rcp( new lBFGS<Real>(L) );
        step_   = "Limited-Memory BFGS";
      }
      else if ( type == Secant_lDFP ) {
        secant_ = Teuchos::rcp( new lDFP<Real>(L) );
        step_   = "Limited-Memory DFP";
      }
      else if ( type == Secant_BarzilaiBorwein ) {
        secant_ = Teuchos::rcp( new BarzilaiBorwein<Real>(BBtype) );
        step_   = "Barzilai-Borwein";
        //maxit_  = 0;
      } 
    }
    else {
      step_ = "Gradient Descent";
    }
  }

  /** \brief Compute step.
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj, AlgorithmState<Real> &algo_state ) {
    if ( useSecant_ ) {
      secant_->applyH(s,*(Step<Real>::state_->gradientVec),x,algo_state.iter);
      //secant_->test(s,x,algo_state.iter);
    }
    else {
      s.set(*(Step<Real>::state_->gradientVec));
    }
    s.scale(-1.0);
    Real gs = (Step<Real>::state_->gradientVec)->dot(s);

    // Perform backtracking line search from Nocedal/Wright.
    it_        = 0;
    Real alpha = alpha0_;
    Real f     = algo_state.value;

    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    xnew->set(x);
    xnew->axpy(alpha, s);
    Real fnew = obj.value(*xnew);

    while ( (fnew > f + c_*alpha*gs) && (it_ < maxit_) ) {
      alpha *= rho_;
      xnew->set(x);
      xnew->axpy(alpha, s);
      fnew = obj.value(*xnew);
      it_++;
    }
    s.scale(alpha);

    // Update step state information
    (Step<Real>::state_->descentVec)->set(s);

    // Update algorithm state information
    algo_state.snorm = s.norm();
    algo_state.value = fnew;
  }

  /** \brief Update step, if successful.
  */
  void update( Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj, AlgorithmState<Real> &algo_state ) {
    // Update iterate
    x.axpy(1.0, s);

    // Compute new gradient
    Teuchos::RCP<Vector<Real> > gp;
    if ( useSecant_ ) {
      gp = x.clone();
      gp->set(*(Step<Real>::state_->gradientVec));
    }
    obj.gradient(*(Step<Real>::state_->gradientVec),x);

    // Update Secant Information
    if ( useSecant_ ) {
      secant_->update(*(Step<Real>::state_->gradientVec),*gp,s,algo_state.snorm);
    }

    // Update algorithm state
    (algo_state.iterateVec)->set(x);
    algo_state.gnorm = (Step<Real>::state_->gradientVec)->norm();
    algo_state.iter++;
  }

  /** \brief Print iterate status.
  */
  std::string print( AlgorithmState<Real> & algo_state, bool printHeader = false ) const  {
    std::stringstream hist;
    if ( algo_state.iter == 0 ) {
      hist << "\n" << step_ << " with Backtracking Linesearch\n"; 
      hist << "  ";
      hist << std::setw(15) << std::left << "iter";  
      hist << std::setw(15) << std::left << "value";
      hist << std::setw(15) << std::left << "gnorm"; 
      hist << std::setw(15) << std::left << "snorm";
      hist << std::setw(10) << std::left << "ls";
      hist << "\n";
      hist << std::setfill('-') << std::setw(80) << "-" << "\n";
    }
    else {
      hist << "  "; 
      hist << std::setw(15) << std::left << algo_state.iter;  
      hist << std::setw(15) << std::left << algo_state.value; 
      hist << std::setw(15) << std::left << algo_state.gnorm; 
      hist << std::setw(15) << std::left << algo_state.snorm; 
      hist << std::setw(10) << std::left << it_;              
      hist << "\n";
    }
    return hist.str();
  }

  // struct StepState (scalars, vectors) map?

  // getState

  // setState

}; // class Step

} // namespace ROL

#endif
