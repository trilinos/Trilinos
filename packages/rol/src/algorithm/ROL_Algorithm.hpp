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

#ifndef ROL_ALGORITHM_H
#define ROL_ALGORITHM_H

#include "ROL_Step.hpp"

/** \class ROL::Algorithm
    \brief Provides an interface to run optimization algorithms.
*/


namespace ROL {


template<class Real>
struct AlgorithmState {
  Real value;
  Teuchos::RCP<Vector<Real> > gradientVec;
  Real gnorm;
  Teuchos::RCP<Vector<Real> > descentVec;
  Real dnorm;
  Teuchos::RCP<Vector<Real> > iterateVec;
  Real inorm;
  int iter;
};


template <class Real>
class DefaultAlgorithm {
private:
  Teuchos::RCP<Step<Real> > step_;
  Teuchos::RCP<StatusTest<Real> > status_;
  Teuchos::RCP<AlgorithmState<Real> > state_;

public:

  virtual ~DefaultAlgorithm() {}

  DefaultAlgorithm(Step<Real> & step, StatusTest<Real> & status) {
    step_   = Teuchos::rcp(&step, false);
    status_ = Teuchos::rcp(&status, false);
    state_  = Teuchos::rcp(new AlgorithmState<Real>);
    state_->iter = 0;
  }

  /** \brief Run algorithm.
  */
  virtual void run( Vector<Real> &x, Objective<Real> &obj ) {
    Teuchos::RCP<Vector<Real> > s = x.clone();
    state_->gradientVec = x.clone();
    state_->descentVec = x.clone();
    state_->iterateVec = x.clone();
    state_->iterateVec->set(x);
    
    while (status_->check(state_)) {
      step_->compute(*s, x, obj);
      step_->update(x, *s, obj);
    }
  }

}; // class DefaultAlgorithm

} // namespace ROL

#endif
