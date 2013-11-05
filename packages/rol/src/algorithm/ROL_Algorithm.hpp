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
#include "ROL_StatusTest.hpp"

/** \class ROL::Algorithm
    \brief Provides an interface to run optimization algorithms.
*/


namespace ROL {

template <class Real>
class DefaultAlgorithm {
private:
  Teuchos::RCP<Step<Real> >           step_;
  Teuchos::RCP<StatusTest<Real> >     status_;
  Teuchos::RCP<AlgorithmState<Real> > state_;

public:

  virtual ~DefaultAlgorithm() {}

  DefaultAlgorithm(Step<Real> & step, StatusTest<Real> & status) {
    step_   = Teuchos::rcp(&step,   false);
    status_ = Teuchos::rcp(&status, false);
    state_  = Teuchos::rcp(new AlgorithmState<Real>);
    state_->iter = 0;
  }

  /** \brief Run algorithm.
  */
  virtual void run( Vector<Real> &x, Objective<Real> &obj ) {
    // Initialize Current Iterate Container 
    if ( state_->iterateVec == Teuchos::null ) {
      state_->iterateVec = x.clone();
      state_->iterateVec->set(x);
    }

    // Initialize Step Container
    Teuchos::RCP<Vector<Real> > s = x.clone();

    // Initialize Step
    step_->initialize(x, obj, *state_);
    std::cout << step_->print(*state_,true);

    // Run Algorithm
    while (status_->check(*state_)) {
      step_->compute(*s, x, obj, *state_);
      step_->update(x, *s, obj, *state_);
      std::cout << step_->print(*state_);
    }
  }

}; // class DefaultAlgorithm

} // namespace ROL

#endif
