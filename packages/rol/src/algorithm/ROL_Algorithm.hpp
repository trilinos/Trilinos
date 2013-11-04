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

template <class Real>
class DefaultAlgorithm {
private:
  Teuchos::RCP<Step<Real> > step_;
  //StatusTest

public:

  virtual ~DefaultAlgorithm() {}

  DefaultAlgorithm(Step<Real> & step) {
    step_ = Teuchos::rcp(&step, false);
  }

  /** \brief Run algorithm.
  */
  virtual void run( Vector<Real> &x, Objective<Real> &obj ) {
    Teuchos::RCP<Vector<Real> > s = x.clone();
    for (int i=0; i<1000; i++) {
      step_->compute(*s, x, obj);
      step_->update(x, *s, obj);
    }
  }

}; // class DefaultAlgorithm

} // namespace ROL

#endif
