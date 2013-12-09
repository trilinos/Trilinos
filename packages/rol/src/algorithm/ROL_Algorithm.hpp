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

  bool printHeader_;

public:

  virtual ~DefaultAlgorithm() {}

  DefaultAlgorithm(Step<Real> & step, StatusTest<Real> & status, bool printHeader = false ) {
    step_   = Teuchos::rcp(&step,   false);
    status_ = Teuchos::rcp(&status, false);
    state_  = Teuchos::rcp(new AlgorithmState<Real>);
    state_->iter = 0;
    printHeader_ = printHeader;
  }

  /** \brief Run algorithm.
  */
  virtual std::vector<std::string> run( Vector<Real> &x, Objective<Real> &obj, bool print = false ) {
    std::vector<std::string> output;

    // Initialize Current Iterate Container 
    if ( this->state_->iterateVec == Teuchos::null ) {
      this->state_->iterateVec = x.clone();
      this->state_->iterateVec->set(x);
    }

    // Initialize Step Container
    Teuchos::RCP<Vector<Real> > s = x.clone();

    // Initialize Step
    this->step_->initialize(x, obj, *this->state_);
    output.push_back(this->step_->print(*this->state_,true));
    if ( print ) {
      std::cout << this->step_->print(*this->state_,true);
    }

    // Run Algorithm
    while (this->status_->check(*state_)) {
      this->step_->compute(*s, x, obj, *this->state_);
      this->step_->update(x, *s, obj, *this->state_);
      output.push_back(this->step_->print(*this->state_,this->printHeader_));
      if ( print ) {
        std::cout << this->step_->print(*this->state_,this->printHeader_);
      }
    }
    return output;
  }

  std::string getIterHeader(void) {
    return this->step_->printHeader();
  }

  std::string getIterInfo(bool withHeader = false) {
    return this->step_->print(*this->state_,withHeader);
  }

}; // class DefaultAlgorithm

} // namespace ROL

#endif
