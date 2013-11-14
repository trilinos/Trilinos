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

#ifndef ROL_STEP_H
#define ROL_STEP_H

#include "ROL_Vector.hpp"

/** \class ROL::Step
    \brief Provides the interface to compute optimization steps.
*/


namespace ROL {

template<class Real>
struct AlgorithmState {
  int  iter;
  int  nfval;
  int  ngrad;
  Real value;
  Real gnorm;
  Real snorm;
  Teuchos::RCP<Vector<Real> > iterateVec;
};

template<class Real>
struct StepState {
  Teuchos::RCP<Vector<Real> > gradientVec;
  Teuchos::RCP<Vector<Real> > descentVec;
};


template <class Real>
class Step {
private:

public:
  Teuchos::RCP<StepState<Real> > state_;

  virtual ~Step() {}

  Step(void) { 
    state_ = Teuchos::rcp( new StepState<Real> );
  }

  Teuchos::RCP<StepState<Real> >& get_state() { return this->state_; }

  /** \brief Initialize step.
  */
  virtual void initialize( const Vector<Real> &x, Objective<Real> &obj, AlgorithmState<Real> &algo_state ) {
    state_->descentVec  = x.clone();
    state_->gradientVec = x.clone();
    obj.gradient(*(state_->gradientVec),x);
    algo_state.ngrad = 1;
    algo_state.gnorm = (state_->gradientVec)->norm();
    algo_state.snorm = 1.e10;
    algo_state.value = obj.value(x);
    algo_state.nfval = 1;
  }

  /** \brief Compute step.
  */
  virtual void compute( Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj, 
                        AlgorithmState<Real> &algo_state ) = 0;

  /** \brief Update step, if successful.
  */
  virtual void update( Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj, 
                       AlgorithmState<Real> &algo_state ) = 0;

  /** \brief Print iterate status.
  */
  virtual std::string print( AlgorithmState<Real> &algo_state, bool printHeader = false ) const = 0;

  // struct StepState (scalars, vectors) map?

  // getState

  // setState

}; // class Step

} // namespace ROL

#endif
