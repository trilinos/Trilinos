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
struct State {
  Real value;
  Teuchos::RCP<Vector<Real> > gradientVec;
  Real gnorm;
  Teuchos::RCP<Vector<Real> > descentVec;
  Real dnorm;
  Teuchos::RCP<Vector<Real> > iterateVec;
  Real inorm;
};


template <class Real>
class Step {
private:
  State<Real> state_;

public:

  virtual ~Step() {}

  /** \brief Compute step.
  */
  virtual void compute( Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj ) = 0;

  /** \brief Update step, if successful.
  */
  virtual void update( Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj ) = 0;

  /** \brief Print iterate status.
  */
  virtual std::string print( bool printHeader = false ) const = 0;

  // struct StepState (scalars, vectors) map?

  // getState

  // setState

}; // class Step

} // namespace ROL

#endif
