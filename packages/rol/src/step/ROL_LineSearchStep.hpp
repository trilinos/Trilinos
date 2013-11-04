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

/** \class ROL::LineSearchStep
    \brief Provides the interface to compute optimization steps
           with line search.
*/


namespace ROL {

template <class Real>
class LineSearchStep : public Step<Real> {
public:

  virtual ~LineSearchStep() {}

  /** \brief Compute step.
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj ) {
    obj.gradient(s, x);
    s.scale(-1.0);

    // Perform backtracking line search from Nocedal/Wright.
    int maxit = 20;
    int it = 0;
    Real rho = 0.5;
    Real alpha = 1.0;
    Real c = 1e-4;
    Real fnew = 0.0;
    Real f = 0.0;
    Teuchos::RCP<Vector<Real> > xnew = x.clone();

    xnew->set(x);
    xnew->axpy(alpha, s);
    fnew = obj.value(*xnew);

    while ( (fnew > f + c*alpha*s.norm()*s.norm()) && (it < maxit) ) {
      alpha *= rho;
      xnew->set(x);
      xnew->axpy(alpha, s);
      fnew = obj.value(*xnew);
      it++;
    }
    s.scale(alpha);

  }

  /** \brief Update step, if successful.
  */
  void update( Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj ) {
    x.axpy(1.0, s);
  }

  /** \brief Print iterate status.
  */
  std::string print( bool printHeader = false ) const  {
    return "LineSearchStep";
  }

  // struct StepState (scalars, vectors) map?

  // getState

  // setState

}; // class Step

} // namespace ROL

#endif
