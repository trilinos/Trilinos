// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_GRADIENTSTEP_H
#define ROL_GRADIENTSTEP_H

#include "ROL_Types.hpp"
#include "ROL_Step.hpp"
#include "ROL_Secant.hpp"

/** @ingroup step_group
    \class ROL::GradientStep
    \brief Provides the interface to compute optimization steps
           with the gradient descent method globalized using line search.
*/

namespace ROL {

template <class Real>
class GradientStep : public Step<Real> {
private:

  int verbosity_;         ///< Verbosity setting
  const bool computeObj_; ///< Allows line search to compute objective

public:

  using Step<Real>::initialize;
  using Step<Real>::compute;
  using Step<Real>::update;

  /** \brief Constructor.

      Constructor to build a GradientStep object.  Algorithmic specifications
      are passed in through a ROL::ParameterList.

      @param[in]     parlist    is a parameter list containing algorithmic specifications
  */
  GradientStep( ROL::ParameterList &parlist, const bool computeObj = true )
    : Step<Real>(), verbosity_(0), computeObj_(computeObj) {
    // Parse ParameterList
    verbosity_ = parlist.sublist("General").get("Print Verbosity",0);
  }

  void compute( Vector<Real> &s, const Vector<Real> &x,
                Objective<Real> &obj, BoundConstraint<Real> &bnd,
                AlgorithmState<Real> &algo_state ) {
    Real one(1);
    ROL::Ptr<StepState<Real> > step_state = Step<Real>::getState();

    // Compute search direction
    s.set((step_state->gradientVec)->dual());
    s.scale(-one);
  }

  void update( Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj, BoundConstraint<Real> &con,
               AlgorithmState<Real> &algo_state ) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    ROL::Ptr<StepState<Real> > step_state = Step<Real>::getState();

    // Update iterate and store step
    algo_state.iter++;
    x.plus(s);
    (step_state->descentVec)->set(s);
    algo_state.snorm = s.norm();

    // Compute new gradient
    obj.update(x,true,algo_state.iter);
    if ( computeObj_ ) {
      algo_state.value = obj.value(x,tol);
      algo_state.nfval++;
    }
    obj.gradient(*(step_state->gradientVec),x,tol);
    algo_state.ngrad++;

    // Update algorithm state
    (algo_state.iterateVec)->set(x);
    algo_state.gnorm = (step_state->gradientVec)->norm();
  }

  std::string printHeader( void ) const {
    std::stringstream hist;

    if( verbosity_>0 ) {
      hist << std::string(109,'-') <<  "\n";
      hist <<  EDescentToString(DESCENT_STEEPEST);
      hist << " status output definitions\n\n";
      hist << "  iter     - Number of iterates (steps taken) \n";
      hist << "  value    - Objective function value \n";
      hist << "  gnorm    - Norm of the gradient\n";
      hist << "  snorm    - Norm of the step (update to optimization vector)\n";
      hist << "  #fval    - Cumulative number of times the objective function was evaluated\n";
      hist << "  #grad    - Number of times the gradient was computed\n";
      hist << std::string(109,'-') << "\n";
    }

    hist << "  ";
    hist << std::setw(6)  << std::left << "iter";
    hist << std::setw(15) << std::left << "value";
    hist << std::setw(15) << std::left << "gnorm";
    hist << std::setw(15) << std::left << "snorm";
    hist << std::setw(10) << std::left << "#fval";
    hist << std::setw(10) << std::left << "#grad";
    hist << "\n";
    return hist.str();
  }
  std::string printName( void ) const {
    std::stringstream hist;
    hist << "\n" << EDescentToString(DESCENT_STEEPEST) << "\n";
    return hist.str();
  }
  std::string print( AlgorithmState<Real> &algo_state, bool print_header = false ) const {
    std::stringstream hist;
    hist << std::scientific << std::setprecision(6);
    if ( algo_state.iter == 0 ) {
      hist << printName();
    }
    if ( print_header ) {
      hist << printHeader();
    }
    if ( algo_state.iter == 0 ) {
      hist << "  ";
      hist << std::setw(6) << std::left << algo_state.iter;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << "\n";
    }
    else {
      hist << "  ";
      hist << std::setw(6)  << std::left << algo_state.iter;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << std::setw(15) << std::left << algo_state.snorm;
      hist << std::setw(10) << std::left << algo_state.nfval;
      hist << std::setw(10) << std::left << algo_state.ngrad;
      hist << "\n";
    }
    return hist.str();
  }
}; // class GradientStep

} // namespace ROL
#endif
