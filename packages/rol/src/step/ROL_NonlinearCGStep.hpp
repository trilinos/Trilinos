// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_NONLINEARCGSTEP_H
#define ROL_NONLINEARCGSTEP_H

#include "ROL_Types.hpp"
#include "ROL_Step.hpp"
#include "ROL_NonlinearCG.hpp"

/** @ingroup step_group
    \class ROL::NonlinearCGStep
    \brief Provides the interface to compute optimization steps
           with nonlinear CG.
*/

namespace ROL {

template <class Real>
class NonlinearCGStep : public Step<Real> {
private:

  ROL::Ptr<NonlinearCG<Real> > nlcg_; ///< NonlinearCG object (used for quasi-Newton)
  ENonlinearCG enlcg_;
  int verbosity_;                         ///< Verbosity setting
  const bool computeObj_;

  std::string ncgName_;

public:

  using Step<Real>::initialize;
  using Step<Real>::compute;
  using Step<Real>::update;

  /** \brief Constructor.

      Constructor to build a NonlinearCGStep object with a user-defined 
      nonlinear CG object.  Algorithmic specifications are passed in through 
      a ROL::ParameterList.

      @param[in]     parlist    is a parameter list containing algorithmic specifications
      @param[in]     nlcg       is a user-defined NonlinearCG object
  */
  NonlinearCGStep( ROL::ParameterList &parlist,
             const ROL::Ptr<NonlinearCG<Real> > &nlcg = ROL::nullPtr,
             const bool computeObj = true )
    : Step<Real>(), nlcg_(nlcg), enlcg_(NONLINEARCG_USERDEFINED),
      verbosity_(0), computeObj_(computeObj) {
    // Parse ParameterList
    verbosity_ = parlist.sublist("General").get("Print Verbosity",0);
    // Initialize secant object
    ROL::ParameterList& Llist = parlist.sublist("Step").sublist("Line Search");
    if ( nlcg == ROL::nullPtr ) {
      ncgName_ = Llist.sublist("Descent Method").get("Nonlinear CG Type","Oren-Luenberger");
      enlcg_
        = StringToENonlinearCG(ncgName_);
      nlcg_ = ROL::makePtr<NonlinearCG<Real>>(enlcg_);
    }
    else {
      ncgName_ = Llist.sublist("Descent Method").get("User Defined Nonlinear CG Name",
                                                     "Unspecified User Define Nonlinear CG Method");
    }
  }

  void compute( Vector<Real> &s, const Vector<Real> &x,
                Objective<Real> &obj, BoundConstraint<Real> &bnd,
                AlgorithmState<Real> &algo_state ) {
    ROL::Ptr<StepState<Real> > step_state = Step<Real>::getState();
    Real one(1);

    // Compute search direction
    nlcg_->run(s,*(step_state->gradientVec),x,obj);
    s.scale(-one);
  }

  void update( Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj, BoundConstraint<Real> &con,
               AlgorithmState<Real> &algo_state ) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    ROL::Ptr<StepState<Real> > step_state = Step<Real>::getState();

    // Update iterate
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
      hist << EDescentToString(DESCENT_NONLINEARCG);
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
    hist << "\n" << ncgName_ << " "
         << EDescentToString(DESCENT_NONLINEARCG) << "\n";
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
}; // class NonlinearCGStep

} // namespace ROL

#endif
