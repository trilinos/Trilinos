// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PROJECTEDSECANTSTEP_H
#define ROL_PROJECTEDSECANTSTEP_H

#include "ROL_Types.hpp"
#include "ROL_Step.hpp"
#include "ROL_Secant.hpp"

/** @ingroup step_group
    \class ROL::ProjectedSecantStep
    \brief Provides the interface to compute optimization steps
           with projected secant method using line search.
*/

namespace ROL {

template <class Real>
class ProjectedSecantStep : public Step<Real> {
private:

  ROL::Ptr<Secant<Real> > secant_; ///< Secant object (used for quasi-Newton)
  ESecant esec_;
  ROL::Ptr<Vector<Real> > d_;      ///< Additional vector storage
  ROL::Ptr<Vector<Real> > gp_;     ///< Additional vector storage
  int verbosity_;                      ///< Verbosity level
  const bool computeObj_;
  bool useProjectedGrad_;              ///< Whether or not to use to the projected gradient criticality measure

public:

  using Step<Real>::initialize;
  using Step<Real>::compute;
  using Step<Real>::update;

  /** \brief Constructor.

      Standard constructor to build a ProjectedSecantStep object.  Algorithmic 
      specifications are passed in through a ROL::ParameterList.

      @param[in]     parlist    is a parameter list containing algorithmic specifications
      @param[in]     secant     is a user-defined secant object
  */
  ProjectedSecantStep( ROL::ParameterList &parlist,
                       const ROL::Ptr<Secant<Real> > &secant = ROL::nullPtr, 
                       const bool computeObj = true )
    : Step<Real>(), secant_(secant), d_(ROL::nullPtr), gp_(ROL::nullPtr),
      verbosity_(0), computeObj_(computeObj), useProjectedGrad_(false) {
    // Parse ParameterList
    ROL::ParameterList& Glist = parlist.sublist("General");
    useProjectedGrad_ = Glist.get("Projected Gradient Criticality Measure", false);
    verbosity_ = parlist.sublist("General").get("Print Verbosity",0);
    // Initialize secant object
    if ( secant == ROL::nullPtr ) {
      esec_ = StringToESecant(parlist.sublist("General").sublist("Secant").get("Type","Limited-Memory BFGS"));
      secant_ = SecantFactory<Real>(parlist);
    }
  }

  void initialize( Vector<Real> &x, const Vector<Real> &s, const Vector<Real> &g, 
                   Objective<Real> &obj, BoundConstraint<Real> &bnd, 
                   AlgorithmState<Real> &algo_state ) {
    Step<Real>::initialize(x,s,g,obj,bnd,algo_state);
    d_  = s.clone();
    gp_ = g.clone();
  }

  void compute( Vector<Real> &s, const Vector<Real> &x,
                Objective<Real> &obj, BoundConstraint<Real> &bnd,
                AlgorithmState<Real> &algo_state ) {
    ROL::Ptr<StepState<Real> > step_state = Step<Real>::getState();
    Real one(1);

    // Compute projected secant step
    // ---> Apply inactive-inactive block of inverse secant to gradient
    gp_->set(*(step_state->gradientVec));
    bnd.pruneActive(*gp_,*(step_state->gradientVec),x,algo_state.gnorm);
    secant_->applyH(s,*gp_);
    bnd.pruneActive(s,*(step_state->gradientVec),x,algo_state.gnorm);
    // ---> Add in active gradient components
    gp_->set(*(step_state->gradientVec));
    bnd.pruneInactive(*gp_,*(step_state->gradientVec),x,algo_state.gnorm);
    s.plus(gp_->dual());
    s.scale(-one);
  }

  void update( Vector<Real> &x, const Vector<Real> &s,
               Objective<Real> &obj, BoundConstraint<Real> &bnd,
               AlgorithmState<Real> &algo_state ) {
    Real tol = std::sqrt(ROL_EPSILON<Real>()), one(1);
    ROL::Ptr<StepState<Real> > step_state = Step<Real>::getState();

    // Update iterate and store previous step
    algo_state.iter++;
    d_->set(x);
    x.plus(s);
    bnd.project(x);
    (step_state->descentVec)->set(x);
    (step_state->descentVec)->axpy(-one,*d_);
    algo_state.snorm = s.norm();

    // Compute new gradient
    gp_->set(*(step_state->gradientVec));
    obj.update(x,true,algo_state.iter);
    if ( computeObj_ ) {
      algo_state.value = obj.value(x,tol);
      algo_state.nfval++;
    }
    obj.gradient(*(step_state->gradientVec),x,tol);
    algo_state.ngrad++;

    // Update Secant Information
    secant_->updateStorage(x,*(step_state->gradientVec),*gp_,s,algo_state.snorm,algo_state.iter+1);

    // Update algorithm state
    (algo_state.iterateVec)->set(x);
    if ( useProjectedGrad_ ) {
      gp_->set(*(step_state->gradientVec));
      bnd.computeProjectedGradient( *gp_, x );
      algo_state.gnorm = gp_->norm();
    }
    else {
      d_->set(x);
      d_->axpy(-one,(step_state->gradientVec)->dual());
      bnd.project(*d_);
      d_->axpy(-one,x);
      algo_state.gnorm = d_->norm();
    }
  }

  std::string printHeader( void ) const {
    std::stringstream hist;

    if( verbosity_>0 ) {
      hist << std::string(109,'-') <<  "\n";
      hist << EDescentToString(DESCENT_SECANT);
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
    hist << "\n" << EDescentToString(DESCENT_SECANT);
    hist << " with " << ESecantToString(esec_) << "\n";
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
}; // class ProjectedSecantStep

} // namespace ROL

#endif
