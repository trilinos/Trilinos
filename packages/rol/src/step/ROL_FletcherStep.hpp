// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_FLETCHERSTEP_H
#define ROL_FLETCHERSTEP_H

#include "ROL_FletcherBase.hpp"
#include "ROL_Step.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_Types.hpp"
#include "ROL_ParameterList.hpp"
/** @ingroup step_group
    \class ROL::FletcherStep
    \brief Provides the interface to compute Fletcher steps.
*/


namespace ROL {

template <class Real>
class FletcherStep : public Step<Real> {
private:
  ROL::Ptr<Step<Real> > step_;
  ROL::Ptr<BoundConstraint<Real> > bnd_;

  ROL::ParameterList parlist_;

  ROL::Ptr<Vector<Real> > x_; 

  // Lagrange multiplier update
  Real penaltyUpdate_;
  bool modifyPenalty_;
  Real maxPenaltyParam_;
  Real minPenaltyParam_;
  // Subproblem information
  bool print_;
  std::string subStep_;

  Real delta_;
  Real deltaMin_;
  Real deltaUpdate_;
  ETrustRegion etr_;

  bool bnd_activated_;

  ROL::Ptr<Vector<Real> > g_;

  int numSuccessSteps_;

  // For printing output
  mutable bool isDeltaChanged_;
  mutable bool isPenaltyChanged_;

  mutable AlgorithmState<Real> tr_algo_state_;

  mutable int stepHeaderLength_; // For formatting

  Real computeProjGradientNorm(const Vector<Real> &g, const Vector<Real> &x,
                               BoundConstraint<Real> &bnd) {
    Real gnorm = 0.;
    // Compute norm of projected gradient
    if (bnd.isActivated()) {
      x_->set(x);
      x_->axpy(-1.,g.dual());
      bnd.project(*x_);
      x_->axpy(-1.,x);
      gnorm = x_->norm();
    }
    else {
      gnorm = g.norm();
    }
    return gnorm;
  }

public:

  using Step<Real>::initialize;
  using Step<Real>::compute;
  using Step<Real>::update;

  ~FletcherStep() {}

  FletcherStep(ROL::ParameterList &parlist)
    : Step<Real>(), bnd_activated_(false), numSuccessSteps_(0),
      isDeltaChanged_(true), isPenaltyChanged_(true), stepHeaderLength_(0) {
    Real zero(0), one(1), two(2), oe8(1.e8), oe1(1.e-1), oem6(1e-6), oem8(1.e-8);

    ROL::ParameterList& sublist = parlist.sublist("Step").sublist("Fletcher");
    Step<Real>::getState()->searchSize = sublist.get("Penalty Parameter",one);
    delta_ = sublist.get("Regularization Parameter",zero);
    deltaMin_ = sublist.get("Min Regularization Parameter",oem8);
    deltaUpdate_ = sublist.get("Regularization Parameter Decrease Factor", oe1);
    // penalty parameters
    penaltyUpdate_ = sublist.get("Penalty Parameter Growth Factor", two);
    modifyPenalty_ = sublist.get("Modify Penalty Parameter", false);
    maxPenaltyParam_ = sublist.get("Maximum Penalty Parameter", oe8);
    minPenaltyParam_ = sublist.get("Minimum Penalty Parameter", oem6);     

    subStep_ = sublist.get("Subproblem Solver", "Trust Region");

    parlist_ = parlist;
  }

  /** \brief Initialize step with equality constraint.
  */
  void initialize( Vector<Real> &x, const Vector<Real> &g, Vector<Real> &l, const Vector<Real> &c,
                   Objective<Real> &obj, Constraint<Real> &con,
                   AlgorithmState<Real> &algo_state ) {
    bnd_ = ROL::makePtr<BoundConstraint<Real>>();
    bnd_->deactivate();
    initialize(x,g,l,c,obj,con,*bnd_,algo_state);
  }

  /** \brief Initialize step with equality and bound constraints.
  */
  void initialize( Vector<Real> &x, const Vector<Real> &g, Vector<Real> &l, const Vector<Real> &c,
                   Objective<Real> &obj, Constraint<Real> &con, BoundConstraint<Real> &bnd,
                   AlgorithmState<Real> &algo_state ) {
    // Determine what kind of step
    bnd_activated_ = bnd.isActivated();

    ROL::ParameterList trlist(parlist_);
    bool inexactFletcher = trlist.sublist("Step").sublist("Fletcher").get("Inexact Solves", false);
    if( inexactFletcher ) {
      trlist.sublist("General").set("Inexact Objective Value", true);
      trlist.sublist("General").set("Inexact Gradient", true);
    }
    if( bnd_activated_ ) {
      trlist.sublist("Step").sublist("Trust Region").set("Subproblem Model", "Coleman-Li");
    }

    if ( subStep_ == "Line Search" ) {
      step_ = makePtr<LineSearchStep<Real>>(trlist);
    }
    else {
      step_ = makePtr<TrustRegionStep<Real>>(trlist);
    }
    etr_ = StringToETrustRegion(parlist_.sublist("Step").sublist("Trust Region").get("Subproblem Solver", "Truncated CG"));

    // Initialize class members
    g_ = g.clone();
    x_ = x.clone();

    // Rest of initialize
    FletcherBase<Real>& fletcher = dynamic_cast<FletcherBase<Real>&>(obj);

    tr_algo_state_.iterateVec = x.clone();
    tr_algo_state_.minIterVec = x.clone();
    tr_algo_state_.lagmultVec = l.clone();

    step_->initialize(x, g, obj, bnd, tr_algo_state_);

    // Initialize step state
    ROL::Ptr<StepState<Real> > state = Step<Real>::getState();
    state->descentVec    = x.clone();
    state->gradientVec   = g.clone();
    state->constraintVec = c.clone();
    // Initialize the algorithm state
    algo_state.nfval = 0;
    algo_state.ncval = 0;
    algo_state.ngrad = 0;

    algo_state.value = fletcher.getObjectiveValue(x);
    algo_state.gnorm = computeProjGradientNorm(*(fletcher.getLagrangianGradient(x)),
                                               x, bnd);
    algo_state.aggregateGradientNorm = tr_algo_state_.gnorm;

    state->constraintVec->set(*(fletcher.getConstraintVec(x)));
    algo_state.cnorm = (state->constraintVec)->norm();
    // Update evaluation counters
    algo_state.ncval = fletcher.getNumberConstraintEvaluations();
    algo_state.nfval = fletcher.getNumberFunctionEvaluations();
    algo_state.ngrad = fletcher.getNumberGradientEvaluations();
  }

  /** \brief Compute step (equality constraint).
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, const Vector<Real> &l,
                Objective<Real> &obj, Constraint<Real> &con, 
                AlgorithmState<Real> &algo_state ) {
    compute(s,x,l,obj,con,*bnd_, algo_state);
  }

  /** \brief Compute step (equality and bound constraints).
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, const Vector<Real> &l,
                Objective<Real> &obj, Constraint<Real> &con, 
                BoundConstraint<Real> &bnd, AlgorithmState<Real> &algo_state ) {
    step_->compute( s, x, obj, bnd, tr_algo_state_ );  
  }

  /** \brief Update step, if successful (equality constraint).
  */
  void update( Vector<Real> &x, Vector<Real> &l, const Vector<Real> &s,
               Objective<Real> &obj, Constraint<Real> &con,
               AlgorithmState<Real> &algo_state ) {
    update(x,l,s,obj,con,*bnd_, algo_state);
  }

  /** \brief Update step, if successful (equality and bound constraints).
  */
  void update( Vector<Real> &x, Vector<Real> &l, const Vector<Real> &s,
               Objective<Real> &obj, Constraint<Real> &con,
               BoundConstraint<Real> &bnd,
               AlgorithmState<Real> &algo_state ) {
    
    // This should be in print, but this will not work there
    isDeltaChanged_ = false;
    isPenaltyChanged_ = false;
    bool modified = false;

    FletcherBase<Real> &fletcher = dynamic_cast<FletcherBase<Real>&>(obj);
    ROL::Ptr<StepState<Real> > fletcherState = Step<Real>::getState();
    const ROL::Ptr<const StepState<Real> > state = step_->getStepState();

    step_->update(x,s,obj,bnd,tr_algo_state_);
    numSuccessSteps_ += (state->flag == 0);

    Real gPhiNorm = tr_algo_state_.gnorm;
    Real cnorm = (fletcherState->constraintVec)->norm();
    bool too_infeasible = cnorm > static_cast<Real>(100.)*gPhiNorm;
    bool too_feasible = cnorm < static_cast<Real>(1e-2)*gPhiNorm;

    if( too_infeasible && !modified && modifyPenalty_ && numSuccessSteps_ > 1 ) {
      Real penaltyParam = Step<Real>::getStepState()->searchSize;
      if( penaltyParam >= maxPenaltyParam_ ) {
        // Penalty parameter too large, exit
        algo_state.flag = true;
      }
      penaltyParam *= penaltyUpdate_;
      penaltyParam = std::min(penaltyParam, maxPenaltyParam_);
      fletcher.setPenaltyParameter(penaltyParam);
      Step<Real>::getState()->searchSize = penaltyParam;
      isPenaltyChanged_ = true;
      modified = true;
    }

    if( too_feasible && !modified && modifyPenalty_ && numSuccessSteps_ > 1 ) {
      Real penaltyParam = Step<Real>::getStepState()->searchSize;
      if( penaltyParam <= minPenaltyParam_ ) {
        // Penalty parameter too small, exit (this is unlikely)
        algo_state.flag = true;
      }
      penaltyParam /= penaltyUpdate_;
      penaltyParam = std::max(penaltyParam, minPenaltyParam_);
      fletcher.setPenaltyParameter(penaltyParam);
      Step<Real>::getState()->searchSize = penaltyParam;
      isPenaltyChanged_ = true;
      modified = true;      
    }

    if( delta_ > deltaMin_ && !modified ) {
      Real deltaNext = delta_ * deltaUpdate_;
      if( gPhiNorm < deltaNext ) {
        delta_ = deltaNext;
        fletcher.setDelta(deltaNext);
        isDeltaChanged_ = true;
        modified = true;
      }
    }

    if( modified ) {
      // Penalty function has been changed somehow, need to recompute
      Real tol = static_cast<Real>(1e-12);
      tr_algo_state_.value = fletcher.value(x, tol);
      fletcher.gradient(*g_, x, tol);

      tr_algo_state_.nfval++;
      tr_algo_state_.ngrad++;
      tr_algo_state_.ncval++;
      tr_algo_state_.minIter = tr_algo_state_.iter;
      tr_algo_state_.minValue = tr_algo_state_.value;
      tr_algo_state_.gnorm = computeProjGradientNorm(*g_, x, bnd);
    }

    // Update the step and store in state
    algo_state.iterateVec->set(x);
    algo_state.iter++;

    fletcherState->descentVec->set(s);
    fletcherState->gradientVec->set(*(fletcher.getLagrangianGradient(x)));
    fletcherState->constraintVec->set(*(fletcher.getConstraintVec(x)));

    // Update objective function value
    algo_state.value = fletcher.getObjectiveValue(x);
    // Update constraint value
    algo_state.cnorm = (fletcherState->constraintVec)->norm();
    // Update the step size
    algo_state.snorm = tr_algo_state_.snorm;
    // Compute gradient of the Lagrangian
    algo_state.gnorm = computeProjGradientNorm(*(fletcherState->gradientVec),
                                               x, bnd);
    // Compute gradient of penalty function
    algo_state.aggregateGradientNorm = tr_algo_state_.gnorm;
    // Update evaluation countersgetConstraintVec
    algo_state.nfval = fletcher.getNumberFunctionEvaluations();
    algo_state.ngrad = fletcher.getNumberGradientEvaluations();
    algo_state.ncval = fletcher.getNumberConstraintEvaluations();
    // Update objective function and constraints
    // fletcher.update(x,true,algo_state.iter);
    // Update multipliers
    algo_state.lagmultVec->set(*(fletcher.getMultiplierVec(x)));
  }

  /** \brief Print iterate header.
  */
  std::string printHeader( void ) const {
    std::stringstream hist;
    if( subStep_ == "Trust Region" ) {
      hist << "  ";
      hist << std::setw(6)  << std::left << "iter";
      hist << std::setw(15) << std::left << "merit";
      hist << std::setw(15) << std::left << "fval";
      hist << std::setw(15) << std::left << "gpnorm";
      hist << std::setw(15) << std::left << "gLnorm";
      hist << std::setw(15) << std::left << "cnorm";
      hist << std::setw(15) << std::left << "snorm";
      hist << std::setw(15) << std::left << "tr_radius";
      hist << std::setw(10) << std::left << "tr_flag";
      if ( etr_ == TRUSTREGION_TRUNCATEDCG && subStep_ == "Trust Region") {
        hist << std::setw(10) << std::left << "iterCG";
        hist << std::setw(10) << std::left << "flagCG";
      }
      hist << std::setw(15) << std::left << "penalty";
      hist << std::setw(15) << std::left << "delta";
      hist << std::setw(10) << std::left << "#fval";
      hist << std::setw(10) << std::left << "#grad";
      hist << std::setw(10) << std::left << "#cval";
      hist << "\n"; 
    }
    else {
      std::string stepHeader = step_->printHeader();
      stepHeaderLength_ = stepHeader.length();
      hist << stepHeader.substr(0, stepHeaderLength_-1);
      hist << std::setw(15) << std::left << "fval";
      hist << std::setw(15) << std::left << "gLnorm";
      hist << std::setw(15) << std::left << "cnorm";
      hist << std::setw(15) << std::left << "penalty";
      hist << std::setw(15) << std::left << "delta";
      hist << std::setw(10) << std::left << "#cval";
      hist << "\n";
    }
    return hist.str();
  }

  /** \brief Print step name.
  */
  std::string printName( void ) const {
    std::stringstream hist;
    hist << "\n" << " Fletcher solver : " << subStep_;
    hist << "\n";
    return hist.str();
  }

  /** \brief Print iterate status.
  */
  std::string print( AlgorithmState<Real> &algo_state, bool pHeader = false ) const {
    std::string stepHist = step_->print( tr_algo_state_, false );
    stepHist.erase(std::remove(stepHist.end()-3, stepHist.end(),'\n'), stepHist.end());
    std::string name = step_->printName();
    size_t pos = stepHist.find(name);
    if ( pos != std::string::npos ) {
      stepHist.erase(pos, name.length());
    }

    std::stringstream hist;
    hist << std::scientific << std::setprecision(6);
    if ( algo_state.iter == 0 ) {
      hist << printName();
    }
    if ( pHeader ) {
      hist << printHeader();
    }

    std::string penaltyString = getValueString( Step<Real>::getStepState()->searchSize, isPenaltyChanged_ );
    std::string deltaString = getValueString( delta_, isDeltaChanged_ );

    if( subStep_ == "Trust Region" ) {
      hist << "  ";
      hist << std::setw(6)  << std::left << algo_state.iter;
      hist << std::setw(15) << std::left << tr_algo_state_.value;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << tr_algo_state_.gnorm;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << std::setw(15) << std::left << algo_state.cnorm;
      hist << std::setw(15) << std::left << stepHist.substr(38,15); // snorm
      hist << std::setw(15) << std::left << stepHist.substr(53,15); // tr_radius
      hist << std::setw(10) << std::left << (algo_state.iter == 0 ? "" : stepHist.substr(88,10)); // tr_flag
      if ( etr_ == TRUSTREGION_TRUNCATEDCG && subStep_ == "Trust Region") {
        hist << std::setw(10) << std::left << (algo_state.iter == 0 ? "" : stepHist.substr(93,10)); // iterCG
        hist << std::setw(10) << std::left << (algo_state.iter == 0 ? "" : stepHist.substr(103,10)); // flagCG
      }
      hist << std::setw(15) << std::left << penaltyString;
      hist << std::setw(15) << std::left << deltaString;
      hist << std::setw(10) << std::left << (algo_state.iter == 0 ? "" : stepHist.substr(68,10)); // #fval
      hist << std::setw(10) << std::left << (algo_state.iter == 0 ? "" : stepHist.substr(78,10)); // #gval
      hist << std::setw(10) << std::left << algo_state.ncval;
      hist << "\n"; 
    } else {
      hist << std::setw(stepHeaderLength_-1) << std::left << stepHist;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << std::setw(15) << std::left << algo_state.cnorm;
      hist << std::setw(15) << std::left << penaltyString;
      hist << std::setw(15) << std::left << deltaString;
      hist << std::setw(10) << std::left << algo_state.ncval;
      hist << "\n";
    }

    return hist.str();
  }

  std::string getValueString( const Real value, const bool print ) const {
    std::stringstream valString;
    valString << std::scientific << std::setprecision(6);
    if( print ) {
      valString << std::setw(15) << std::left << value;
    } else {
      valString << std::setw(15) << "";
    }
    return valString.str();
  }

  /** \brief Compute step for bound constraints; here only to satisfy the
             interface requirements, does nothing, needs refactoring.
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj,
                        BoundConstraint<Real> &con,
                        AlgorithmState<Real> &algo_state ) {}

  /** \brief Update step, for bound constraints; here only to satisfy the
             interface requirements, does nothing, needs refactoring.
  */
  void update( Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj,
                       BoundConstraint<Real> &con,
                       AlgorithmState<Real> &algo_state ) {}

}; // class FletcherStep

} // namespace ROL

#endif
