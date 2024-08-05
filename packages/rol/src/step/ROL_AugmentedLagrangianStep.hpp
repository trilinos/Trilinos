// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_AUGMENTEDLAGRANGIANSTEP_H
#define ROL_AUGMENTEDLAGRANGIANSTEP_H

#include "ROL_AugmentedLagrangian.hpp"
#include "ROL_Types.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_ParameterList.hpp"

// Step (bound constrained or unconstrained) includes
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_PrimalDualActiveSetStep.hpp"
#include "ROL_MoreauYosidaPenaltyStep.hpp"
#include "ROL_BundleStep.hpp"
#include "ROL_InteriorPointStep.hpp"

// StatusTest includes
#include "ROL_StatusTest.hpp"
#include "ROL_BundleStatusTest.hpp"

/** @ingroup step_group
    \class ROL::AugmentedLagrangianStep
    \brief Provides the interface to compute augmented Lagrangian steps.

    The Augmented Lagrangian algorithm is used to solve Type-EB problems.
    This algorithm solves the scaled problem
    \f[
        \min_{x} w_J J(x) \quad\text{subject to}\quad
           w_c c(x) = 0,\quad \ell \le x \le u
    \f]
    for some positive consants \f$w_J,\, w_c\f$.  These constants are
    either input by the user or automatically estimated.  To solve this
    scaled problem, the Augmented Lagrangian algorithm minimizes the
    so-called Augmented Lagrangian functional
    \f[
        L(x,\lambda,r) := w_J J(x) + w_c \langle \lambda, c(x)\rangle_{X^*,X}
              + \frac{w_c^2 r}{2} \|c(x)\|_X^2
    \f]
    subject to the bound constraints \f$\ell \le x \le u\f$.  The multiplier
    estimate \f$\lambda\f$ is updated as
    \f[
       \lambda \leftarrow \lambda + r w_c c(x).
    \f]
    The penalty parameter \f$r>0\f$ is also updated based on the progress
    of the algorithm.  The initial penalty parameter is either input by
    the user or automatically computed.

    User Input Parameters:
    bool    Step -> Augmented Lagrangian -> Use Default initial Penalty Parameter
      Use automatically determined initial penalty parameter.  Default: true

    Real    Step -> Augmented Lagrangian -> Initial Penalty Parameter
      Initial penalty parameter.  Default: 10

    Real    Step -> Augmented Lagrangian -> Use Scaled Augmented Lagrangian
      Use Augmented Lagrangian scaled by the reciprocal of the penalty parameter.  Default: false

    Real    Step -> Augmented Lagrangian -> Penalty Parameter Reciprocal Lower Bound
      Minimum penalty parameter reciprocal for tolerance updates.  Default: 0.1

    Real    Step -> Augmented Lagrangian -> Penalty Parameter Growth Factor
      Rate of growth for penalty parameter.  Default: 10

    Real    Step -> Augmented Lagrangian -> Maximum Penalty Parameter
      Maximum penalty parameter size.  Default: 1e8

    Real    Step -> Augmented Lagrangian -> Optimality Tolerance Update Exponent
      Rate at which to update optimality tolerance.  Default: 1

    Real    Step -> Augmented Lagrangian -> Optimality Tolerance Decrease Exponent
      Rate at which to decrease optimality tolerance.  Default: 1

    Real    Step -> Augmented Lagrangian -> Initial Optimality Tolerance
      Initial tolerance for optimality.  Default: 1

    Real    Step -> Augmented Lagrangian -> Feasibility Tolerance Update Exponent
      Rate at which to update feasibility tolerance.  Default: 0.1

    Real    Step -> Augmented Lagrangian -> Feasibility Tolerance Decrease Exponent
      Rate at which to decrease feasibility tolerance.  Default: 0.9

    Real    Step -> Augmented Lagrangian -> Initial Feasibility Tolerance
      Initial tolerance for equality constraint feasibility.  Default: 1

    bool    Step -> Augmented Lagrangian -> Print Intermediate Optimization History
      Print iteration history for subproblem solve.  Default: false

    int     Step -> Augmented Lagrangian -> Subproblem Iteration Limit
      Subproblem iteration limit.  Default: 1000

    string  Step -> Augmented Lagrangian -> Subproblem Step Type
      Subproblem (bound constrained) solver type. Default: Trust Region

    bool    Step -> Augmented Lagrangian -> Use Default Problem Scaling
      Use automatic constraint and objective scaling.  Default: true

    Real    Step -> Augmented Lagrangian -> Objective Scaling
      Positive scaling constant for objective. Default: 1

    Real    Step -> Augmented Lagrangian -> Constraint Scaling
      Positive scaling constant for constraint. Default: 1

    int     General -> Print Verbosity
      Print additional information to screen for debugging purposes.  Default: 0
*/


namespace ROL {

template<class Real>
class MoreauYosidaPenaltyStep;

template<class Real>
class InteriorPointStep;

template <class Real>
class AugmentedLagrangianStep : public Step<Real> {
private:
  ROL::Ptr<StatusTest<Real>>      status_;
  ROL::Ptr<Step<Real>>            step_;
  ROL::Ptr<Algorithm<Real>>       algo_;
  ROL::Ptr<Vector<Real>>          x_; 
  ROL::Ptr<BoundConstraint<Real>> bnd_;

  ROL::ParameterList parlist_;
  // Lagrange multiplier update
  bool useDefaultInitPen_;
  bool scaleLagrangian_;
  Real minPenaltyReciprocal_;
  Real minPenaltyLowerBound_;
  Real penaltyUpdate_;
  Real maxPenaltyParam_;
  // Optimality tolerance update
  Real optIncreaseExponent_;
  Real optDecreaseExponent_;
  Real optToleranceInitial_;
  Real optTolerance_;
  // Feasibility tolerance update
  Real feasIncreaseExponent_;
  Real feasDecreaseExponent_;
  Real feasToleranceInitial_;
  Real feasTolerance_;
  // Subproblem information
  bool print_;
  int maxit_;
  int subproblemIter_;
  std::string subStep_;
  Real outerOptTolerance_;
  Real outerFeasTolerance_;
  Real outerStepTolerance_;
  // Scaling information
  bool useDefaultScaling_;
  Real fscale_;
  Real cscale_;
  // Verbosity flag
  int verbosity_;

  Real computeGradient(Vector<Real> &g, const Vector<Real> &x,
                       const Real mu, Objective<Real> &obj,
                       BoundConstraint<Real> &bnd) {
    AugmentedLagrangian<Real> &augLag
      = dynamic_cast<AugmentedLagrangian<Real>&>(obj);
    Real gnorm = 0., tol = std::sqrt(ROL_EPSILON<Real>());
    augLag.gradient(g,x,tol);
    if ( scaleLagrangian_ ) {
      g.scale(mu);
    }
    // Compute norm of projected gradient
    if (bnd.isActivated()) {
      x_->set(x);
      x_->axpy(static_cast<Real>(-1),g.dual());
      bnd.project(*x_);
      x_->axpy(static_cast<Real>(-1),x);
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

  ~AugmentedLagrangianStep() {}

  AugmentedLagrangianStep(ROL::ParameterList &parlist)
    : Step<Real>(), algo_(ROL::nullPtr),
      x_(ROL::nullPtr), parlist_(parlist), subproblemIter_(0) {
    Real one(1), p1(0.1), p9(0.9), ten(1.e1), oe8(1.e8), oem8(1.e-8);
    ROL::ParameterList& sublist = parlist.sublist("Step").sublist("Augmented Lagrangian");
    useDefaultInitPen_    = sublist.get("Use Default Initial Penalty Parameter",true);
    Step<Real>::getState()->searchSize = sublist.get("Initial Penalty Parameter",ten);
    // Multiplier update parameters
    scaleLagrangian_      = sublist.get("Use Scaled Augmented Lagrangian",          false);
    minPenaltyLowerBound_ = sublist.get("Penalty Parameter Reciprocal Lower Bound", p1);
    minPenaltyReciprocal_ = p1;
    penaltyUpdate_        = sublist.get("Penalty Parameter Growth Factor",          ten);
    maxPenaltyParam_      = sublist.get("Maximum Penalty Parameter",                oe8);
    // Optimality tolerance update
    optIncreaseExponent_ = sublist.get("Optimality Tolerance Update Exponent",    one);
    optDecreaseExponent_ = sublist.get("Optimality Tolerance Decrease Exponent",  one);
    optToleranceInitial_ = sublist.get("Initial Optimality Tolerance",            one);
    // Feasibility tolerance update    
    feasIncreaseExponent_ = sublist.get("Feasibility Tolerance Update Exponent",   p1);
    feasDecreaseExponent_ = sublist.get("Feasibility Tolerance Decrease Exponent", p9);
    feasToleranceInitial_ = sublist.get("Initial Feasibility Tolerance",           one);
    // Subproblem information
    print_   = sublist.get("Print Intermediate Optimization History", false);
    maxit_   = sublist.get("Subproblem Iteration Limit",              1000);
    subStep_ = sublist.get("Subproblem Step Type",                    "Trust Region");
    parlist_.sublist("Step").set("Type",subStep_);
    parlist_.sublist("Status Test").set("Iteration Limit",maxit_);
    // Verbosity setting
    verbosity_          = parlist.sublist("General").get("Print Verbosity", 0);
    print_              = (verbosity_ > 0 ? true : print_);
    // Outer iteration tolerances
    outerFeasTolerance_ = parlist.sublist("Status Test").get("Constraint Tolerance", oem8);
    outerOptTolerance_  = parlist.sublist("Status Test").get("Gradient Tolerance",   oem8);
    outerStepTolerance_ = parlist.sublist("Status Test").get("Step Tolerance",       oem8);
    // Scaling
    useDefaultScaling_  = sublist.get("Use Default Problem Scaling", true);
    fscale_             = sublist.get("Objective Scaling", 1.0);
    cscale_             = sublist.get("Constraint Scaling", 1.0);
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
    Real one(1), TOL(1.e-2);
    AugmentedLagrangian<Real> &augLag
      = dynamic_cast<AugmentedLagrangian<Real>&>(obj);
    // Initialize step state
    ROL::Ptr<StepState<Real> > state = Step<Real>::getState();
    state->descentVec    = x.clone();
    state->gradientVec   = g.clone();
    state->constraintVec = c.clone();
    // Initialize additional storage
    x_ = x.clone();
    // Initialize the algorithm state
    algo_state.nfval = 0;
    algo_state.ncval = 0;
    algo_state.ngrad = 0;
    // Project x onto the feasible set
    if ( bnd.isActivated() ) {
      bnd.project(x);
    }
    // Update objective and constraint.
    augLag.update(x,true,algo_state.iter);
    if (useDefaultScaling_) {
      fscale_ = one/std::max(one,augLag.getObjectiveGradient(x)->norm());
      try {
        Real tol = std::sqrt(ROL_EPSILON<Real>());
        Ptr<Vector<Real>> ji = x.clone();
        Real maxji(0), normji(0);
        for (int i = 0; i < c.dimension(); ++i) {
          con.applyAdjointJacobian(*ji,*c.basis(i),x,tol);
          normji = ji->norm();
          maxji  = std::max(normji,maxji);
        }
        cscale_ = one/std::max(one,maxji);
      }
      catch (std::exception &e) {
        cscale_ = one;
      }
    }
    augLag.setScaling(fscale_,cscale_);
    algo_state.value = augLag.getObjectiveValue(x);
    algo_state.gnorm = computeGradient(*(state->gradientVec),x,state->searchSize,obj,bnd);
    augLag.getConstraintVec(*(state->constraintVec),x);
    algo_state.cnorm = (state->constraintVec)->norm();
    if (useDefaultInitPen_) {
      Step<Real>::getState()->searchSize
        = std::max(static_cast<Real>(1e-8),std::min(static_cast<Real>(10)*
            std::max(one,std::abs(fscale_*algo_state.value))
              /std::max(one,std::pow(cscale_*algo_state.cnorm,2)),
            static_cast<Real>(1e-2)*maxPenaltyParam_));
    }
    // Update evaluation counters
    algo_state.ncval += augLag.getNumberConstraintEvaluations();
    algo_state.nfval += augLag.getNumberFunctionEvaluations();
    algo_state.ngrad += augLag.getNumberGradientEvaluations();
    // Initialize intermediate stopping tolerances
    minPenaltyReciprocal_ = std::min(one/state->searchSize,minPenaltyLowerBound_);
    optTolerance_  = std::max<Real>(TOL*outerOptTolerance_,
                              optToleranceInitial_*std::pow(minPenaltyReciprocal_,optDecreaseExponent_));
    optTolerance_  = std::min<Real>(optTolerance_,TOL*algo_state.gnorm);
    feasTolerance_ = std::max<Real>(TOL*outerFeasTolerance_,
                              feasToleranceInitial_*std::pow(minPenaltyReciprocal_,feasDecreaseExponent_));
    if (verbosity_ > 0) {
      std::cout << std::endl;
      std::cout << "Augmented Lagrangian Initialize" << std::endl;
      std::cout << "Objective Scaling:  " << fscale_ << std::endl;
      std::cout << "Constraint Scaling: " << cscale_ << std::endl;
      std::cout << std::endl;
    }
  }

  /** \brief Compute step (equality constraint).
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, const Vector<Real> &l,
                Objective<Real> &obj, Constraint<Real> &con, 
                AlgorithmState<Real> &algo_state ) {
    compute(s,x,l,obj,con,*bnd_,algo_state);
  }

  /** \brief Compute step (equality and bound constraints).
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, const Vector<Real> &l,
                Objective<Real> &obj, Constraint<Real> &con, 
                BoundConstraint<Real> &bnd, AlgorithmState<Real> &algo_state ) {
    Real one(1);
    //AugmentedLagrangian<Real> &augLag
    //  = dynamic_cast<AugmentedLagrangian<Real>&>(obj);
    parlist_.sublist("Status Test").set("Gradient Tolerance",optTolerance_);
    parlist_.sublist("Status Test").set("Step Tolerance",1.e-6*optTolerance_);
    Ptr<Objective<Real>> penObj;
    if (subStep_ == "Bundle") {
      step_   = makePtr<BundleStep<Real>>(parlist_);
      status_ = makePtr<BundleStatusTest<Real>>(parlist_);
      penObj  = makePtrFromRef(obj);
    }
    else if (subStep_ == "Line Search") {
      step_   = makePtr<LineSearchStep<Real>>(parlist_);
      status_ = makePtr<StatusTest<Real>>(parlist_);
      penObj  = makePtrFromRef(obj);
    }
    else if (subStep_ == "Moreau-Yosida Penalty") {
      step_   = makePtr<MoreauYosidaPenaltyStep<Real>>(parlist_);
      status_ = makePtr<StatusTest<Real>>(parlist_);
      Ptr<Objective<Real>> raw_obj = makePtrFromRef(obj);
      penObj = ROL::makePtr<MoreauYosidaPenalty<Real>>(raw_obj,bnd_,x,parlist_);
    }
    else if (subStep_ == "Primal Dual Active Set") {
      step_   = makePtr<PrimalDualActiveSetStep<Real>>(parlist_);
      status_ = makePtr<StatusTest<Real>>(parlist_);
      penObj  = makePtrFromRef(obj);
    }
    else if (subStep_ == "Trust Region") {
      step_   = makePtr<TrustRegionStep<Real>>(parlist_);
      status_ = makePtr<StatusTest<Real>>(parlist_);
      penObj  = makePtrFromRef(obj);
    }
    else if (subStep_ == "Interior Point") {
      step_   = makePtr<InteriorPointStep<Real>>(parlist_);
      status_ = makePtr<StatusTest<Real>>(parlist_);
      Ptr<Objective<Real>> raw_obj = makePtrFromRef(obj);
      penObj = ROL::makePtr<InteriorPoint::PenalizedObjective<Real>>(raw_obj,bnd_,x,parlist_);
    }
    else {
      throw Exception::NotImplemented(">>> ROL::AugmentedLagrangianStep: Incompatible substep type!"); 
    }
    algo_ = makePtr<Algorithm<Real>>(step_,status_,false);
    //algo_ = ROL::makePtr<Algorithm<Real>>(subStep_,parlist_,false);
    x_->set(x);
    if ( bnd.isActivated() ) {
      //algo_->run(*x_,augLag,bnd,print_);
      algo_->run(*x_,*penObj,bnd,print_);
    }
    else {
      //algo_->run(*x_,augLag,print_);
      algo_->run(*x_,*penObj,print_);
    }
    s.set(*x_); s.axpy(-one,x);
    subproblemIter_ = (algo_->getState())->iter;
  }

  /** \brief Update step, if successful (equality constraint).
  */
  void update( Vector<Real> &x, Vector<Real> &l, const Vector<Real> &s,
               Objective<Real> &obj, Constraint<Real> &con,
               AlgorithmState<Real> &algo_state ) {
    update(x,l,s,obj,con,*bnd_,algo_state);
  }

  /** \brief Update step, if successful (equality and bound constraints).
  */
  void update( Vector<Real> &x, Vector<Real> &l, const Vector<Real> &s,
               Objective<Real> &obj, Constraint<Real> &con,
               BoundConstraint<Real> &bnd,
               AlgorithmState<Real> &algo_state ) {
    Real one(1), oem2(1.e-2);
    AugmentedLagrangian<Real> &augLag
      = dynamic_cast<AugmentedLagrangian<Real>&>(obj);
    ROL::Ptr<StepState<Real> > state = Step<Real>::getState();
    state->SPiter = subproblemIter_;
    // Update the step and store in state
    x.plus(s);
    algo_state.iterateVec->set(x);
    state->descentVec->set(s);
    algo_state.snorm = s.norm();
    algo_state.iter++;
    // Update objective function value
    obj.update(x);
    algo_state.value = augLag.getObjectiveValue(x);
    // Update constraint value
    augLag.getConstraintVec(*(state->constraintVec),x);
    algo_state.cnorm = (state->constraintVec)->norm();
    // Compute gradient of the augmented Lagrangian
    algo_state.gnorm  = computeGradient(*(state->gradientVec),x,state->searchSize,obj,bnd);
    algo_state.gnorm /= std::min(fscale_,cscale_);
    // Update evaluation counters
    algo_state.nfval += augLag.getNumberFunctionEvaluations();
    algo_state.ngrad += augLag.getNumberGradientEvaluations();
    algo_state.ncval += augLag.getNumberConstraintEvaluations();
    // Update objective function and constraints
    augLag.update(x,true,algo_state.iter);
    // Update multipliers
    minPenaltyReciprocal_ = std::min(one/state->searchSize,minPenaltyLowerBound_);
    if ( cscale_*algo_state.cnorm < feasTolerance_ ) {
      l.axpy(state->searchSize*cscale_,(state->constraintVec)->dual());
      if ( algo_->getState()->statusFlag == EXITSTATUS_CONVERGED ) {
        optTolerance_  = std::max(oem2*outerOptTolerance_,
                         optTolerance_*std::pow(minPenaltyReciprocal_,optIncreaseExponent_));
      }
      feasTolerance_ = std::max(oem2*outerFeasTolerance_,
                       feasTolerance_*std::pow(minPenaltyReciprocal_,feasIncreaseExponent_));
      // Update Algorithm State
      algo_state.snorm += state->searchSize*cscale_*algo_state.cnorm;
      algo_state.lagmultVec->set(l);
    }
    else {
      state->searchSize = std::min(penaltyUpdate_*state->searchSize,maxPenaltyParam_);
      optTolerance_     = std::max(oem2*outerOptTolerance_,
                          optToleranceInitial_*std::pow(minPenaltyReciprocal_,optDecreaseExponent_));
      feasTolerance_    = std::max(oem2*outerFeasTolerance_,
                          feasToleranceInitial_*std::pow(minPenaltyReciprocal_,feasDecreaseExponent_));
    }
    augLag.reset(l,state->searchSize);
  }

  /** \brief Print iterate header.
  */
  std::string printHeader( void ) const {
    std::stringstream hist;

    if(verbosity_>0) {
      hist << std::string(114,'-') << std::endl;
      hist << "Augmented Lagrangian status output definitions" << std::endl << std::endl;
      hist << "  iter    - Number of iterates (steps taken)"            << std::endl;
      hist << "  fval    - Objective function value"                    << std::endl;
      hist << "  cnorm   - Norm of the constraint violation"            << std::endl;
      hist << "  gLnorm  - Norm of the gradient of the Lagrangian"      << std::endl;
      hist << "  snorm   - Norm of the step"                            << std::endl;
      hist << "  penalty - Penalty parameter"                           << std::endl;
      hist << "  feasTol - Feasibility tolerance"                       << std::endl;
      hist << "  optTol  - Optimality tolerance"                        << std::endl;
      hist << "  #fval   - Number of times the objective was computed"  << std::endl;
      hist << "  #grad   - Number of times the gradient was computed"   << std::endl;
      hist << "  #cval   - Number of times the constraint was computed" << std::endl;
      hist << "  subIter - Number of iterations to solve subproblem"    << std::endl;
      hist << std::string(114,'-') << std::endl;
    }
    hist << "  ";
    hist << std::setw(6)  << std::left << "iter";
    hist << std::setw(15) << std::left << "fval";
    hist << std::setw(15) << std::left << "cnorm";
    hist << std::setw(15) << std::left << "gLnorm";
    hist << std::setw(15) << std::left << "snorm";
    hist << std::setw(10) << std::left << "penalty";
    hist << std::setw(10) << std::left << "feasTol";
    hist << std::setw(10) << std::left << "optTol";
    hist << std::setw(8)  << std::left << "#fval";
    hist << std::setw(8)  << std::left << "#grad";
    hist << std::setw(8)  << std::left << "#cval";
    hist << std::setw(8)  << std::left << "subIter";
    hist << std::endl;
    return hist.str();
  }

  /** \brief Print step name.
  */
  std::string printName( void ) const {
    std::stringstream hist;
    hist << std::endl << " Augmented Lagrangian Solver";
    hist << std::endl;
    hist << "Subproblem Solver: " << subStep_ << std::endl;
    return hist.str();
  }

  /** \brief Print iterate status.
  */
  std::string print( AlgorithmState<Real> &algo_state, bool pHeader = false ) const {
    std::stringstream hist;
    hist << std::scientific << std::setprecision(6);
    if ( algo_state.iter == 0 ) {
      hist << printName();
    }
    if ( pHeader ) {
      hist << printHeader();
    }
    if ( algo_state.iter == 0 ) {
      hist << "  ";
      hist << std::setw(6)  << std::left << algo_state.iter;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.cnorm;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << std::setw(15) << std::left << " ";
      hist << std::scientific << std::setprecision(2);
      hist << std::setw(10) << std::left << Step<Real>::getStepState()->searchSize;
      hist << std::setw(10) << std::left << std::max(feasTolerance_,outerFeasTolerance_);
      hist << std::setw(10) << std::left << std::max(optTolerance_,outerOptTolerance_);
      hist << std::endl;
    }
    else {
      hist << "  ";
      hist << std::setw(6)  << std::left << algo_state.iter;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.cnorm;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << std::setw(15) << std::left << algo_state.snorm;
      hist << std::scientific << std::setprecision(2);
      hist << std::setw(10) << std::left << Step<Real>::getStepState()->searchSize;
      hist << std::setw(10) << std::left << feasTolerance_;
      hist << std::setw(10) << std::left << optTolerance_;
      hist << std::scientific << std::setprecision(6);
      hist << std::setw(8) << std::left << algo_state.nfval;
      hist << std::setw(8) << std::left << algo_state.ngrad;
      hist << std::setw(8) << std::left << algo_state.ncval;
      hist << std::setw(8) << std::left << subproblemIter_;
      hist << std::endl;
    }
    return hist.str();
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

}; // class AugmentedLagrangianStep

} // namespace ROL

#endif
