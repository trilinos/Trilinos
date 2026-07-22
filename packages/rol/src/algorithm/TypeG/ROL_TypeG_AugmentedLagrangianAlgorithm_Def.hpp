// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEG_AUGMENTEDLAGRANGIANALGORITHM_DEF_H
#define ROL_TYPEG_AUGMENTEDLAGRANGIANALGORITHM_DEF_H

#include "ROL_TypeB_AlgorithmFactory.hpp"

namespace ROL {
namespace TypeG {

template<typename Real>
AugmentedLagrangianAlgorithm<Real>::AugmentedLagrangianAlgorithm( ParameterList &list, const Ptr<Secant<Real>> &secant )
  : TypeG::Algorithm<Real>::Algorithm(), secant_(secant), list_(list), subproblemIter_(0) {
  // Set status test
  status_->reset();
  status_->add(makePtr<ConstraintStatusTest<Real>>(list));

  Real one(1), p1(0.1), p9(0.9), ten(1.e1), oe8(1.e8), oem8(1.e-8);
  ParameterList& sublist = list.sublist("Step").sublist("Augmented Lagrangian");
  useDefaultInitPen_     = sublist.get("Use Default Initial Penalty Parameter", true);
  state_->searchSize     = sublist.get("Initial Penalty Parameter",             ten);
  // Multiplier update parameters
  scaleLagrangian_      = sublist.get("Use Scaled Augmented Lagrangian",          false);
  minPenaltyLowerBound_ = sublist.get("Penalty Parameter Reciprocal Lower Bound", p1);
  penaltyUpdate_        = sublist.get("Penalty Parameter Growth Factor",          ten);
  maxPenaltyParam_      = sublist.get("Maximum Penalty Parameter",                oe8);
  minPenaltyReciprocal_ = p1;
  // Optimality tolerance update
  optIncreaseExponent_ = sublist.get("Optimality Tolerance Update Exponent",   one);
  optDecreaseExponent_ = sublist.get("Optimality Tolerance Decrease Exponent", one);
  optToleranceInitial_ = sublist.get("Initial Optimality Tolerance",           one);
  // Feasibility tolerance update    
  feasIncreaseExponent_ = sublist.get("Feasibility Tolerance Update Exponent",   p1);
  feasDecreaseExponent_ = sublist.get("Feasibility Tolerance Decrease Exponent", p9);
  feasToleranceInitial_ = sublist.get("Initial Feasibility Tolerance",           one);
  // Subproblem information
  print_         = sublist.get("Print Intermediate Optimization History", false);
  maxit_         = sublist.get("Subproblem Iteration Limit",              1000);
  subStep_       = sublist.get("Subproblem Step Type",                    "Trust Region");
  HessianApprox_ = sublist.get("Level of Hessian Approximation",          0); 
  list_.sublist("Step").set("Type",subStep_);
  list_.sublist("Status Test").set("Iteration Limit",maxit_);
  list_.sublist("Status Test").set("Use Relative Tolerances",false);
  // Verbosity setting
  verbosity_          = list.sublist("General").get("Output Level", 0);
  printHeader_        = verbosity_ > 2;
  print_              = (verbosity_ > 2 ? true : print_);
  list_.sublist("General").set("Output Level",(print_ ? verbosity_ : 0));
  // Outer iteration tolerances
  outerFeasTolerance_ = list.sublist("Status Test").get("Constraint Tolerance",    oem8);
  outerOptTolerance_  = list.sublist("Status Test").get("Gradient Tolerance",      oem8);
  outerStepTolerance_ = list.sublist("Status Test").get("Step Tolerance",          oem8);
  useRelTol_          = list.sublist("Status Test").get("Use Relative Tolerances", false);
  // Scaling
  useDefaultScaling_  = sublist.get("Use Default Problem Scaling", true);
  fscale_             = sublist.get("Objective Scaling",           one);
  cscale_             = sublist.get("Constraint Scaling",          one);
}

template<typename Real>
void AugmentedLagrangianAlgorithm<Real>::initialize( Vector<Real>                       &x,
                                                     const Vector<Real>                 &g,
                                                     const Vector<Real>                 &l,
                                                     const Vector<Real>                 &c,
                                                     AugmentedLagrangianObjective<Real> &alobj,
                                                     BoundConstraint<Real>              &bnd,
                                                     Constraint<Real>                   &con,
                                                     std::ostream                       &outStream ) {
  hasPolyProj_ = true;
  if (proj_ == nullPtr) {
    proj_ = makePtr<PolyhedralProjection<Real>>(makePtrFromRef(bnd));
    hasPolyProj_ = false;
  }
  proj_->project(x,outStream);

  const Real one(1), TOL(1.e-2);
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  TypeG::Algorithm<Real>::initialize(x,g,l,c);

  // Initialize the algorithm state
  state_->nfval = 0;
  state_->ncval = 0;
  state_->ngrad = 0;

  // Compute objective value
  alobj.update(x,UpdateType::Initial,state_->iter);
  state_->value = alobj.getObjectiveValue(x,tol);
  alobj.gradient(*state_->gradientVec,x,tol);

  // Compute constraint violation
  state_->constraintVec->set(*alobj.getConstraintVec(x,tol));
  state_->cnorm = state_->constraintVec->norm();

  // Update evaluation counters
  state_->ncval += alobj.getNumberConstraintEvaluations();
  state_->nfval += alobj.getNumberFunctionEvaluations();
  state_->ngrad += alobj.getNumberGradientEvaluations();

  // Compute problem scaling
  if (useDefaultScaling_) {
    fscale_ = one/std::max(one,alobj.getObjectiveGradient(x,tol)->norm());
    try {
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
  alobj.setScaling(fscale_,cscale_);

  // Compute gradient of the lagrangian
  x.axpy(-one,state_->gradientVec->dual());
  proj_->project(x,outStream);
  x.axpy(-one/std::min(fscale_,cscale_),*state_->iterateVec);
  state_->gnorm = x.norm();
  x.set(*state_->iterateVec);

  // Compute initial penalty parameter
  if (useDefaultInitPen_) {
    const Real oem8(1e-8), oem2(1e-2), two(2), ten(10);
    state_->searchSize = std::max(oem8,
      std::min(ten*std::max(one,std::abs(fscale_*state_->value))
      / std::max(one,std::pow(cscale_*state_->cnorm,two)),oem2*maxPenaltyParam_));
  }
  // Initialize intermediate stopping tolerances
  if (useRelTol_) outerOptTolerance_ *= state_->gnorm;
  minPenaltyReciprocal_ = std::min(one/state_->searchSize,minPenaltyLowerBound_);
  optTolerance_  = std::max<Real>(TOL*outerOptTolerance_,
                            optToleranceInitial_*std::pow(minPenaltyReciprocal_,optDecreaseExponent_));
  optTolerance_  = std::min<Real>(optTolerance_,TOL*state_->gnorm);
  feasTolerance_ = std::max<Real>(TOL*outerFeasTolerance_,
                            feasToleranceInitial_*std::pow(minPenaltyReciprocal_,feasDecreaseExponent_));

  // Set data
  alobj.reset(l,state_->searchSize);

  if (verbosity_ > 1) {
    outStream << std::endl;
    outStream << "Augmented Lagrangian Initialize"            << std::endl;
    outStream << "Objective Scaling:  " << fscale_            << std::endl;
    outStream << "Constraint Scaling: " << cscale_            << std::endl;
    outStream << "Penalty Parameter:  " << state_->searchSize << std::endl;
    outStream << std::endl;
  }
}

template<typename Real>
void AugmentedLagrangianAlgorithm<Real>::run( Vector<Real>          &x,
                                              const Vector<Real>    &g,
                                              Objective<Real>       &obj,
                                              BoundConstraint<Real> &bnd,
                                              Constraint<Real>      &econ,
                                              Vector<Real>          &emul,
                                              const Vector<Real>    &eres,
                                              std::ostream          &outStream ) {
  const Real one(1), oem2(1e-2);
  Real tol(std::sqrt(ROL_EPSILON<Real>()));
  // Initialize augmented Lagrangian data
  AugmentedLagrangianObjective<Real> alobj(makePtrFromRef(obj),makePtrFromRef(econ),
                                           state_->searchSize,g,eres,emul,
                                           scaleLagrangian_,HessianApprox_);
  initialize(x,g,emul,eres,alobj,bnd,econ,outStream);
  Ptr<TypeB::Algorithm<Real>> algo;

  // Output
  if (verbosity_ > 0) writeOutput(outStream,true);

  while (status_->check(*state_)) {
    // Solve unconstrained augmented Lagrangian subproblem
    list_.sublist("Status Test").set("Gradient Tolerance",optTolerance_);
    list_.sublist("Status Test").set("Step Tolerance",1.e-6*optTolerance_);
    algo = TypeB::AlgorithmFactory<Real>(list_,secant_);
    if (hasPolyProj_) algo->run(x,g,alobj,bnd,*proj_->getLinearConstraint(),
                                *proj_->getMultiplier(),*proj_->getResidual(),
                                outStream);
    else              algo->run(x,g,alobj,bnd,outStream);
    subproblemIter_ = algo->getState()->iter;

    // Compute step
    state_->stepVec->set(x);
    state_->stepVec->axpy(-one,*state_->iterateVec);
    state_->snorm = state_->stepVec->norm();

    // Update iteration information
    state_->iter++;
    state_->iterateVec->set(x);
    state_->value = alobj.getObjectiveValue(x,tol);
    state_->constraintVec->set(*alobj.getConstraintVec(x,tol));
    state_->cnorm = state_->constraintVec->norm();
    alobj.gradient(*state_->gradientVec,x,tol);
    if (scaleLagrangian_) {
      state_->gradientVec->scale(state_->searchSize);
    }
    x.axpy(-one/std::min(fscale_,cscale_),state_->gradientVec->dual());
    proj_->project(x,outStream);
    x.axpy(-one,*state_->iterateVec);
    state_->gnorm = x.norm();
    x.set(*state_->iterateVec);
    //alobj.update(x,UpdateType::Accept,state_->iter);

    // Update evaluation counters
    state_->nfval += alobj.getNumberFunctionEvaluations();
    state_->ngrad += alobj.getNumberGradientEvaluations();
    state_->ncval += alobj.getNumberConstraintEvaluations();

    // Update multipliers
    if ( algo->getState()->statusFlag == EXITSTATUS_CONVERGED ) {
     minPenaltyReciprocal_ = std::min(one/state_->searchSize,minPenaltyLowerBound_);
     if ( cscale_*state_->cnorm < feasTolerance_ ) {
       emul.axpy(state_->searchSize*cscale_,state_->constraintVec->dual());
       optTolerance_  = std::max(oem2*outerOptTolerance_,
                        optTolerance_*std::pow(minPenaltyReciprocal_,optIncreaseExponent_));
       feasTolerance_ = std::max(oem2*outerFeasTolerance_,
                        feasTolerance_*std::pow(minPenaltyReciprocal_,feasIncreaseExponent_));
       // Update Algorithm State
       state_->snorm += state_->searchSize*cscale_*state_->cnorm;
       state_->lagmultVec->set(emul);
     }
     else {
       state_->searchSize = std::min(penaltyUpdate_*state_->searchSize,maxPenaltyParam_);
       optTolerance_      = std::max(oem2*outerOptTolerance_,
                            optToleranceInitial_*std::pow(minPenaltyReciprocal_,optDecreaseExponent_));
       feasTolerance_     = std::max(oem2*outerFeasTolerance_,
                            feasToleranceInitial_*std::pow(minPenaltyReciprocal_,feasDecreaseExponent_));
     }
     alobj.reset(emul,state_->searchSize);
    }

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,printHeader_);
  }
  if (verbosity_ > 0) TypeG::Algorithm<Real>::writeExitStatus(outStream);
}

template<typename Real>
void AugmentedLagrangianAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  if(verbosity_>1) {
    os << std::string(114,'-') << std::endl;
    os << "Augmented Lagrangian status output definitions" << std::endl << std::endl;
    os << "  iter    - Number of iterates (steps taken)"            << std::endl;
    os << "  fval    - Objective function value"                    << std::endl;
    os << "  cnorm   - Norm of the constraint violation"            << std::endl;
    os << "  gLnorm  - Norm of the gradient of the Lagrangian"      << std::endl;
    os << "  snorm   - Norm of the step"                            << std::endl;
    os << "  penalty - Penalty parameter"                           << std::endl;
    os << "  feasTol - Feasibility tolerance"                       << std::endl;
    os << "  optTol  - Optimality tolerance"                        << std::endl;
    os << "  #fval   - Number of times the objective was computed"  << std::endl;
    os << "  #grad   - Number of times the gradient was computed"   << std::endl;
    os << "  #cval   - Number of times the constraint was computed" << std::endl;
    os << "  subIter - Number of iterations to solve subproblem"    << std::endl;
    os << std::string(114,'-') << std::endl;
  }
  os << "  ";
  os << std::setw(6)  << std::left << "iter";
  os << std::setw(15) << std::left << "fval";
  os << std::setw(15) << std::left << "cnorm";
  os << std::setw(15) << std::left << "gLnorm";
  os << std::setw(15) << std::left << "snorm";
  os << std::setw(10) << std::left << "penalty";
  os << std::setw(10) << std::left << "feasTol";
  os << std::setw(10) << std::left << "optTol";
  os << std::setw(8)  << std::left << "#fval";
  os << std::setw(8)  << std::left << "#grad";
  os << std::setw(8)  << std::left << "#cval";
  os << std::setw(8)  << std::left << "subIter";
  os << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void AugmentedLagrangianAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::endl << "Augmented Lagrangian Solver (Type G, General Constraints)";
  os << std::endl;
  os << "Subproblem Solver: " << subStep_ << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void AugmentedLagrangianAlgorithm<Real>::writeOutput( std::ostream& os, const bool print_header ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::scientific << std::setprecision(6);
  if ( state_->iter == 0 ) writeName(os);
  if ( print_header )      writeHeader(os);
  if ( state_->iter == 0 ) {
    os << "  ";
    os << std::setw(6)  << std::left << state_->iter;
    os << std::setw(15) << std::left << state_->value;
    os << std::setw(15) << std::left << state_->cnorm;
    os << std::setw(15) << std::left << state_->gnorm;
    os << std::setw(15) << std::left << "---";
    os << std::scientific << std::setprecision(2);
    os << std::setw(10) << std::left << state_->searchSize;
    os << std::setw(10) << std::left << std::max(feasTolerance_,outerFeasTolerance_);
    os << std::setw(10) << std::left << std::max(optTolerance_,outerOptTolerance_);
    os << std::scientific << std::setprecision(6);
    os << std::setw(8) << std::left << state_->nfval;
    os << std::setw(8) << std::left << state_->ngrad;
    os << std::setw(8) << std::left << state_->ncval;
    os << std::setw(8) << std::left << "---";
    os << std::endl;
  }
  else {
    os << "  ";
    os << std::setw(6)  << std::left << state_->iter;
    os << std::setw(15) << std::left << state_->value;
    os << std::setw(15) << std::left << state_->cnorm;
    os << std::setw(15) << std::left << state_->gnorm;
    os << std::setw(15) << std::left << state_->snorm;
    os << std::scientific << std::setprecision(2);
    os << std::setw(10) << std::left << state_->searchSize;
    os << std::setw(10) << std::left << feasTolerance_;
    os << std::setw(10) << std::left << optTolerance_;
    os << std::scientific << std::setprecision(6);
    os << std::setw(8) << std::left << state_->nfval;
    os << std::setw(8) << std::left << state_->ngrad;
    os << std::setw(8) << std::left << state_->ncval;
    os << std::setw(8) << std::left << subproblemIter_;
    os << std::endl;
  }
  os.flags(osFlags);
}

} // namespace TypeG
} // namespace ROL

#endif
