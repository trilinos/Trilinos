// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_TYPEE_AUGMENTEDLAGRANGIANALGORITHM_DEF_H
#define ROL_TYPEE_AUGMENTEDLAGRANGIANALGORITHM_DEF_H

#include "ROL_TypeU_AlgorithmFactory.hpp"

namespace ROL {
namespace TypeE {

template<typename Real>
AugmentedLagrangianAlgorithm<Real>::AugmentedLagrangianAlgorithm( ParameterList &list )
  : TypeE::Algorithm<Real>::Algorithm(), list_(list), subproblemIter_(0) {
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
                                                     Constraint<Real>                   &con,
                                                     std::ostream                       &outStream ) {
  const Real one(1), TOL(1.e-2);
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  TypeE::Algorithm<Real>::initialize(x,g,l,c);

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
  state_->gnorm = state_->gradientVec->norm()/std::min(fscale_,cscale_);

  // Compute initial penalty parameter
  if (useRelTol_) outerOptTolerance_  *= state_->gnorm;
  if (useDefaultInitPen_) {
    const Real oem8(1e-8), oem2(1e-2), two(2), ten(10);
    state_->searchSize = std::max(oem8,
      std::min(ten*std::max(one,std::abs(fscale_*state_->value))
      / std::max(one,std::pow(cscale_*state_->cnorm,two)),oem2*maxPenaltyParam_));
  }
  // Initialize intermediate stopping tolerances
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
    outStream << "Augmented Lagrangian Initialize" << std::endl;
    outStream << "Objective Scaling:  " << fscale_ << std::endl;
    outStream << "Constraint Scaling: " << cscale_ << std::endl;
    outStream << std::endl;
  }
}

template<typename Real>
void AugmentedLagrangianAlgorithm<Real>::run( Vector<Real>       &x,
                                              const Vector<Real> &g,
                                              Objective<Real>    &obj,
                                              Constraint<Real>   &econ,
                                              Vector<Real>       &emul,
                                              const Vector<Real> &eres,
                                              std::ostream       &outStream ) {
  const Real one(1), oem2(1e-2);
  Real tol(std::sqrt(ROL_EPSILON<Real>()));
  // Initialize augmented Lagrangian data
  AugmentedLagrangianObjective<Real> alobj(makePtrFromRef(obj),makePtrFromRef(econ),
                                           state_->searchSize,g,eres,emul,
                                           scaleLagrangian_,HessianApprox_);
  initialize(x,g,emul,eres,alobj,econ,outStream);
  Ptr<TypeU::Algorithm<Real>> algo;

  // Output
  if (verbosity_ > 0) writeOutput(outStream,true);

  while (status_->check(*state_)) {
    // Solve unconstrained augmented Lagrangian subproblem
    list_.sublist("Status Test").set("Gradient Tolerance",optTolerance_);
    list_.sublist("Status Test").set("Step Tolerance",1.e-6*optTolerance_);
    algo = TypeU::AlgorithmFactory<Real>(list_);
    algo->run(x,g,alobj,outStream);
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
    state_->gnorm = state_->gradientVec->norm()/std::min(fscale_,cscale_);
    //alobj.update(x,UpdateType::Accept,state_->iter);

    // Update evaluation counters
    state_->nfval += alobj.getNumberFunctionEvaluations();
    state_->ngrad += alobj.getNumberGradientEvaluations();
    state_->ncval += alobj.getNumberConstraintEvaluations();

    // Update multipliers
    minPenaltyReciprocal_ = std::min(one/state_->searchSize,minPenaltyLowerBound_);
    if ( cscale_*state_->cnorm < feasTolerance_ ) {
      emul.axpy(state_->searchSize*cscale_,state_->constraintVec->dual());
      if ( algo->getState()->statusFlag == EXITSTATUS_CONVERGED ) {
        optTolerance_  = std::max(oem2*outerOptTolerance_,
                         optTolerance_*std::pow(minPenaltyReciprocal_,optIncreaseExponent_));
      }
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

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,printHeader_);
  }
  emul.scale(cscale_);
  if (verbosity_ > 0) TypeE::Algorithm<Real>::writeExitStatus(outStream);
}

template<typename Real>
void AugmentedLagrangianAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::stringstream hist;
  if(verbosity_>1) {
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
  os << hist.str();
}

template<typename Real>
void AugmentedLagrangianAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::stringstream hist;
  hist << std::endl << "Augmented Lagrangian Solver (Type E, Equality Constraints)";
  hist << std::endl;
  hist << "Subproblem Solver: " << subStep_ << std::endl;
  os << hist.str();
}

template<typename Real>
void AugmentedLagrangianAlgorithm<Real>::writeOutput( std::ostream& os, const bool print_header ) const {
  std::stringstream hist;
  hist << std::scientific << std::setprecision(6);
  if ( state_->iter == 0 ) writeName(os);
  if ( print_header )      writeHeader(os);
  if ( state_->iter == 0 ) {
    hist << "  ";
    hist << std::setw(6)  << std::left << state_->iter;
    hist << std::setw(15) << std::left << state_->value;
    hist << std::setw(15) << std::left << state_->cnorm;
    hist << std::setw(15) << std::left << state_->gnorm;
    hist << std::setw(15) << std::left << "---";
    hist << std::scientific << std::setprecision(2);
    hist << std::setw(10) << std::left << state_->searchSize;
    hist << std::setw(10) << std::left << std::max(feasTolerance_,outerFeasTolerance_);
    hist << std::setw(10) << std::left << std::max(optTolerance_,outerOptTolerance_);
    hist << std::scientific << std::setprecision(6);
    hist << std::setw(8) << std::left << state_->nfval;
    hist << std::setw(8) << std::left << state_->ngrad;
    hist << std::setw(8) << std::left << state_->ncval;
    hist << std::setw(8) << std::left << "---";
    hist << std::endl;
  }
  else {
    hist << "  ";
    hist << std::setw(6)  << std::left << state_->iter;
    hist << std::setw(15) << std::left << state_->value;
    hist << std::setw(15) << std::left << state_->cnorm;
    hist << std::setw(15) << std::left << state_->gnorm;
    hist << std::setw(15) << std::left << state_->snorm;
    hist << std::scientific << std::setprecision(2);
    hist << std::setw(10) << std::left << state_->searchSize;
    hist << std::setw(10) << std::left << feasTolerance_;
    hist << std::setw(10) << std::left << optTolerance_;
    hist << std::scientific << std::setprecision(6);
    hist << std::setw(8) << std::left << state_->nfval;
    hist << std::setw(8) << std::left << state_->ngrad;
    hist << std::setw(8) << std::left << state_->ncval;
    hist << std::setw(8) << std::left << subproblemIter_;
    hist << std::endl;
  }
  os << hist.str();
}

} // namespace TypeE
} // namespace ROL

#endif
