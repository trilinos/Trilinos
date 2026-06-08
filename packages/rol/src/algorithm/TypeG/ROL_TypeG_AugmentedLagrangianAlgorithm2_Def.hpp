// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEG_AUGMENTEDLAGRANGIANALGORITHM2_DEF_H
#define ROL_TYPEG_AUGMENTEDLAGRANGIANALGORITHM2_DEF_H

#include <cmath>

#include "ROL_AugmentedLagrangianObjective2.hpp"


namespace ROL {

template<class Real> class Solver;

namespace TypeG {

template <typename Real> class Algorithm;

template <typename Real>
inline Ptr<TypeG::Algorithm<Real>> AlgorithmFactory(
    ParameterList &parlist,
    const Ptr<Secant<Real>> &secant);

template<typename Real>
AugmentedLagrangianAlgorithm2<Real>::AugmentedLagrangianAlgorithm2( ParameterList &list, const Ptr<Secant<Real>> &secant )
  : TypeG::Algorithm<Real>::Algorithm(), secant_(secant), list_(list), subproblemIter_(0) {
  // Set status test
  status_->reset();
  status_->add(makePtr<ConstraintStatusTest<Real>>(list));

  Real one(1), p1(0.1), p9(0.9), ten(1.e1), oe8(1.e8), oem8(1.e-8);
  ParameterList& sublist = list.sublist("Step").sublist("Augmented Lagrangian");

  // Penalty parameter quantities
  useDefaultInitPen_     = sublist.get("Use Default Initial Penalty Parameter", true);
  state_->searchSize     = sublist.get("Initial Penalty Parameter",             ten);
  penaltyUpdate_         = sublist.get("Penalty Parameter Growth Factor",       ten);
  maxPenaltyParam_       = sublist.get("Maximum Penalty Parameter",             oe8);
  minPenaltyReciprocal_  = p1; 

  // Dual feasibility parameters
  alphat_ = sublist.get("Feasibility Tolerance Update Exponent",    p1);
  betat_  = sublist.get("Feasibility Tolerance Decrease Exponent",  p9);
  theta_  = sublist.get("Penalty Parameter Reciprocal Lower Bound", p1); 
  tau0_   = sublist.get("Initial Dual Feasibility Tolerance",       one);

  // Subproblem information
  useDefaultInitTol_  = sublist.get("Use Default Initial Subproblem Tolerances", false);
  epsilon_            = sublist.get("Initial Optimality Tolerance",              1e-4);
  delta_              = sublist.get("Initial Feasibility Tolerance",             1e-4);
  maxit_              = sublist.get("Subproblem Iteration Limit",                1000);
  bool print          = sublist.get("Print Intermediate Optimization History",   false);
  subStep_            = sublist.get("Subproblem Step Type",             "Trust Region");
  list_.sublist("Step").set("Type",subStep_);
  list_.sublist("Status Test").set("Iteration Limit", maxit_);
  list_.sublist("Status Test").set("Use Relative Tolerances",false);

  // Optimality parameters
  optIncreaseExponent_  = sublist.get("Optimality Tolerance Update Exponent",   one);
  optDecreaseExponent_  = sublist.get("Optimality Tolerance Decrease Exponent", one);

  // Verbosity 
  verbosity_    = list.sublist("General").get("Output Level", 0);
  printHeader_  = verbosity_ > 0;
  print         = (verbosity_ > 2 ? true : print);
  list_.sublist("General").set("Output Level",(print ? verbosity_ - 2 : 0));

  // Outer iteration tolerances
  outerFeasTolerance_ = list.sublist("Status Test").get("Constraint Tolerance",    oem8);
  outerOptTolerance_  = list.sublist("Status Test").get("Gradient Tolerance",      oem8);
  outerStepTolerance_ = list.sublist("Status Test").get("Step Tolerance",          oem8);
  useRelTol_          = list.sublist("Status Test").get("Use Relative Tolerances", false);

  // Augmented Lagrangian parameters
  useDefaultScaling_  = sublist.get("Use Default Problem Scaling",     true);
  scaleLagrangian_    = sublist.get("Use Scaled Augmented Lagrangian", false);
  fscale_             = sublist.get("Objective Scaling",               one);
  cscale_             = sublist.get("Constraint Scaling",              one);
  hessianApprox_      = sublist.get("Level of Hessian Approximation",  0);

  // Scaling

}

template<typename Real>
void AugmentedLagrangianAlgorithm2<Real>::initialize( Vector<Real>                        &x,
                                                     const Vector<Real>                   &g,
                                                     AugmentedLagrangianObjective2<Real>  &alobj,
                                                     std::ostream                         &outStream ) {

  if (useRelTol_) {
    outStream << "Warning: \"Use Relative Tolerances\" parameter is unsupported!" << std::endl;
    useRelTol_ = false;
  }

  const Real one(1), TOL(1.e-2);
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  // > TypeG::Algorithm<Real>::initialize(x,g,l,c);
  if (state_->iterateVec == nullPtr) {
    state_->iterateVec = x.clone();
  }
  state_->iterateVec->set(x);
  // > if (state_->lagmultVec == nullPtr) {
  // >   state_->lagmultVec = mul.clone();
  // > }
  // > state_->lagmultVec->set(mul);
  if (state_->stepVec == nullPtr) {
    state_->stepVec = x.clone();
  }
  state_->stepVec->zero();
  if (state_->gradientVec == nullPtr) {
    state_->gradientVec = g.clone();
  }
  state_->gradientVec->set(g);
  // > if (state_->constraintVec == nullPtr) {
  // >   state_->constraintVec = c.clone();
  // > }
  // > state_->constraintVec->zero();
  if (state_->minIterVec == nullPtr) {
    state_->minIterVec = x.clone();
  }
  state_->minIterVec->set(x);
  state_->minIter = state_->iter;
  state_->minValue = state_->value;

  // Initialize the algorithm state
  state_->nfval = 0;
  state_->ncval = 0;
  state_->ngrad = 0;

  // Compute objective value
  alobj.update(x,UpdateType::Initial,state_->iter);
  state_->value = alobj.getObjectiveValue(x,tol);
  // > alobj.gradient(*state_->gradientVec,x,tol);

  // ========================================================================
  // lagmultVec, constraintVec, and gradientVec are unused state_ elements.
  // ========================================================================

  if (useDefaultScaling_) {
    fscale_ = one/std::max(one,alobj.getObjectiveGradient(x,tol)->norm());
    alobj.setScaling(fscale_);
  }

  unsigned numberPenalties = alobj.getNumberConstraints();

  feasibilities_.resize(numberPenalties);

  // Compute constraint violation
  for (unsigned i = 0; i < numberPenalties; ++i)
    feasibilities_[i] = alobj.feasibility(x,tol,i);
  state_->cnorm = 0;
  if (!feasibilities_.empty())
    state_->cnorm = *std::max_element(feasibilities_.begin(),feasibilities_.end());

  // Update evaluation counters
  state_->ncval += alobj.getNumberConstraintEvaluations();
  state_->nfval += alobj.getNumberFunctionEvaluations();
  state_->ngrad += alobj.getNumberGradientEvaluations();

  // Compute problem scaling
  std::vector<Real> constraintScalings(numberPenalties);
  Ptr<Vector<Real>> constraintVec;
  // > TODO
  // >if (useDefaultScaling_) {
  // >  Ptr<Vector<Real>> ji = x.clone();
  // >  // Helper function
  // >  auto getConstraintScale = [&] (const Vector<Real> &c) -> Real {
  // >    Real cscale(1);
  // >    try {
  // >      Real maxji(0), normji(0);
  // >      for (int i = 0; i < c.dimension(); ++i) {
  // >        c.applyAdjointJacobian(*ji,*c.basis(i),x,tol);
  // >        normji = ji->norm();
  // >        maxji  = std::max(normji,maxji);
  // >      }
  // >      cscale = one/std::max(one,maxji);
  // >    }
  // >    catch (std::exception &e) {}
  // >    return cscale;
  // >  };
  // >  for (unsigned k = 0; k < numberPenalties; ++k) {
  // >    constraintVec = alobj.getConstraintVec(x,tol,k);
  // >    constraintScalings[k] = getConstraintScale(*constraintVec);
  // >    alobj.setScaling(constraintScalings[k],k);
  // >  }
  // >  if (!constraintScalings.empty())
  // >    cscale_ = *std::max_element(constraintScalings.begin(),constraintScalings.end());
  // >}
  // >else {
  if (!useDefaultScaling_) {
    for (unsigned i = 0; i < numberPenalties; ++i) {
      constraintScalings[i] = cscale_;
      alobj.setScaling(constraintScalings[i],i);
    }
  }

  // > // Compute gradient of the lagrangian
  // > x.axpy(-one,state_->gradientVec->dual());
  // > proj_->project(x,outStream);
  // > x.axpy(-one/std::min(fscale_,cscale_),*state_->iterateVec);
  // > state_->gnorm = x.norm();
  // > x.set(*state_->iterateVec);
  // > TODO: problem.getSubproblemStationarityMeasure(*state_->iterateVec,x); which is used below to set tolerances

  Real temp;

  // Initial penalty parameters
  if (useDefaultInitPen_) {
    const Real oem8(1e-8), oem2(1e-2);
    state_->searchSize = 0;
    for (unsigned i = 0; i < numberPenalties; ++i) {
      temp = alobj.dualNorm(x,tol,i);
      if (temp <= oem8) temp = 1;
      temp = std::max(oem8, oem2*std::abs(fscale_*state_->value)/temp);
      alobj.setPenaltyParameter(temp,i);
      state_->searchSize = std::max(state_->searchSize,temp);
    }
  }
  else {
    for (unsigned i = 0; i < numberPenalties; ++i) {
      Real max_penalty = 0;
      temp = list_.sublist("Step").sublist("Augmented Lagrangian").sublist(group_names_[i]).get("Initial Penalty Parameter",state_->searchSize);
      alobj.setPenaltyParameter(temp,i);
      max_penalty = std::max(max_penalty,temp);
      state_->searchSize = max_penalty;
    }
  }

  // Penalty updates
  for (unsigned i = 0; i < numberPenalties; ++i) {
    temp = list_.sublist("Step").sublist("Augmented Lagrangian").sublist(group_names_[i]).get("Penalty Parameter Growth Factor",penaltyUpdate_);
    penalty_growthf_.push_back(temp);
  }

  // Subproblem stopping
  if (useDefaultInitTol_) {
    temp = 0;
    for (unsigned i = 0; i < numberPenalties; ++i) 
      temp += alobj.getPenaltyParameter(i);
    temp = std::min(minPenaltyReciprocal_,1./temp);
    epsilon_ = std::max(TOL*outerOptTolerance_, epsilon_*std::pow(temp,optDecreaseExponent_));
    delta_   = std::max(TOL*outerFeasTolerance_,delta_  *std::pow(temp,optDecreaseExponent_));
  }
  
  alobj.reset();

  if (verbosity_ > 1) {
    outStream << std::endl;
    outStream << "Augmented Lagrangian Initialize"                    << std::endl;
    outStream << "Objective Scaling:  "         << fscale_            << std::endl;
    outStream << "Maximum Constraint Scaling: " << cscale_            << std::endl;
    outStream << "Maximum Penalty Parameter:  " << state_->searchSize << std::endl;
    outStream << std::endl;
  }
}

template<typename Real>
void AugmentedLagrangianAlgorithm2<Real>::run( Problem<Real> &problem,
                                               std::ostream  &outStream ) {

  // ========================================================================
  // STEP 1: Get augmented Lagrangian objective and subproblem
  // ========================================================================

  problem.finalize(false,verbosity_>1,outStream,false);
  Vector<Real> &x = *problem.getPrimalOptimizationVector();

  auto makePenalty = [&] (Ptr<ConstraintData<Real>>& constraint_data)
                           -> Ptr<AugmentedLagrangianPenalty<Real>> {
    Ptr<AugmentedLagrangianPenalty<Real>> penalty = makePtr<AugmentedLagrangianPenalty<Real>>(
                                                      constraint_data->constraint,
                                                      constraint_data->projection,
                                                      Real(1.0),
                                                      x.dual(),
                                                      *constraint_data->residual,
                                                      *constraint_data->multiplier,
                                                      hessianApprox_);
    penalty->setMultiplier(*constraint_data->multiplier);
    return penalty;
  };

  group_names_ = problem.getGroupNames();

  Ptr<BoundConstraint<Real>> bound_constraint = problem.getBoundConstraint();
  bool has_bound_constraint = (bound_constraint != nullPtr);
  std::vector<std::string> ungrouped_linear_equality_constraints = problem.getUngroupedLinearEqualityConstraintNames();
  std::vector<std::string> ungrouped_equality_constraints = problem.getUngroupedEqualityConstraintNames();

  Ptr<ConstraintData<Real>> constraint_data;

  // Construct augmented Lagrangian objective

  std::vector<Ptr<AugmentedLagrangianPenalty<Real>>> penalties;
  for (const std::string& name : group_names_) {
    constraint_data = problem.getGroupConstraintData(name);
    penalties.push_back(makePenalty(constraint_data));
  }
  // Penalizes the bound constraint
  if (ungrouped_equality_constraints.size() > 0 && has_bound_constraint) {
    Ptr<Projection<Real>> projection = makePtr<PolyhedralProjection<Real>>(bound_constraint);
    Ptr<Vector<Real>> b = x.clone();
    b->zero();
    Ptr<LinearOperator<Real>> A = makePtr<IdentityOperator<Real>>();
    Ptr<Constraint<Real>> constraint = makePtr<LinearConstraint<Real>>(A,b);
    constraint_data = makePtr<ConstraintData<Real>>(constraint,x.dual().clone(),makePtrFromRef(x),nullPtr,projection);
    penalties.push_back(makePenalty(constraint_data));
    group_names_.push_back("Bounds");
  }
  Ptr<Objective<Real>> objective = problem.getObjective();
  Ptr<AugmentedLagrangianObjective2<Real>> alobj = makePtr<AugmentedLagrangianObjective2<Real>>(objective,penalties,x.dual(),false);

  Ptr<Problem<Real>> subproblem = makePtr<Problem<Real>>(alobj,makePtrFromRef(x));

  if (ungrouped_equality_constraints.size() == 0 && has_bound_constraint)
    subproblem->addBoundConstraint(bound_constraint);

  for (const auto& name : ungrouped_linear_equality_constraints) {
    constraint_data = problem.getConstraintData(name);
    subproblem->addLinearConstraint(name,constraint_data->constraint,constraint_data->multiplier,constraint_data->residual);
  }

  for (const auto& name : ungrouped_equality_constraints) {
    constraint_data = problem.getConstraintData(name);
    subproblem->addConstraint(name,constraint_data->constraint,constraint_data->multiplier,constraint_data->residual);
  }

  if (verbosity_ > 1) {
    outStream << std::endl << "Subproblem Finalize:" << std::endl;
    subproblem->finalize(false,verbosity_>0,outStream,true);
  }

  // ========================================================================
  // STEP 2: Configure
  // ========================================================================

  initialize(x,x.dual(),*alobj,outStream);

  // Constants
  const Real one(1);
  Real tol(std::sqrt(ROL_EPSILON<Real>()));

  // Additional parameters
  unsigned maxSubproblemFails = 2;
  int k0 = 100;
  Real nu = 1.e6;
  Real gamma = 0.49;

  state_->gnorm = ROL_INF<Real>();

  // Output
  if (verbosity_ > 0) writeOutput(outStream,true);

  // Data

  unsigned numberPenalties = alobj->getNumberConstraints();

  std::vector<Real> dualResiduals(numberPenalties);

  std::vector<Real> dualTolerances(numberPenalties); // \tau
  for (unsigned i = 0; i < numberPenalties; ++i)
    dualTolerances[i] = tau0_*std::pow(theta_,alphat_);

  EExitStatus statusFlag;
  bool isSubproblemConverged = false;
  Real penaltyParameter;
  bool isForcedUpdate = false;

  Real theta;

  // ========================================================================
  // STEP 3: Run algorithm
  // ========================================================================

  while (status_->check(*state_)) {
    // Solve augmented Lagrangian subproblem
    list_.sublist("Status Test").set("Gradient Tolerance",epsilon_);
    list_.sublist("Status Test").set("Constraint Tolerance",delta_);
    //list_.sublist("Status Test").set("Step Tolerance",1.e-6*optTolerance_);
    list_.sublist("Status Test").set("Step Tolerance",1.e-14);
    Solver<Real> solver(subproblem,list_);
    isSubproblemConverged = false;
    for (unsigned i = 0; i < maxSubproblemFails; ++i) {
      solver.reset();
      alobj->update(x,UpdateType::Initial,state_->iter);
      solver.solve(outStream);
      statusFlag = solver.getAlgorithmState()->statusFlag;
      isSubproblemConverged = (statusFlag == EXITSTATUS_CONVERGED) || (statusFlag == EXITSTATUS_USERDEFINED);
      if (isSubproblemConverged)
        break;
    }

    if (!isSubproblemConverged) {
      if (verbosity_ > 0)
        outStream << "Warning: Augmented Lagrangian subproblem failed to converge!" << std::endl;
    }

    subproblemIter_ = solver.getAlgorithmState()->iter;

    // Compute step
    state_->stepVec->set(x);
    state_->stepVec->axpy(-one,*state_->iterateVec);
    state_->snorm = state_->stepVec->norm();

    // Update iteration information
    state_->iter++;
    state_->iterateVec->set(x);
    state_->value = alobj->getObjectiveValue(x,tol);
    for (unsigned i = 0; i < numberPenalties; ++i)
      feasibilities_[i] = alobj->feasibility(x,tol,i);
    if (!feasibilities_.empty())
      state_->cnorm = *std::max_element(feasibilities_.begin(),feasibilities_.end());
    state_->gnorm = solver.getAlgorithmState()->gnorm;
    //alobj.update(x,UpdateType::Accept,state_->iter);

    // Update evaluation counters
    state_->nfval += alobj->getNumberFunctionEvaluations();
    state_->ngrad += alobj->getNumberGradientEvaluations();
    state_->ncval += alobj->getNumberConstraintEvaluations();

    // ========================================================================
    // ALESQP Updates (Algorithm 4.1)
    // ========================================================================

    for (unsigned i = 0; i < numberPenalties; ++i)
      dualResiduals[i] = alobj->dualResidual(x,tol,i);

    // Forced updates for large iterations (line 10)
    isForcedUpdate = false;
    if (state_->iter > k0) {
      for (unsigned i = 0; i < numberPenalties; ++i) {
        penaltyParameter = alobj->getPenaltyParameter(i);
        if (alobj->getScaling(i)*dualResiduals[i] > penaltyParameter*dualTolerances[i])
          isForcedUpdate = true;
      }
    }

    // Update (line 13)
    isUpdated_ = false;
    for (unsigned i = 0; i < numberPenalties; ++i) {
      penaltyParameter = alobj->getPenaltyParameter(i);
      // std::cout << "Dual Residuals " << i << ": " << dualResiduals[i] << " Scaling: " << alobj->getScaling(i) << " Tolerance: " << dualTolerances[i] << std::endl;
      if (isForcedUpdate || alobj->getScaling(i)*dualResiduals[i] > penaltyParameter*dualTolerances[i]) {
        penaltyParameter  *= penalty_growthf_[i];
        penaltyParameter   = std::min(penaltyParameter,maxPenaltyParam_);
        state_->searchSize = std::max(state_->searchSize,penaltyParameter);
        theta              = std::min(one/penaltyParameter,theta_);
        dualTolerances[i]  = tau0_*std::pow(theta,alphat_);
        // dualTolerances[i]  = std::max(dualTolerances[i],oem2*outerFeasTolerance_);   // ROL convention
        alobj->setPenaltyParameter(penaltyParameter,i);
        isUpdated_ = true;
      }
      else {
        theta              = std::min(one/penaltyParameter,theta_);
        dualTolerances[i] *= std::pow(theta,betat_);
        // dualTolerances[i]  = std::max(dualTolerances[i],oem2*outerFeasTolerance_); // ROL convention
      }
      if (alobj->getScaling(i)*dualResiduals[i] <= nu*std::pow(penaltyParameter,gamma)) {
        alobj->updateMultiplier(x,tol,i);
        // outStream << "multiplier update" << std::endl;
      }
    }
    if (isUpdated_) {
      epsilon_ = 0.9*epsilon_;
      delta_   = 0.9*delta_;
    }
    else {
      epsilon_ = 0.25*epsilon_;
      delta_   = 0.25*delta_;
    }

    alobj->reset();

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,false);
  }
  if (verbosity_ > 0) {
    outStream << std::endl;
    TypeG::Algorithm<Real>::writeExitStatus(outStream);
  }
}


template<typename Real>
void AugmentedLagrangianAlgorithm2<Real>::run( Vector<Real>         &x,
                                              const Vector<Real>    &g,
                                              Objective<Real>       &obj,
                                              BoundConstraint<Real> &bnd,
                                              Constraint<Real>      &econ,
                                              Vector<Real>          &emul,
                                              const Vector<Real>    &eres,
                                              std::ostream          &outStream ) {
  throw Exception::NotImplemented(">>> ROL::TypeG""AugmentedLagrangianAlgorithm2::run: Overload Not Implemented!");
}

template<typename Real>
void AugmentedLagrangianAlgorithm2<Real>::writeHeader( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  if(verbosity_ > 1 && state_->iter == 0) {
    os << std::string(114,'-') << std::endl;
    os << "Augmented Lagrangian status output definitions" << std::endl << std::endl;
    os << "  iter       - Number of iterates (steps taken)"            << std::endl;
    os << "  fval       - Objective function value"                    << std::endl;
    os << "  max infeas - Largest infeasibility"                       << std::endl;
    os << "  optimality - Norm of the subproblem optimality measure"   << std::endl;
    os << "  snorm      - Norm of the step"                            << std::endl;
    os << "  max pen    - Largest penalty parameter"                   << std::endl;
    os << "  feasTol    - Subproblem Feasibility tolerance"            << std::endl;
    os << "  optTol     - Subproblem Optimality tolerance"             << std::endl;
    os << "  #fval      - Number of times the objective was computed"  << std::endl;
    os << "  #grad      - Number of times the gradient was computed"   << std::endl;
    os << "  #cval      - Number of times the constraint was computed" << std::endl;
    os << "  subIter    - Number of iterations to solve subproblem"    << std::endl;
    os << std::string(114,'-') << std::endl;
  }
  if (verbosity_ > 1) os << std::endl;
  os << "  ";
  os << std::setw(6)  << std::left << "iter";
  os << std::setw(15) << std::left << "fval";
  os << std::setw(15) << std::left << "max infeas";
  os << std::setw(15) << std::left << "optimality";
  os << std::setw(15) << std::left << "snorm";
  os << std::setw(10) << std::left << "max pen";
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
void AugmentedLagrangianAlgorithm2<Real>::writeName( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::endl << "Augmented Lagrangian Solver (Type G, General Constraints)";
  os << std::endl;
  os << "Subproblem Solver: " << subStep_ << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void AugmentedLagrangianAlgorithm2<Real>::writeOutput( std::ostream& os, const bool print_header ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::scientific << std::setprecision(6);
  if ( state_->iter == 0 ) writeName(os);
  if (  print_header || verbosity_ > 1 ) writeHeader(os);
  if ( state_->iter == 0 ) {
    os << "  ";
    os << std::setw(6)  << std::left << state_->iter;
    os << std::setw(15) << std::left << state_->value;
    os << std::setw(15) << std::left << state_->cnorm;
    os << std::setw(15) << std::left << "---";
    os << std::setw(15) << std::left << "---";
    os << std::scientific << std::setprecision(2);
    os << std::setw(10) << std::left << state_->searchSize;
    os << std::setw(10) << std::left << delta_;
    os << std::setw(10) << std::left << epsilon_;
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
    os << std::setw(10) << std::left << delta_;
    os << std::setw(10) << std::left << epsilon_;
    os << std::scientific << std::setprecision(6);
    os << std::setw(8) << std::left << state_->nfval;
    os << std::setw(8) << std::left << state_->ngrad;
    os << std::setw(8) << std::left << state_->ncval;
    os << std::setw(8) << std::left << subproblemIter_;
    os << std::endl;
  }
  if (verbosity_ > 1) {
    os << std::endl;
    os << "    Penalty Update?  " << (isUpdated_ ? "Yes" : "No") << "  |  Feasibilities:  ";
    os << std::scientific << std::setprecision(6);
    for (unsigned i = 0; i < feasibilities_.size(); ++i) {
      std::string name = group_names_[i];
      os << name << ": "  << std::setw(8) << std::left << feasibilities_[i] << "  ";
    }
    os << std::endl;
  }
  os.flags(osFlags);
}

} // namespace TypeG
} // namespace ROL

#include "ROL_Solver.hpp"

#endif
