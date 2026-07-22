// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEG_STABILIZEDLCLALGORITHM_DEF_H
#define ROL_TYPEG_STABILIZEDLCLALGORITHM_DEF_H

#include "ROL_TypeB_AlgorithmFactory.hpp"
#include "ROL_Bounds.hpp"

namespace ROL {
namespace TypeG {

template<typename Real>
StabilizedLCLAlgorithm<Real>::StabilizedLCLAlgorithm( ParameterList &list, const Ptr<Secant<Real>> &secant )
  : TypeG::Algorithm<Real>::Algorithm(), secant_(secant), list_(list), subproblemIter_(0) {
  // Set status test
  status_->reset();
  status_->add(makePtr<ConstraintStatusTest<Real>>(list));

  Real one(1), p1(0.1), p9(0.9), ten(1.e1), oe8(1.e8), oem8(1.e-8);
  ParameterList& sublist = list.sublist("Step").sublist("Stabilized LCL");
  useDefaultInitPen_     = sublist.get("Use Default Initial Penalty Parameter", true);
  state_->searchSize     = sublist.get("Initial Penalty Parameter",             ten);
  sigma_                 = sublist.get("Initial Elastic Penalty Parameter",     ten*ten);
  // Multiplier update parameters
  scaleLagrangian_      = sublist.get("Use Scaled Stabilized LCL",          false);
  penaltyUpdate_        = sublist.get("Penalty Parameter Growth Factor",          ten);
  maxPenaltyParam_      = sublist.get("Maximum Penalty Parameter",                oe8);
  sigmaMax_             = sublist.get("Maximum Elastic Penalty Parameter",        oe8);
  sigmaUpdate_          = sublist.get("Elastic Penalty Parameter Growth Rate",    ten);
  // Optimality tolerance update
  optIncreaseExponent_ = sublist.get("Optimality Tolerance Increase Exponent", one);
  optDecreaseExponent_ = sublist.get("Optimality Tolerance Decrease Exponent", one);
  optToleranceInitial_ = sublist.get("Initial Optimality Tolerance",           one);
  // Feasibility tolerance update    
  feasIncreaseExponent_ = sublist.get("Feasibility Tolerance Increase Exponent", p9);
  feasDecreaseExponent_ = sublist.get("Feasibility Tolerance Decrease Exponent", p1);
  feasToleranceInitial_ = sublist.get("Initial Feasibility Tolerance",           one);
  // Subproblem information
  maxit_         = sublist.get("Subproblem Iteration Limit",              1000);
  subStep_       = sublist.get("Subproblem Step Type",                    "Trust Region");
  HessianApprox_ = sublist.get("Level of Hessian Approximation",          0); 
  list_.sublist("Step").set("Type",subStep_);
  list_.sublist("Status Test").set("Iteration Limit",maxit_);
  // Verbosity setting
  verbosity_          = list.sublist("General").get("Output Level", 0);
  printHeader_        = verbosity_ > 2;
  bool print          = verbosity_ > 2;
  list_.sublist("General").set("Output Level",(print ? verbosity_ : 0));
  // Outer iteration tolerances
  outerFeasTolerance_ = list.sublist("Status Test").get("Constraint Tolerance", oem8);
  outerOptTolerance_  = list.sublist("Status Test").get("Gradient Tolerance",   oem8);
  outerStepTolerance_ = list.sublist("Status Test").get("Step Tolerance",       oem8);
  // Scaling
  useDefaultScaling_  = sublist.get("Use Default Problem Scaling", true);
  fscale_             = sublist.get("Objective Scaling",           one);
  cscale_             = sublist.get("Constraint Scaling",          one);
}

template<typename Real>
void StabilizedLCLAlgorithm<Real>::initialize( Vector<Real>           &x,
                                               const Vector<Real>     &g,
                                               const Vector<Real>     &l,
                                               const Vector<Real>     &c,
                                               ElasticObjective<Real> &alobj,
                                               BoundConstraint<Real>  &bnd,
                                               Constraint<Real>       &con,
                                               std::ostream           &outStream ) {
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
  alobj.getAugmentedLagrangian()->update(x,UpdateType::Initial,state_->iter);
  state_->value = alobj.getObjectiveValue(x,tol);
  alobj.getAugmentedLagrangian()->gradient(*state_->gradientVec,x,tol);

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
  optTolerance_  = std::max<Real>(TOL*outerOptTolerance_,
                            optToleranceInitial_); ///(one+std::pow(state_->searchSize,optDecreaseExponent_)));
  //optTolerance_  = std::min<Real>(optTolerance_,TOL*state_->gnorm);
  feasTolerance_ = std::max<Real>(TOL*outerFeasTolerance_,
                            feasToleranceInitial_); ///(one+std::pow(state_->searchSize,feasDecreaseExponent_)));

  // Set data
  alobj.reset(l,state_->searchSize,sigma_);

  if (verbosity_ > 1) {
    outStream << std::endl;
    outStream << "Stabilized LCL Initialize"                  << std::endl;
    outStream << "Objective Scaling:  " << fscale_            << std::endl;
    outStream << "Constraint Scaling: " << cscale_            << std::endl;
    outStream << "Penalty Parameter:  " << state_->searchSize << std::endl;
    outStream << std::endl;
  }
}

template<typename Real>
void StabilizedLCLAlgorithm<Real>::run( Problem<Real> &problem,
                                        std::ostream  &outStream ) {
  if (problem.getProblemType() == TYPE_EB) {
    problem.edit();
    problem.finalize(true,verbosity_>3,outStream); // Lump linear and nonlinear constraints
    run(*problem.getPrimalOptimizationVector(),
        *problem.getDualOptimizationVector(),
        *problem.getObjective(),
        *problem.getBoundConstraint(),
        *problem.getConstraint(),
        *problem.getMultiplierVector(),
        *problem.getResidualVector(),
        outStream);
    problem.finalizeIteration();
  }
  else {
    throw Exception::NotImplemented(">>> ROL::TypeG::Algorithm::run : Optimization problem is not Type G!");
  }
}

template<typename Real>
void StabilizedLCLAlgorithm<Real>::run( Vector<Real>          &x,
                                        const Vector<Real>    &g,
                                        Objective<Real>       &obj,
                                        BoundConstraint<Real> &bnd,
                                        Constraint<Real>      &econ,
                                        Vector<Real>          &emul,
                                        const Vector<Real>    &eres,
                                        std::ostream          &outStream ) {
  const Real one(1), oem2(1e-2);
  Real tol(std::sqrt(ROL_EPSILON<Real>())), cnorm(0), lnorm;;
  // Initialize augmented Lagrangian data
  ElasticObjective<Real> alobj(makePtrFromRef(obj),makePtrFromRef(econ),
                               state_->searchSize,sigma_,g,eres,emul,
                               scaleLagrangian_,HessianApprox_);
  initialize(x,g,emul,eres,alobj,bnd,econ,outStream);
  // Define Elastic Subproblem
  Ptr<Vector<Real>>  u = eres.clone(),  v = eres.clone(), c = eres.clone();
  Ptr<Vector<Real>> gu = emul.clone(), gv = emul.clone(), l = emul.clone();
  Ptr<Vector<Real>> s = x.clone(), gs = g.clone(), cvec = eres.clone();
  Ptr<ElasticLinearConstraint<Real>> lcon
    = makePtr<ElasticLinearConstraint<Real>>(makePtrFromRef(x),
                                             makePtrFromRef(econ),
                                             makePtrFromRef(eres));
  std::vector<Ptr<Vector<Real>>> vecList = {s,u,v};
  Ptr<PartitionedVector<Real>> xp  = makePtr<PartitionedVector<Real>>(vecList);
  Ptr<PartitionedVector<Real>> gxp = makePtr<PartitionedVector<Real>>({gs,gu,gv});
  Ptr<Vector<Real>>            lb  = u->clone(); lb->zero();
  std::vector<Ptr<BoundConstraint<Real>>> bndList(3);
  bndList[0] = makePtrFromRef(bnd);
  bndList[1] = makePtr<Bounds<Real>>(*lb,true);
  bndList[2] = makePtr<Bounds<Real>>(*lb,true);
  Ptr<BoundConstraint<Real>>  xbnd
    = makePtr<BoundConstraint_Partitioned<Real>>(bndList,vecList);
  ParameterList ppa_list;
  if (c->dimension() == 1)
    ppa_list.sublist("General").sublist("Polyhedral Projection").set("Type","Dai-Fletcher");
  else
    ppa_list.sublist("General").sublist("Polyhedral Projection").set("Type","Semismooth Newton");
  Problem<Real> elc(makePtrFromRef(alobj),xp,gxp);
  elc.addBoundConstraint(xbnd);
  elc.addLinearConstraint("ElasticLinearConstraint",lcon,l,c);
  elc.setProjectionAlgorithm(ppa_list);
  elc.finalize(false,verbosity_>2,outStream);
 
  // Initialize subproblem algorithm
  Ptr<TypeB::Algorithm<Real>> algo;
  
  // Output
  if (verbosity_ > 0) writeOutput(outStream,true);

  while (status_->check(*state_)) {
    lcon->setAnchor(state_->iterateVec);
    if (verbosity_ > 3) elc.check(true,outStream);

    // Solve linearly constrained augmented Lagrangian subproblem
    list_.sublist("Status Test").set("Gradient Tolerance",optTolerance_);
    list_.sublist("Status Test").set("Step Tolerance",1.e-6*optTolerance_);
    algo = TypeB::AlgorithmFactory<Real>(list_,secant_);
    algo->run(elc,outStream);
    x.set(*xp->get(0));

    // Update evaluation counters
    subproblemIter_ = algo->getState()->iter;
    state_->nfval += alobj.getNumberFunctionEvaluations();
    state_->ngrad += alobj.getNumberGradientEvaluations();
    state_->ncval += alobj.getNumberConstraintEvaluations();

    // Compute step
    state_->stepVec->set(x);
    state_->stepVec->axpy(-one,*state_->iterateVec);
    state_->snorm = state_->stepVec->norm();

    // Update iteration information
    state_->iter++;
    cvec->set(*alobj.getConstraintVec(x,tol));
    cnorm = cvec->norm();
    if ( cscale_*cnorm < feasTolerance_ ) {
      // Update iteration information
      state_->iterateVec->set(x);
      state_->value = alobj.getObjectiveValue(x,tol);
      state_->constraintVec->set(*cvec);
      state_->cnorm = cnorm;

      // Update multipliers
      emul.axpy(static_cast<Real>(-1),*elc.getPolyhedralProjection()->getMultiplier());
      emul.axpy(state_->searchSize*cscale_,state_->constraintVec->dual());

      alobj.getAugmentedLagrangian()->gradient(*state_->gradientVec,x,tol);
      if (scaleLagrangian_) state_->gradientVec->scale(state_->searchSize);
      econ.applyAdjointJacobian(*gs,*elc.getPolyhedralProjection()->getMultiplier(),x,tol);
      state_->gradientVec->axpy(-cscale_,*gs);
      x.axpy(-one/std::min(fscale_,cscale_),state_->gradientVec->dual());
      proj_->project(x,outStream);
      x.axpy(-one,*state_->iterateVec);
      state_->gnorm = x.norm();
      x.set(*state_->iterateVec);

      // Update subproblem information
      lnorm  = elc.getPolyhedralProjection()->getMultiplier()->norm();
      sigma_ = std::min(one+lnorm,sigmaMax_)/(one+state_->searchSize);
      if ( algo->getState()->statusFlag == EXITSTATUS_CONVERGED ) {
        optTolerance_  = std::max(oem2*outerOptTolerance_,
                         optTolerance_/(one + std::pow(state_->searchSize,optIncreaseExponent_)));
      }
      feasTolerance_ = std::max(oem2*outerFeasTolerance_,
                       feasTolerance_/(one + std::pow(state_->searchSize,feasIncreaseExponent_)));

      // Update Algorithm State
      state_->snorm += lnorm + state_->searchSize*cscale_*state_->cnorm;
      state_->lagmultVec->set(emul);
    }
    else {
      // Update subproblem information
      state_->searchSize = std::min(penaltyUpdate_*state_->searchSize,maxPenaltyParam_);
      sigma_            /= sigmaUpdate_;
      optTolerance_      = std::max(oem2*outerOptTolerance_,
                           optToleranceInitial_/(one + std::pow(state_->searchSize,optDecreaseExponent_)));
      feasTolerance_     = std::max(oem2*outerFeasTolerance_,
                           feasToleranceInitial_/(one + std::pow(state_->searchSize,feasDecreaseExponent_)));
    }
    alobj.reset(emul,state_->searchSize,sigma_);

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,printHeader_);
  }
  if (verbosity_ > 0) TypeG::Algorithm<Real>::writeExitStatus(outStream);
}

template<typename Real>
void StabilizedLCLAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  if(verbosity_>1) {
    os << std::string(114,'-') << std::endl;
    os << "Stabilized LCL status output definitions" << std::endl << std::endl;
    os << "  iter    - Number of iterates (steps taken)"            << std::endl;
    os << "  fval    - Objective function value"                    << std::endl;
    os << "  cnorm   - Norm of the constraint violation"            << std::endl;
    os << "  gLnorm  - Norm of the gradient of the Lagrangian"      << std::endl;
    os << "  snorm   - Norm of the step"                            << std::endl;
    os << "  penalty - Penalty parameter"                           << std::endl;
    os << "  sigma   - Elastic Penalty parameter"                   << std::endl;
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
  os << std::setw(10) << std::left << "sigma";
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
void StabilizedLCLAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::endl << "Stabilized LCL Solver (Type G, General Constraints)";
  os << std::endl;
  os << "Subproblem Solver: " << subStep_ << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void StabilizedLCLAlgorithm<Real>::writeOutput( std::ostream& os, const bool print_header ) const {
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
    os << std::setw(10) << std::left << sigma_;
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
    os << std::setw(10) << std::left << sigma_;
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
