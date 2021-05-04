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

#ifndef ROL_TYPEB_PRIMALDUALACTIVESETALGORITHM_DEF_HPP
#define ROL_TYPEB_PRIMALDUALACTIVESETALGORITHM_DEF_HPP

namespace ROL {
namespace TypeB {

template<typename Real>
PrimalDualActiveSetAlgorithm<Real>::PrimalDualActiveSetAlgorithm(ParameterList           &list,
                                                                 const Ptr<Secant<Real>> &secant)
  : secant_(secant), esec_(SECANT_USERDEFINED),
    neps_(-ROL_EPSILON<Real>()), itol_(std::sqrt(ROL_EPSILON<Real>())), hasPoly_(true) {
  // Set status test
  status_->reset();
  status_->add(makePtr<StatusTest<Real>>(list));

  if ( secant_ == nullPtr ) {
    secantName_ = list.sublist("General").sublist("Secant").get("Type","Limited-Memory BFGS");
    esec_       = StringToESecant(secantName_);
    secant_     = SecantFactory<Real>(list);
  }
  else {
    secantName_ = list.sublist("General").sublist("Secant").get("User Defined Secant Name",
                                                                "Unspecified User Defined Secant Method");
  }
  useSecantHessVec_ = list.sublist("General").sublist("Secant").get("Use as Hessian",        false);
  useSecantPrecond_ = list.sublist("General").sublist("Secant").get("Use as Preconditioner", false);

  // Algorithmic parameters
  maxit_ = list.sublist("Step").sublist("Primal Dual Active Set").get("Iteration Limit",             10);
  stol_  = list.sublist("Step").sublist("Primal Dual Active Set").get("Relative Step Tolerance",     1e-8);
  gtol_  = list.sublist("Step").sublist("Primal Dual Active Set").get("Relative Gradient Tolerance", 1e-6);
  scale_ = list.sublist("Step").sublist("Primal Dual Active Set").get("Dual Scaling",                1.0);

  // Krylov parameters
  atolKrylov_  = list.sublist("General").sublist("Krylov").get("Absolute Tolerance", 1e-4);
  rtolKrylov_  = list.sublist("General").sublist("Krylov").get("Relative Tolerance", 1e-2);
  maxitKrylov_ = list.sublist("General").sublist("Krylov").get("Iteration Limit",    100);

  verbosity_    = list.sublist("General").get("Output Level",                     0);
  writeHeader_  = verbosity_ > 2;
}

template<typename Real>
void PrimalDualActiveSetAlgorithm<Real>::initialize(Vector<Real>          &x,
                                                    const Vector<Real>    &g,
                                                    Objective<Real>       &obj,
                                                    BoundConstraint<Real> &bnd,
                                                    std::ostream &outStream) {
  // Initialize projection operator
  if (proj_ == nullPtr) {
    proj_ = makePtr<PolyhedralProjection<Real>>(makePtrFromRef(bnd));
    hasPoly_ = false;
  }
  // Create Krylov solver
  if (hasPoly_) {
    ParameterList list;
    list.sublist("General").sublist("Krylov").set("Absolute Tolerance", atolKrylov_);
    list.sublist("General").sublist("Krylov").set("Relative Tolerance", rtolKrylov_);
    list.sublist("General").sublist("Krylov").set("Iteration Limit",   maxitKrylov_);
    krylovName_ = "GMRES";
    krylov_ = makePtr<GMRES<Real>>(list);
  }
  else {
    krylovName_ = "CR";
    krylov_ = makePtr<ConjugateResiduals<Real>>(atolKrylov_,rtolKrylov_,maxitKrylov_);
  }
  ekv_ = StringToEKrylov(krylovName_);
  // Initialize data
  const Real one(1);
  TypeB::Algorithm<Real>::initialize(x,g);
  // Update approximate gradient and approximate objective function.
  Real ftol = std::sqrt(ROL_EPSILON<Real>());
  proj_->project(x,outStream);
  state_->iterateVec->set(x);
  obj.update(x,UpdateType::Initial,state_->iter);
  state_->value = obj.value(x,ftol); state_->nfval++;
  obj.gradient(*state_->gradientVec,x,ftol); state_->ngrad++;
  state_->stepVec->set(x);
  state_->stepVec->axpy(-one,state_->gradientVec->dual());
  proj_->project(*state_->stepVec,outStream);
  state_->stepVec->axpy(-one,x);
  state_->gnorm = state_->stepVec->norm();
  state_->snorm = ROL_INF<Real>();
}

template<typename Real>
void PrimalDualActiveSetAlgorithm<Real>::run( Vector<Real>          &x,
                                              const Vector<Real>    &g, 
                                              Objective<Real>       &obj,
                                              BoundConstraint<Real> &bnd,
                                              std::ostream          &outStream ) {
  const Real zero(0), one(1);
  // Initialize PDAS data
  initialize(x,g,obj,bnd,outStream);
  Ptr<Vector<Real>> xlam   = x.clone(), xtmp = x.clone(), As = x.clone();
  Ptr<Vector<Real>> lambda = x.clone(), gtmp = g.clone();
  Ptr<Vector<Real>> pwa    = x.clone(), dwa  = g.clone();
  Ptr<Vector<Real>> mu, b, r;
  if (hasPoly_) {
    mu  = proj_->getMultiplier()->clone(); mu->zero();
    b   = proj_->getResidual()->clone();
    r   = proj_->getResidual()->clone();
  }
  lambda->set(state_->gradientVec->dual());
  lambda->scale(-one);
  Real xnorm(0), snorm(0), rnorm(0), tol(std::sqrt(ROL_EPSILON<Real>()));

  Ptr<LinearOperator<Real>> hessian, precond;

  // Output
  if (verbosity_ > 0) writeOutput(outStream,true);

  while (status_->check(*state_)) {
    totalKrylov_ = 0;
    state_->stepVec->zero();
    xnorm = x.norm();
    if (hasPoly_) {
      proj_->getLinearConstraint()->value(*r,x,tol);
    }
    for ( iter_ = 0; iter_ < maxit_; iter_++ ) {
      /********************************************************************/
      // MODIFY ITERATE VECTOR TO CHECK ACTIVE SET
      /********************************************************************/
      xlam->set(x);                              // xlam = x0
      xlam->axpy(scale_,*lambda);                // xlam = x0 + c*lambda
      /********************************************************************/
      // PROJECT x ONTO PRIMAL DUAL FEASIBLE SET
      /********************************************************************/
      As->zero();                                // As   = 0
   
      xtmp->set(*bnd.getUpperBound());           // xtmp = u
      xtmp->axpy(-one,*state_->iterateVec);      // xtmp = u - x
      bnd.pruneUpperInactive(*xtmp,*xlam,neps_); // xtmp = A(u - x)
      As->plus(*xtmp);                           // As  += A(u - x)

      xtmp->set(*bnd.getLowerBound());           // xtmp = l
      xtmp->axpy(-one,*state_->iterateVec);      // xtmp = l - x
      bnd.pruneLowerInactive(*xtmp,*xlam,neps_); // xtmp = A(l - x)
      As->plus(*xtmp);                           // As  += A(l - x)
      /********************************************************************/
      // APPLY HESSIAN TO ACTIVE COMPONENTS OF s AND REMOVE INACTIVE
      /********************************************************************/
      if ( useSecantHessVec_ && secant_ != nullPtr ) { // gtmp = H*As
        secant_->applyB(*gtmp,*As);
      }
      else {
        obj.hessVec(*gtmp,*As,*state_->iterateVec,itol_);
      }
      if (hasPoly_) {
        proj_->getLinearConstraint()->applyJacobian(*b,*As,*state_->iterateVec,tol);
	b->plus(*r);
        b->scale(-one);
      }
      gtmp->plus(*state_->gradientVec);       // gtmp = g + H*As + ...
      gtmp->scale(-one);                      // gtmp = -(g + H*As + ...)
      bnd.pruneActive(*gtmp,*xlam,neps_);     // gtmp = -I(g + H*As + ...)
      /********************************************************************/
      // SOLVE REDUCED NEWTON SYSTEM 
      /********************************************************************/
      if (hasPoly_) {
        rnorm = std::sqrt(gtmp->dot(*gtmp)+b->dot(*b));
      }
      else {
        rnorm = gtmp->norm();
      }
      if ( rnorm > zero ) {
        if (hasPoly_) {
          // Initialize Hessian and preconditioner
          hessian = makePtr<HessianPDAS_Poly>(makePtrFromRef(obj),makePtrFromRef(bnd),
                proj_->getLinearConstraint(),state_->iterateVec,xlam,neps_,secant_,
                useSecantHessVec_,pwa,dwa);
          precond = makePtr<PrecondPDAS_Poly>(makePtrFromRef(obj),makePtrFromRef(bnd),
                state_->iterateVec,xlam,neps_,secant_,useSecantPrecond_,dwa);
          PartitionedVector<Real> rhs(std::vector<Ptr<Vector<Real>>>({gtmp,b}));
          PartitionedVector<Real> sol(std::vector<Ptr<Vector<Real>>>({state_->stepVec,mu}));
          krylov_->run(sol,*hessian,rhs,*precond,iterKrylov_,flagKrylov_);
        }
        else {
          // Initialize Hessian and preconditioner
          hessian = makePtr<HessianPDAS>(makePtrFromRef(obj),makePtrFromRef(bnd),
                state_->iterateVec,xlam,neps_,secant_,useSecantHessVec_,pwa);
          precond = makePtr<PrecondPDAS>(makePtrFromRef(obj),makePtrFromRef(bnd),
                state_->iterateVec,xlam,neps_,secant_,useSecantPrecond_,dwa);
          krylov_->run(*state_->stepVec,*hessian,*gtmp,*precond,iterKrylov_,flagKrylov_);
        }
        totalKrylov_ += iterKrylov_;
        bnd.pruneActive(*state_->stepVec,*xlam,neps_);     // s <- Is
      }
      state_->stepVec->plus(*As);                          // s = Is + As
      /********************************************************************/
      // UPDATE STEP AND MULTIPLIER 
      /********************************************************************/
      x.set(*state_->iterateVec);
      x.plus(*state_->stepVec);
      snorm = state_->stepVec->norm();
      // Compute gradient of Lagrangian for QP
      if ( useSecantHessVec_ && secant_ != nullPtr ) {
        secant_->applyB(*gtmp,*state_->stepVec);
      }
      else {
        obj.hessVec(*gtmp,*state_->stepVec,*state_->iterateVec,itol_);
      }
      if (hasPoly_) {
        proj_->getLinearConstraint()->applyAdjointJacobian(*dwa,*mu,*state_->iterateVec,tol);
        gtmp->plus(*dwa);
      }
      gtmp->plus(*state_->gradientVec);
      gtmp->scale(-one);
      lambda->set(gtmp->dual());
      // Compute criticality measure  
      xtmp->set(x);
      xtmp->plus(*lambda);
      bnd.project(*xtmp);
      xtmp->axpy(-one,x);
      // Update multiplier
      bnd.pruneInactive(*lambda,*xlam,neps_);
      // Check stopping conditions
      if ( xtmp->norm() < gtol_*state_->gnorm ) {
        flag_ = 0;
        break;
      }
      if ( snorm < stol_*xnorm ) {
        flag_ = 2;
        break;
      }
    }
    if ( iter_ == maxit_ ) {
      flag_ = 1;
    }
    else {
      iter_++;
    }

    // Update algorithm state
    state_->iter++;
    state_->iterateVec->set(x);
    feasible_ = bnd.isFeasible(x);
    state_->snorm = snorm;
    obj.update(x,UpdateType::Accept,state_->iter);
    state_->value = obj.value(x,tol); state_->nfval++;
    
    if ( secant_ != nullPtr ) {
      gtmp->set(*state_->gradientVec);
    }
    obj.gradient(*state_->gradientVec,x,tol); state_->ngrad++;
    xtmp->set(x); xtmp->axpy(-one,state_->gradientVec->dual());
    proj_->project(*xtmp,outStream);
    xtmp->axpy(-one,x);
    state_->gnorm = xtmp->norm();
    if ( secant_ != nullPtr ) {
      secant_->updateStorage(x,*state_->gradientVec,*gtmp,*state_->stepVec,state_->snorm,state_->iter+1);
    }

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,writeHeader_);
  }
  if (verbosity_ > 0) TypeB::Algorithm<Real>::writeExitStatus(outStream);
}

template<typename Real>
void PrimalDualActiveSetAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::stringstream hist;
  if (verbosity_ > 1) {
    hist << std::string(114,'-') << std::endl;
    if (!useSecantHessVec_) {
      hist << "Primal Dual Active Set Newton's Method";
    }
    else {
      hist << "Primal Dual Active Set Quasi-Newton Method with " << secantName_ << " Hessian approximation";
    }
    hist << " status output definitions" << std::endl << std::endl;
    hist << "  iter       - Number of iterates (steps taken)" << std::endl;
    hist << "  value      - Objective function value" << std::endl;
    hist << "  gnorm      - Norm of the gradient" << std::endl;
    hist << "  snorm      - Norm of the step (update to optimization vector)" << std::endl;
    hist << "  #fval      - Cumulative number of times the objective function was evaluated" << std::endl;
    hist << "  #grad      - Cumulative number of times the gradient was computed" << std::endl;
    if (maxit_ > 1) {
      hist << "  iterPDAS   - Number of Primal Dual Active Set iterations" << std::endl << std::endl;
      hist << "  flagPDAS   - Primal Dual Active Set flag" << std::endl;
      hist << "  iterK      - Number of Krylov iterations" << std::endl << std::endl;
    }
    else {
      hist << "  iterK      - Number of Krylov iterations" << std::endl << std::endl;
      hist << "  flagK      - Krylov flag" << std::endl;
      for( int flag = CG_FLAG_SUCCESS; flag != CG_FLAG_UNDEFINED; ++flag ) {
        hist << "    " << NumberToString(flag) << " - "
             << ECGFlagToString(static_cast<ECGFlag>(flag)) << std::endl;
      }            
    }
    hist << "  feasible - Is iterate feasible?" << std::endl;
    hist << std::string(114,'-') << std::endl;
  }

  hist << "  ";
  hist << std::setw(6)  << std::left << "iter";
  hist << std::setw(15) << std::left << "value";
  hist << std::setw(15) << std::left << "gnorm";
  hist << std::setw(15) << std::left << "snorm";
  hist << std::setw(10) << std::left << "#fval";
  hist << std::setw(10) << std::left << "#grad";
  if (maxit_ > 1) {
    hist << std::setw(10) << std::left << "iterPDAS";
    hist << std::setw(10) << std::left << "flagPDAS";
    hist << std::setw(10) << std::left << "iterK";
  }
  else {
    hist << std::setw(10) << std::left << "iterK";
    hist << std::setw(10) << std::left << "flagK";
  }
  hist << std::setw(10) << std::left << "feasible";
  hist << std::endl;
  os << hist.str();
}

template<typename Real>
void PrimalDualActiveSetAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::stringstream hist;
  if (!useSecantHessVec_) {
    hist << std::endl << "Primal Dual Active Set Newton's Method (Type B, Bound Constraints)" << std::endl;
  }
  else {
    hist << std::endl << "Primal Dual Active Set Quasi-Newton Method with "
         << secantName_ << " Hessian approximation" << std::endl;
  }
  os << hist.str();
}

template<typename Real>
void PrimalDualActiveSetAlgorithm<Real>::writeOutput( std::ostream& os, bool write_header ) const {
  std::stringstream hist;
  hist << std::scientific << std::setprecision(6);
  if ( state_->iter == 0 ) writeName(os);
  if ( write_header )      writeHeader(os);
  if ( state_->iter == 0 ) {
    hist << "  ";
    hist << std::setw(6)  << std::left << state_->iter;
    hist << std::setw(15) << std::left << state_->value;
    hist << std::setw(15) << std::left << state_->gnorm;
    hist << std::setw(15) << std::left << "---";
    hist << std::setw(10) << std::left << state_->nfval;
    hist << std::setw(10) << std::left << state_->ngrad;
    hist << std::setw(10) << std::left << "---";
    hist << std::setw(10) << std::left << "---";
    if (maxit_ > 1) {
      hist << std::setw(10) << std::left << "---";
    }
    if ( feasible_ ) {
      hist << std::setw(10) << std::left << "YES";
    }
    else {
      hist << std::setw(10) << std::left << "NO";
    }
    hist << std::endl;
  }
  else {
    hist << "  ";
    hist << std::setw(6)  << std::left << state_->iter;
    hist << std::setw(15) << std::left << state_->value;
    hist << std::setw(15) << std::left << state_->gnorm;
    hist << std::setw(15) << std::left << state_->snorm;
    hist << std::setw(10) << std::left << state_->nfval;
    hist << std::setw(10) << std::left << state_->ngrad;
    if (maxit_ > 1) {
      hist << std::setw(10) << std::left << iter_;
      hist << std::setw(10) << std::left << flag_;
      hist << std::setw(10) << std::left << totalKrylov_;
    }
    else {
      hist << std::setw(10) << std::left << iterKrylov_;
      hist << std::setw(10) << std::left << flagKrylov_;
    }
    if ( feasible_ ) {
      hist << std::setw(10) << std::left << "YES";
    }
    else {
      hist << std::setw(10) << std::left << "NO";
    }
    hist << std::endl;
  }
  os << hist.str();
}

} // namespace TypeB
} // namespace ROL

#endif
