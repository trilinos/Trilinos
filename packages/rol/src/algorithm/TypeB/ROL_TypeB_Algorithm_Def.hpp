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

#ifndef ROL_TYPEB_ALGORITHM_DEF_H
#define ROL_TYPEB_ALGORITHM_DEF_H

#include "ROL_SlacklessObjective.hpp"
//#include "ROL_ConstraintManager.hpp"

namespace ROL {
namespace TypeB {

template<typename Real>
Algorithm<Real>::Algorithm()
  : status_(makePtr<CombinedStatusTest<Real>>()),
    state_(makePtr<AlgorithmState<Real>>()),
    proj_(nullPtr) {
  status_->reset();
  status_->add(makePtr<StatusTest<Real>>());
}

template<typename Real>
void Algorithm<Real>::initialize(const Vector<Real> &x, const Vector<Real> &g) {
  if (state_->iterateVec == nullPtr) {
    state_->iterateVec = x.clone();
  }
  state_->iterateVec->set(x);
  if (state_->stepVec == nullPtr) {
    state_->stepVec = x.clone();
  }
  state_->stepVec->zero();
  if (state_->gradientVec == nullPtr) {
    state_->gradientVec = g.clone();
  }
  state_->gradientVec->set(g);
  if (state_->minIterVec == nullPtr) {
    state_->minIterVec = x.clone();
  }
  state_->minIterVec->set(x);
  state_->minIter = state_->iter;
  state_->minValue = state_->value;
}

template<typename Real>
Real Algorithm<Real>::optimalityCriterion(const Vector<Real> &x,
                                           const Vector<Real> &g,
                                           Vector<Real> &primal,
                                           std::ostream &outStream) const {
  const Real one(1);
  primal.set(x);
  primal.axpy(-one,g.dual());
  proj_->project(primal,outStream); state_->nproj++;
  primal.axpy(-one,x);
  return primal.norm();
}

template<typename Real>
void Algorithm<Real>::setStatusTest(const Ptr<StatusTest<Real>> &status,
                                    bool combineStatus) {
  if (!combineStatus) { // Do not combine status tests
    status_->reset();
  }
  status_->add(status); // Add user-defined StatusTest
}

template<typename Real>
void Algorithm<Real>::run( Problem<Real> &problem,
                           std::ostream  &outStream ) {
  if (problem.getProblemType() == TYPE_B) {
    proj_ = problem.getPolyhedralProjection();
    run(*problem.getPrimalOptimizationVector(),
        *problem.getDualOptimizationVector(),
        *problem.getObjective(),
        *problem.getBoundConstraint(),
        outStream);
    problem.finalizeIteration();
  }
  else {
    throw Exception::NotImplemented(">>> ROL::Algorithm::TypeB::run : Optimization problem is not Type B!");
  }
}

template<typename Real>
void Algorithm<Real>::run( Vector<Real>          &x,
                           Objective<Real>       &obj,
                           BoundConstraint<Real> &bnd,
                           std::ostream          &outStream ) {
  run(x,x.dual(),obj,bnd,outStream);
}

template<typename Real>
void Algorithm<Real>::run( Vector<Real>          &x,
                           Objective<Real>       &obj,
                           BoundConstraint<Real> &bnd,
                           Constraint<Real>      &linear_econ,
                           Vector<Real>          &linear_emul,
                           std::ostream          &outStream ) {
  run(x,x.dual(),obj,bnd,linear_econ,linear_emul,linear_emul.dual(),outStream);
}

template<typename Real>
void Algorithm<Real>::run( Vector<Real>          &x,
                           const Vector<Real>    &g,
                           Objective<Real>       &obj,
                           BoundConstraint<Real> &bnd,
                           Constraint<Real>      &linear_econ,
                           Vector<Real>          &linear_emul,
                           const Vector<Real>    &linear_eres,
                           std::ostream          &outStream ) {
  ParameterList parlist;
  proj_ = PolyhedralProjectionFactory<Real>(x,g,makePtrFromRef(bnd),makePtrFromRef(linear_econ),linear_emul,linear_eres,parlist);
  run(x,g,obj,bnd,outStream);
}

template<typename Real>
void Algorithm<Real>::run( Vector<Real>          &x,
                           Objective<Real>       &obj,
                           BoundConstraint<Real> &bnd,
                           Constraint<Real>      &linear_icon,
                           Vector<Real>          &linear_imul,
                           BoundConstraint<Real> &linear_ibnd,
                           std::ostream          &outStream ) {
  run(x,x.dual(),obj,bnd,linear_icon,linear_imul,linear_ibnd,linear_imul.dual(),outStream);
}

template<typename Real>
void Algorithm<Real>::run( Vector<Real>          &x,
                           const Vector<Real>    &g,
                           Objective<Real>       &obj,
                           BoundConstraint<Real> &bnd,
                           Constraint<Real>      &linear_icon,
                           Vector<Real>          &linear_imul,
                           BoundConstraint<Real> &linear_ibnd,
                           const Vector<Real>    &linear_ires,
                           std::ostream          &outStream ) {
  Ptr<Vector<Real>> gp = g.clone(), irp = linear_ires.clone();
  Problem<Real> problem(makePtrFromRef(obj),
                        makePtrFromRef(x),gp);
  problem.addBoundConstraint(makePtrFromRef(bnd));
  problem.addLinearConstraint("LinearInequalityConstraint",
                              makePtrFromRef(linear_icon),
                              makePtrFromRef(linear_imul),
                              makePtrFromRef(linear_ibnd),irp);
  problem.finalize(false,false,outStream);
  run(problem,outStream);
//  ConstraintManager<Real> cm(makePtrFromRef(linear_icon),makePtrFromRef(linear_imul),
//                             makePtrFromRef(linear_ibnd),makePtrFromRef(x),
//                             makePtrFromRef(bnd));
//  Ptr<Constraint<Real>>      linear_econ = cm.getConstraint();
//  Ptr<Vector<Real>>          linear_emul = cm.getMultiplier();
//  Ptr<Vector<Real>>          xvec = cm.getOptVector();
//  Ptr<BoundConstraint<Real>> xbnd = cm.getBoundConstraint();
//  Ptr<Objective<Real>>       sobj = makePtr<SlacklessObjective<Real>>(makePtrFromRef(obj));
//  //run(*xvec,xvec->dual(),*sobj,*xbnd,*linear_econ,*linear_emul,linear_emul->dual(),outStream);
//  Ptr<Vector<Real>>          xdual = xvec->dual().clone();
//  run(*xvec,*xdual,*sobj,*xbnd,*linear_econ,*linear_emul,linear_emul->dual(),outStream);
}

template<typename Real>
void Algorithm<Real>::run( Vector<Real>          &x,
                           Objective<Real>       &obj,
                           BoundConstraint<Real> &bnd,
                           Constraint<Real>      &linear_econ,
                           Vector<Real>          &linear_emul,
                           Constraint<Real>      &linear_icon,
                           Vector<Real>          &linear_imul,
                           BoundConstraint<Real> &linear_ibnd,
                           std::ostream          &outStream ) {
  run(x,x.dual(),obj,bnd,linear_econ,linear_emul,linear_emul.dual(),
             linear_icon,linear_imul,linear_ibnd,linear_imul.dual(),outStream);
}

template<typename Real>
void Algorithm<Real>::run( Vector<Real>          &x,
                           const Vector<Real>    &g,
                           Objective<Real>       &obj,
                           BoundConstraint<Real> &bnd,
                           Constraint<Real>      &linear_econ,
                           Vector<Real>          &linear_emul,
                           const Vector<Real>    &linear_eres,
                           Constraint<Real>      &linear_icon,
                           Vector<Real>          &linear_imul,
                           BoundConstraint<Real> &linear_ibnd,
                           const Vector<Real>    &linear_ires,
                           std::ostream          &outStream ) {
  Ptr<Vector<Real>> gp = g.clone(), erp = linear_eres.clone(), irp = linear_ires.clone();
  Problem<Real> problem(makePtrFromRef(obj),
                        makePtrFromRef(x),gp);
  problem.addBoundConstraint(makePtrFromRef(bnd));
  problem.addLinearConstraint("LinearEqualityConstraint",
                              makePtrFromRef(linear_econ),
                              makePtrFromRef(linear_emul),erp);
  problem.addLinearConstraint("LinearInequalityConstraint",
                              makePtrFromRef(linear_icon),
                              makePtrFromRef(linear_imul),
                              makePtrFromRef(linear_ibnd),irp);
  problem.finalize(false,false,outStream);
  run(problem,outStream);
  //std::vector<Ptr<Constraint<Real>>> cvec;
  //cvec.push_back(makePtrFromRef(linear_econ));
  //cvec.push_back(makePtrFromRef(linear_icon));
  //std::vector<Ptr<Vector<Real>>> lvec;
  //lvec.push_back(makePtrFromRef(linear_emul));
  //lvec.push_back(makePtrFromRef(linear_imul));
  //std::vector<Ptr<BoundConstraint<Real>>> bvec;
  //bvec.push_back(nullPtr);
  //bvec.push_back(makePtrFromRef(linear_ibnd));
  //ConstraintManager<Real> cm(cvec,lvec,bvec,makePtrFromRef(x),makePtrFromRef(bnd));
  //Ptr<Constraint<Real>>      linear_con = cm.getConstraint();
  //Ptr<Vector<Real>>          linear_mul = cm.getMultiplier();
  //Ptr<Vector<Real>>          xvec = cm.getOptVector();
  //Ptr<BoundConstraint<Real>> xbnd = cm.getBoundConstraint();
  //Ptr<Objective<Real>>       xobj = makePtr<SlacklessObjective<Real>>(makePtrFromRef(obj));
  ////run(*xvec,xvec->dual(),*xobj,*xbnd,*linear_con,*linear_mul,linear_mul->dual(),outStream);
  //Ptr<Vector<Real>>          xdual = xvec->dual().clone();
  //run(*xvec,*xdual,*xobj,*xbnd,*linear_con,*linear_mul,linear_mul->dual(),outStream);
}

template<typename Real>
void Algorithm<Real>::writeHeader( std::ostream& os ) const {
  std::stringstream hist;
  hist << "  ";
  hist << std::setw(6)  << std::left << "iter";
  hist << std::setw(15) << std::left << "value";
  hist << std::setw(15) << std::left << "gnorm";
  hist << std::setw(15) << std::left << "snorm";
  hist << std::setw(10) << std::left << "#fval";
  hist << std::setw(10) << std::left << "#grad";
  hist << std::endl;
  os << hist.str();
}

template<typename Real>
void Algorithm<Real>::writeName( std::ostream &os ) const {
  throw Exception::NotImplemented(">>> ROL::TypeU::Algorithm::writeName() is not implemented!");
}

template<typename Real>
void Algorithm<Real>::writeOutput( std::ostream& os, bool write_header ) const {
  std::stringstream hist;
  hist << std::scientific << std::setprecision(6);
  if ( write_header ) writeHeader(os);
  if ( state_->iter == 0 ) {
    hist << "  ";
    hist << std::setw(6)  << std::left << state_->iter;
    hist << std::setw(15) << std::left << state_->value;
    hist << std::setw(15) << std::left << state_->gnorm;
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
    hist << std::endl;
  }
  os << hist.str();
}

template<typename Real>
void Algorithm<Real>::writeExitStatus( std::ostream& os ) const {
  std::stringstream hist;
  hist << "Optimization Terminated with Status: ";
  hist << EExitStatusToString(state_->statusFlag);
  hist << std::endl;
  os << hist.str();
}

template<typename Real>
//Ptr<const AlgorithmState<Real>>& Algorithm<Real>::getState() const {
Ptr<const AlgorithmState<Real>> Algorithm<Real>::getState() const {
  return state_;
}

template<typename Real>
void Algorithm<Real>::reset() {
  state_->reset();
}

} // namespace TypeB
} // namespace ROL

#endif
