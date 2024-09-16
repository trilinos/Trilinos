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


#ifndef ROL_REDUCED_CONSTRAINT_SIMOPT_DEF_H
#define ROL_REDUCED_CONSTRAINT_SIMOPT_DEF_H

namespace ROL {

template<class Real>
void Reduced_Constraint_SimOpt<Real>::solve_state_equation(const Vector<Real> &z, Real &tol) {
  if (!isUpdated_) {
    // Update equality constraint with new Opt variable.
    if (newUpdate_) conRed_->update_2(z,updateType_,updateIter_);
    else            conRed_->update_2(z,updateFlag_,updateIter_);
  }
  // Check if state has been computed.
  bool isComputed = storage_ ? stateStore_->get(*state_,Constraint<Real>::getParameter()) : false;
  // Solve state equation if not done already.
  if (!isComputed || !storage_) {
    // Solve state equation.
    conRed_->solve(*dualadjoint_,*state_,z,tol);
    nstat_++;
    // Store state.
    if (storage_) stateStore_->set(*state_,Constraint<Real>::getParameter());
  }
  if (!isUpdated_) {
    // Update equality constraint with new Sim variable.
    if (newUpdate_) conRed_->update_1(*state_,updateType_,updateIter_);
    else            conRed_->update_1(*state_,updateFlag_,updateIter_);
    // Update full objective function.
    if (newUpdate_) conVal_->update(*state_,z,updateType_,updateIter_);
    else            conVal_->update(*state_,z,updateFlag_,updateIter_);
    isUpdated_ = true;
  }
}

template<class Real>
void Reduced_Constraint_SimOpt<Real>::solve_adjoint_equation(const Vector<Real> &w, const Vector<Real> &z, Real &tol) { 
  // Check if adjoint has been computed.
  bool isComputed = storage_ ? adjointStore_->get(*adjoint_,Constraint<Real>::getParameter()) : false;
  // Solve adjoint equation if not done already.
  if (!isComputed || !storage_) {
    // Evaluate the full gradient wrt u
    conVal_->applyAdjointJacobian_1(*dualstate_,w,*state_,z,tol);
    // Solve adjoint equation
    conRed_->applyInverseAdjointJacobian_1(*adjoint_,*dualstate_,*state_,z,tol);
    adjoint_->scale(static_cast<Real>(-1));
    nadjo_++;
    // Store adjoint
    if (storage_) adjointStore_->set(*adjoint_,Constraint<Real>::getParameter());
  }
}

template<class Real>
void Reduced_Constraint_SimOpt<Real>::solve_state_sensitivity(const Vector<Real> &v, const Vector<Real> &z, Real &tol) {
  // Solve state sensitivity equation
  conRed_->applyJacobian_2(*dualadjoint_,v,*state_,z,tol);
  dualadjoint_->scale(static_cast<Real>(-1));
  conRed_->applyInverseJacobian_1(*state_sens_,*dualadjoint_,*state_,z,tol);
  nssen_++;
}

template<class Real>
void Reduced_Constraint_SimOpt<Real>::solve_adjoint_sensitivity(const Vector<Real> &w, const Vector<Real> &v, const Vector<Real> &z, Real &tol) {
  // Evaluate full hessVec in the direction (s,v)
  conVal_->applyAdjointHessian_11(*dualstate_,w,*state_sens_,*state_,z,tol);
  conVal_->applyAdjointHessian_21(*dualstate1_,w,v,*state_,z,tol);
  dualstate_->plus(*dualstate1_);
  // Apply adjoint Hessian of constraint
  conRed_->applyAdjointHessian_11(*dualstate1_,*adjoint_,*state_sens_,*state_,z,tol);
  dualstate_->plus(*dualstate1_);
  conRed_->applyAdjointHessian_21(*dualstate1_,*adjoint_,v,*state_,z,tol);
  dualstate_->plus(*dualstate1_);
  // Solve adjoint sensitivity equation
  dualstate_->scale(static_cast<Real>(-1));
  conRed_->applyInverseAdjointJacobian_1(*adjoint_sens_,*dualstate_,*state_,z,tol);
  nasen_++;
}

template<class Real>
Reduced_Constraint_SimOpt<Real>::Reduced_Constraint_SimOpt(
    const ROL::Ptr<Constraint_SimOpt<Real>> &conVal,
    const ROL::Ptr<Constraint_SimOpt<Real>> &conRed,
    const ROL::Ptr<VectorController<Real>> &stateStore,
    const ROL::Ptr<Vector<Real>> &state,
    const ROL::Ptr<Vector<Real>> &control,
    const ROL::Ptr<Vector<Real>> &adjoint,
    const ROL::Ptr<Vector<Real>> &residual,
    bool storage,
    bool useFDhessVec)
  : conVal_(       conVal ),
    conRed_(       conRed ),
    stateStore_(   stateStore ),
    adjointStore_( ROL::makePtr<VectorController<Real>>() ),
    state_(        state->clone() ),
    adjoint_(      adjoint->clone() ),
    residual_(     residual->clone() ),
    state_sens_(   state->clone() ),
    adjoint_sens_( adjoint->clone() ),
    dualstate_(    state->dual().clone() ),
    dualstate1_(   state->dual().clone() ),
    dualadjoint_(  adjoint->dual().clone() ),
    dualcontrol_(  control->dual().clone() ), 
    dualresidual_( residual->dual().clone() ),
    storage_(storage), useFDhessVec_(useFDhessVec),
    nupda_(0), nvalu_(0), njaco_(0), najac_(0), nhess_(0),
    nstat_(0), nadjo_(0), nssen_(0), nasen_(0),
    updateFlag_(true), updateIter_(0), updateType_(UpdateType::Initial),
    newUpdate_(false), isUpdated_(true) {}

template<class Real>
Reduced_Constraint_SimOpt<Real>::Reduced_Constraint_SimOpt(
    const ROL::Ptr<Constraint_SimOpt<Real>> &conVal,
    const ROL::Ptr<Constraint_SimOpt<Real>> &conRed,
    const ROL::Ptr<VectorController<Real>> &stateStore,
    const ROL::Ptr<Vector<Real>> &state,
    const ROL::Ptr<Vector<Real>> &control,
    const ROL::Ptr<Vector<Real>> &adjoint,
    const ROL::Ptr<Vector<Real>> &residual,
    const ROL::Ptr<Vector<Real>> &dualstate,
    const ROL::Ptr<Vector<Real>> &dualcontrol,
    const ROL::Ptr<Vector<Real>> &dualadjoint,
    const ROL::Ptr<Vector<Real>> &dualresidual,
    bool storage,
    bool useFDhessVec)
  : conVal_(       conVal ),
    conRed_(       conRed ),
    stateStore_(   stateStore ),
    adjointStore_( ROL::makePtr<VectorController<Real>>() ),
    state_(        state->clone() ),
    adjoint_(      adjoint->clone() ),
    residual_(     residual->clone() ),
    state_sens_(   state->clone() ),
    adjoint_sens_( adjoint->clone() ),
    dualstate_(    dualstate->clone() ),
    dualstate1_(   dualstate->clone() ),
    dualadjoint_(  dualadjoint->clone() ),
    dualcontrol_(  dualcontrol->clone() ), 
    dualresidual_( dualresidual->clone() ),
    storage_(storage), useFDhessVec_(useFDhessVec),
    nupda_(0), nvalu_(0), njaco_(0), najac_(0), nhess_(0),
    nstat_(0), nadjo_(0), nssen_(0), nasen_(0),
    updateFlag_(true), updateIter_(0), updateType_(UpdateType::Initial),
    newUpdate_(false), isUpdated_(true) {}

template<class Real>
void Reduced_Constraint_SimOpt<Real>::summarize(std::ostream &stream, const Ptr<BatchManager<Real>> &bman) const {
  int nupda(0), nvalu(0), njaco(0), najac(0), nhess(0), nstat(0), nadjo(0), nssen(0), nasen(0);
  if (bman == nullPtr) {
    nupda = nupda_;
    nvalu = nvalu_;
    njaco = njaco_;
    najac = najac_;
    nhess = nhess_;
    nstat = nstat_;
    nadjo = nadjo_;
    nssen = nssen_;
    nasen = nasen_;
  }
  else {
    auto sumAll = [bman](int val) {
      Real global(0), local(val);
      bman->sumAll(&local,&global,1);
      return static_cast<int>(global);
    };
    nupda = sumAll(nupda_);
    nvalu = sumAll(nvalu_);
    njaco = sumAll(njaco_);
    najac = sumAll(najac_);
    nhess = sumAll(nhess_);
    nstat = sumAll(nstat_);
    nadjo = sumAll(nadjo_);
    nssen = sumAll(nssen_);
    nasen = sumAll(nasen_);
  }
  stream << std::endl;
  stream << std::string(80,'=') << std::endl;
  stream << "  ROL::Reduced_Objective_SimOpt::summarize" << std::endl;
  stream << "    Number of calls to update:               " << nupda << std::endl;
  stream << "    Number of calls to value:                " << nvalu << std::endl;
  stream << "    Number of calls to applyJacobian:        " << njaco << std::endl;
  stream << "    Number of calls to applyAdjointJacobian: " << najac << std::endl;
  stream << "    Number of calls to hessvec:              " << nhess << std::endl;
  stream << "    Number of state solves:                  " << nstat << std::endl;
  stream << "    Number of adjoint solves:                " << nadjo << std::endl;
  stream << "    Number of state sensitivity solves:      " << nssen << std::endl;
  stream << "    Number of adjoint sensitivity solves:    " << nasen << std::endl;
  stream << std::string(80,'=') << std::endl;
  stream << std::endl;
}

template<class Real>
void Reduced_Constraint_SimOpt<Real>::reset() {
  nupda_ = 0; nvalu_ = 0; njaco_ = 0; najac_ = 0; nhess_ = 0;
  nstat_ = 0; nadjo_ = 0; nssen_ = 0; nasen_ = 0;
}

template<class Real>
void Reduced_Constraint_SimOpt<Real>::update( const Vector<Real> &z, bool flag, int iter ) {
  nupda_++;
  updateFlag_ = flag;
  updateIter_ = iter;
  stateStore_->constraintUpdate(true);
  adjointStore_->constraintUpdate(flag);
}
template<class Real>
void Reduced_Constraint_SimOpt<Real>::update( const Vector<Real> &z, UpdateType type, int iter ) {
  nupda_++;
  isUpdated_  = false;
  newUpdate_  = true;
  updateType_ = type;
  updateIter_ = iter;
  stateStore_->objectiveUpdate(type);
  adjointStore_->objectiveUpdate(type);
}

template<class Real>
void Reduced_Constraint_SimOpt<Real>::value( Vector<Real> &c, const Vector<Real> &z, Real &tol ) {
  nvalu_++;
  // Solve state equation
  solve_state_equation(z,tol);
  // Get constraint value
  conVal_->value(c,*state_,z,tol);
}

template<class Real>
void Reduced_Constraint_SimOpt<Real>::applyJacobian( Vector<Real> &jv, const Vector<Real> &v,
                    const Vector<Real> &z, Real &tol ) {
  njaco_++;
  // Solve state equation.
  solve_state_equation(z,tol);
  // Solve state sensitivity equation.
  solve_state_sensitivity(v,z,tol);
  // Apply Sim Jacobian to state sensitivity.
  conVal_->applyJacobian_1(*residual_,*state_sens_,*state_,z,tol);
  // Apply Opt Jacobian to vector.
  conVal_->applyJacobian_2(jv,v,*state_,z,tol);
  jv.plus(*residual_);
}

template<class Real>
void Reduced_Constraint_SimOpt<Real>::applyAdjointJacobian( Vector<Real> &ajw, const Vector<Real> &w,
                           const Vector<Real> &z, Real &tol ) {
  najac_++;
  // Solve state equation
  solve_state_equation(z,tol);
  // Solve adjoint equation
  solve_adjoint_equation(w,z,tol);
  // Evaluate the full gradient wrt z
  conVal_->applyAdjointJacobian_2(*dualcontrol_,w,*state_,z,tol);
  // Build gradient
  conRed_->applyAdjointJacobian_2(ajw,*adjoint_,*state_,z,tol);
  ajw.plus(*dualcontrol_);
}

template<class Real>
void Reduced_Constraint_SimOpt<Real>::applyAdjointHessian( Vector<Real> &ahwv, const Vector<Real> &w,
                          const Vector<Real> &v, const Vector<Real> &z,
                          Real &tol ) {
  nhess_++;
  if ( useFDhessVec_ ) {
    Constraint<Real>::applyAdjointHessian(ahwv,w,v,z,tol);
  }
  else {
    // Solve state equation
    solve_state_equation(z,tol);
    // Solve adjoint equation
    solve_adjoint_equation(w,z,tol);
    // Solve state sensitivity equation
    solve_state_sensitivity(v,z,tol);
    // Solve adjoint sensitivity equation
    solve_adjoint_sensitivity(w,v,z,tol);
    // Build hessVec
    conRed_->applyAdjointJacobian_2(ahwv,*adjoint_sens_,*state_,z,tol);
    conVal_->applyAdjointHessian_12(*dualcontrol_,w,*state_sens_,*state_,z,tol);
    ahwv.plus(*dualcontrol_);
    conVal_->applyAdjointHessian_22(*dualcontrol_,w,v,*state_,z,tol);
    ahwv.plus(*dualcontrol_);
    conRed_->applyAdjointHessian_12(*dualcontrol_,*adjoint_,*state_sens_,*state_,z,tol);
    ahwv.plus(*dualcontrol_);
    conRed_->applyAdjointHessian_22(*dualcontrol_,*adjoint_,v,*state_,z,tol);
    ahwv.plus(*dualcontrol_);
  }
}

}

#endif
