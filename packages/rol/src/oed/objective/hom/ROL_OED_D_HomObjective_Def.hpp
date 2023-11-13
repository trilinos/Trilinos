// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package Copyright (2014)
//               Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive license
// for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the contributors may
// be used to endorse or promote products derived from this software without
// specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
// EVENT SHALL SANDIA CORPORATION OR THE CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
// THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers: Drew Kouri   (dpkouri@sandia.gov) and
// Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER


#ifndef ROL_OED_D_HOM_OBJECTIVE_DEF_HPP
#define ROL_OED_D_HOM_OBJECTIVE_DEF_HPP

namespace ROL {
namespace OED {
namespace Hom {

template<typename Real>
D_Objective<Real>::D_Objective( const Ptr<BilinearConstraint<Real>> &con,
             const Ptr<Vector<Real>>      &state,
             bool storage)
  : dim_(state->dimension()), isDetComputed_(false) {
  setConstraint(con);
  setStorage(storage);
  setUpdate(false);
  initialize(state);
}

template<typename Real>
void D_Objective<Real>::update(const Vector<Real> &z, UpdateType type, int iter) {
  ObjectiveBase<Real,std::vector<Real>>::update(z,type,iter);
  getConstraint()->update_2(z,type,iter);
  isDetComputed_ = (type!=UpdateType::Accept ? false : isDetComputed_);
}

template<typename Real>
Real D_Objective<Real>::value( const Vector<Real> &z, Real &tol ) {
  if (!isDetComputed_) {
    det_ = getConstraint()->logDeterminant(z);
    isDetComputed_ = true;
  }
  return -det_;
}

template<typename Real>
void D_Objective<Real>::gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
  if (p_ == nullPtr) p_ = g.clone();
  g.zero();
  for (int i = 0; i < dim_; ++i) {
    // Solve state equation
    std::vector<Real> param = {static_cast<Real>(i)};
    solve_state_equation(param,z,tol);
    Vector<Real> &state = *getState();
    getConstraint()->applyAdjointJacobian_2(*p_,*state.basis(i),state,z,tol);
    g.plus(*p_);
  }
  g.scale(static_cast<Real>(-1));
}

template<typename Real>
void D_Objective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  if (p_ == nullPtr) p_ = hv.clone();
  hv.zero();
  for (int i = 0; i < dim_; ++i) {
    std::vector<Real> param = {static_cast<Real>(i)};
    // Solve state equation
    solve_state_equation(param,z,tol);
    Vector<Real> &state = *getState();
    // Solve state sensitivity equation
    solve_state_sensitivity(v,z,tol);
    Vector<Real> &sens = *getStateSens();
    getConstraint()->applyAdjointJacobian_2(*p_,*state.basis(i),sens,z,tol);
    hv.plus(*p_);
  }
  hv.scale(static_cast<Real>(-1));
}

} // END Hom Namespace
} // END OED Namespace
} // END ROL Namespace

#endif
