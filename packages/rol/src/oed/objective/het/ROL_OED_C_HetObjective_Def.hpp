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


#ifndef ROL_OED_C_HET_OBJECTIVE_DEF_HPP
#define ROL_OED_C_HET_OBJECTIVE_DEF_HPP

namespace ROL {
namespace OED {
namespace Het {

template<typename Real>
C_Objective<Real>::C_Objective( const Ptr<BilinearConstraint<Real>>   &con,
             const Ptr<QuadraticObjective<Real>> &obj,
             const Ptr<Vector<Real>>             &state,
             bool storage) {
  ObjectiveBase<Real,std::vector<Real>>::setConstraint(con);
  ObjectiveBase<Real,std::vector<Real>>::setObjective(obj);
  ObjectiveBase<Real,std::vector<Real>>::setStorage(storage);
  ObjectiveBase<Real,std::vector<Real>>::initialize(state);
}

template<typename Real>
Real C_Objective<Real>::value( const Vector<Real> &z, Real &tol ) {
  // Solve state equation
  ObjectiveBase<Real,std::vector<Real>>::solve_state_equation(param_,z,tol);
  Vector<Real> &state = *ObjectiveBase<Real,std::vector<Real>>::getState();
  // Get objective function value
  return ObjectiveBase<Real,std::vector<Real>>::getObjective()->value(state,z,tol);
}

template<typename Real>
void C_Objective<Real>::gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
  if (dp_ == nullPtr) dp_ = g.clone();  
  // Solve state equation
  ObjectiveBase<Real,std::vector<Real>>::solve_state_equation(param_,z,tol);
  Vector<Real> &state = *ObjectiveBase<Real,std::vector<Real>>::getState();
  // Solve adjoint equation
  ObjectiveBase<Real,std::vector<Real>>::solve_adjoint_equation(param_,z,tol);
  Vector<Real> &adjoint = *ObjectiveBase<Real,std::vector<Real>>::getAdjoint();
  // Build gradient
  ObjectiveBase<Real,std::vector<Real>>::getObjective()->gradient_2(*dp_,state,z,tol);
  ObjectiveBase<Real,std::vector<Real>>::getConstraint()->applyAdjointJacobian_2(g,adjoint,state,z,tol);
  g.plus(*dp_);
}

template<typename Real>
void C_Objective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  if (dp_ == nullPtr) dp_ = hv.clone();  
  // Solve state equation
  ObjectiveBase<Real,std::vector<Real>>::solve_state_equation(param_,z,tol);
  Vector<Real> &state = *ObjectiveBase<Real,std::vector<Real>>::getState();
  // Solve adjoint equation
  ObjectiveBase<Real,std::vector<Real>>::solve_adjoint_equation(param_,z,tol);
  Vector<Real> &adjoint = *ObjectiveBase<Real,std::vector<Real>>::getAdjoint();
  // Solve state sensitivity equation
  ObjectiveBase<Real,std::vector<Real>>::solve_state_sensitivity(v,z,tol);
  Vector<Real> &stateSens = *ObjectiveBase<Real,std::vector<Real>>::getStateSens();
  // Solve adjoint sensitivity equation
  ObjectiveBase<Real,std::vector<Real>>::solve_adjoint_sensitivity(v,z,tol);
  Vector<Real> &adjointSens = *ObjectiveBase<Real,std::vector<Real>>::getAdjointSens();
  // Build hessVec
  ObjectiveBase<Real,std::vector<Real>>::getConstraint()->applyAdjointJacobian_2(hv,adjointSens,state,z,tol);
  ObjectiveBase<Real,std::vector<Real>>::getObjective()->hessVec_21(*dp_,stateSens,state,z,tol);
  hv.plus(*dp_);
  ObjectiveBase<Real,std::vector<Real>>::getConstraint()->applyAdjointHessian_12(*dp_,adjoint,stateSens,state,z,tol);
  hv.plus(*dp_);
}

} // END Het Namespace
} // END OED Namespace
} // END ROL Namespace

#endif
