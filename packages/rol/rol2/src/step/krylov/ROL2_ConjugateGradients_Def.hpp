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

#pragma once
#ifndef ROL2_CONJUGATEGRADIENTS_DEF_HPP
#define ROL2_CONJUGATEGRADIENTS_DEF_HPP


namespace ROL2 {

template<class Real>
Real ConjugateGradients<Real>::run(       Vector<Real>&         x, 
                                          LinearOperator<Real>& A, 
                                    const Vector<Real>&         b, 
                                          LinearOperator<Real>& M,
                                          int&                  iter, 
                                          int&                  flag ) {
  if ( !isInitialized_ ) {
    r_  = b.clone();
    v_  = x.clone();
    p_  = x.clone();
    Ap_ = b.clone();
    isInitialized_ = true;
  }

  Real rnorm = b.norm();
  Real rtol = std::min(Krylov<Real>::getAbsoluteTolerance(),
                       Krylov<Real>::getRelativeTolerance()*rnorm);
  Real itol = std::sqrt(ROL_EPSILON<Real>);

  x.zero();
  r_->set(b);

  M.applyInverse(*v_, *r_, itol);
  p_->set(*v_);

  iter = 0;
  flag = 0;

  Real kappa(0), beta(0), alpha(0), tmp(0), zero(0);
  Real gv = v_->apply(*r_);

  for (iter = 0; iter < static_cast<int>(Krylov<Real>::getMaximumIteration()); iter++) {

    if ( useInexact_ ) itol = rtol/(Krylov<Real>::getMaximumIteration() * rnorm);

    A.apply(*Ap_, *p_, itol);
    kappa = p_->apply(*Ap_);

    if ( kappa <= zero ) {
      flag = 2;
      break;
    }

    alpha = gv/kappa;

    x.axpy(alpha,*p_);

    r_->axpy(-alpha,*Ap_);
    rnorm = r_->norm();

    if ( rnorm < rtol ) break;

    itol = default_tolerance<Real>();
    M.applyInverse(*v_, *r_, itol);
    tmp  = gv;
    gv   = v_->apply(*r_);
    beta = gv/tmp;

    p_->scale(beta);
    p_->plus(*v_);
  }

  if ( iter == (int)Krylov<Real>::getMaximumIteration() ) flag = 1;
  else iter++;
  return rnorm;
}

} // namespace ROL2 

#endif // ROL2_CONJUGATE_GRADIENTS_DEF_HPP
