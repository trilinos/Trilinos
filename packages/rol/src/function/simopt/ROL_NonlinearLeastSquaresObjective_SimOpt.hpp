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

#ifndef ROL_NONLINEARLEASTSQUARESOBJECTIVE_SIMOPT_H
#define ROL_NONLINEARLEASTSQUARESOBJECTIVE_SIMOPT_H

#include "ROL_Objective.hpp"
#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_Types.hpp"

/** @ingroup func_group
    \class ROL::NonlinearLeastSquaresObjective_SimOpt
    \brief Provides the interface to evaluate nonlinear least squares objective
           functions.

    ROL's nonlinear least squares objective function interface constructs the
    the nonlinear least squares objective function associated with the equality
    constraint \f$c(u,z)=0\f$.  That is, given \f$z\f$,
    \f[
       J(u) = \langle \mathfrak{R} c(u,z),c(u,z) \rangle_{\mathcal{C}^*,\mathcal{C}}
    \f]
    where \f$c:\mathcal{U}\times\mathcal{Z}\to\mathcal{C}\f$ and
    \f$\mathfrak{R}\in\mathcal{L}(\mathcal{C},\mathcal{C}^*)\f$ denotes the
    Riesz map from \f$\mathcal{C}\f$ into \f$\mathcal{C}^*\f$.

    ---
*/


namespace ROL {

template <class Real>
class Constraint_SimOpt;

template <class Real>
class NonlinearLeastSquaresObjective_SimOpt : public Objective<Real> {
private:
  const ROL::Ptr<Constraint_SimOpt<Real> > con_;
  const bool GaussNewtonHessian_;

  ROL::Ptr<Vector<Real> > c1_, c2_, cdual_, udual_, z_;

public:
  /** \brief Constructor. 

      This function constructs a nonlinear least squares objective function. 
      @param[in]          con   is the nonlinear equation to be solved.
      @param[in]          vec   is a constraint space vector used for cloning.
      @param[in]          GHN   is a flag dictating whether or not to use the Gauss-Newton Hessian.
  */
  NonlinearLeastSquaresObjective_SimOpt(const ROL::Ptr<Constraint_SimOpt<Real> > &con,
                                        const Vector<Real> &uvec,
                                        const Vector<Real> &zvec,
                                        const Vector<Real> &cvec,
                                        const bool GNH = false)
    : con_(con), GaussNewtonHessian_(GNH) {
    c1_ = cvec.clone();
    c2_ = cvec.clone();
    z_  = zvec.clone(); z_->set(zvec);
    cdual_ = cvec.dual().clone();
    udual_ = uvec.dual().clone();
  }

  void update( const Vector<Real> &u, bool flag = true, int iter = -1 ) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    con_->update_1(u,flag,iter);
    con_->value(*c1_,u,*z_,tol);
    cdual_->set(c1_->dual());
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    Real half(0.5);
    return half*(c1_->dot(*cdual_));
  }

  void gradient( Vector<Real> &g, const Vector<Real> &u, Real &tol ) {
    con_->applyAdjointJacobian_1(g,*cdual_,u,*z_,tol);
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, Real &tol ) {
    con_->applyJacobian_1(*c2_,v,u,*z_,tol);
    con_->applyAdjointJacobian_1(hv,c2_->dual(),u,*z_,tol);
    if ( !GaussNewtonHessian_ ) {
      con_->applyAdjointHessian_11(*udual_,*cdual_,v,u,*z_,tol);
      hv.plus(*udual_);
    }
  }

  void precond( Vector<Real> &pv, const Vector<Real> &v, const Vector<Real> &u, Real &tol ) {
    con_->applyInverseAdjointJacobian_1(*cdual_,v,u,*z_,tol);
    con_->applyInverseJacobian_1(pv,cdual_->dual(),u,*z_,tol);
  }

// Definitions for parametrized (stochastic) equality constraints
public:
  void setParameter(const std::vector<Real> &param) {
    Objective<Real>::setParameter(param);
    con_->setParameter(param);
  }
};

} // namespace ROL

#endif
