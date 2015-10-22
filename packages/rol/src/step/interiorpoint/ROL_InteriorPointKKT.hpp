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

#ifndef ROL_InteriorPointKKT_H
#define ROL_InteriorPointKKT_H

#include "ROL_DiagonalOperator.hpp"
#include "ROL_Elementwise_Function.hpp"
#include "ROL_EqualityConstraint.hpp"
#include "ROL_Objective.hpp"


namespace ROL {



/** @ingroup func_group
    \class ROL::InteriorPointKKT
    \brief KKT Hessian for an InteriorPoint as appears in eq (19.12)
           in Nocedal and Wright 2nd edition.
    \details Application of the KKT Hessian to a partitioned vector

    ---
*/

template<class Real>
class InteriorPointKKT : public LinearOpertor<Real> {

  typedef Vector<Real>             V;  
  typedef PartitionedVector<Real>  PV;
  typedef Objective<Real>          OBJ;
  typedef EqualityConstraint<Real> CON;

  typedef typename PV::size_type   size_type;

private:

  const Teuchos::RCP<OBJ> obj_;    // Objective function
  const Teuchos::RCP<CON> incon_;  // Inequality Constraint
  const Teuchos::RCP<CON> eqcon_;  // Equality Constraint

  const Teuchos::RCP<V> pv_;       // Partitioned vector storing components
  const Teuchos::RCP<V> dpv_;      // Storage partitioned vector for storing updates

  const ROL::Elementwise::Divide<Real> div_;  // Elementwise division operation

  Teuchos::RCP<const V> x_;
  Teuchos::RCP<const V> s_;
  Teuchos::RCP<const V> y_;
  Teuchos::RCP<const V> z_;
   
  Teuchos::RCP<const V> vx_;
  Teuchos::RCP<const V> vs_;
  Teuchos::RCP<const V> vy_;    
  Teuchos::RCP<const V> vz_;   

  Teuchos::RCP<V> Hvx_;
  Teuchos::RCP<V> Hvs_;
  Teuchos::RCP<V> Hvy_;  
  Teuchos::RCP<V> Hvz_; 
 
  Teuchos::RCP<V> dx_;  // Extra storage for optimization vector
  Teuchos::RCP<V> ds_;  // Extra storage for slack vector

  Teuchos::RCP<LinearOperator> Sigma_;

  Teuchos::RCP<V> sigma_;    // z/s 

  const static size_type OPT   = 0;  // Optimization vector
  const static size_type SLACK = 1;  // Slack vector
  const static size_type EQUAL = 2;  // Lagrange multipliers for equality constraint
  const static size_type INEQ  = 3;  // Lagrange multipliers for inequality constraint

public:


  InteriorPointKKT( const Teuchos::RCP<OBJ> &obj,
                    const Teuchos::RCP<CON> &incon,
                    const Teuchos::RCP<CON> &eqcon,
                    const Teuchos::RCP<const V> &ipvec
                     ) : obj_(obj), incon_(incon), eqcon_(eqcon) {

    const PV &v = dyn_cast<const PV>(*ipvec);
    dx_ = v.get(OPT).clone();
    ds_ = v.get(SLACK).clone();
  }


  void update( const V &x, bool flag = true, int iter = -1 ) {
 
    using Teuchos::RCP;
    using Teuchos::dyn_cast;

    const PV &xpv = dyn_cast<const PV>(x);

    // Access subvectors
    x_ = xpv.get(OPT);
    s_ = xpv.get(SLACK);
    y_ = xpv.get(EQUAL);
    z_ = xpv.get(INEQ);

    // Update ratio and operator
    ds_->set(*z_);
    ds_->applyBinary(div_,*s_); // sigma = z/s elementwise    
    Sigma_->update(*ds,flag,iter); 

  }


  void apply( V &Hv, const V &v, Real &tol ) {
 
    using Teuchos::RCP;
    using Teuchos::dyn_cast;

    const PV &vpv = dyn_cast<const PV>(v);
    PV &Hvpv = dyn_cast<PV>(Hv);

    // Access subvectors
    vx_ = vpv.get(OPT);
    vs_ = vpv.get(SLACK);
    vy_ = vpv.get(EQUAL);
    vz_ = vpv.get(INEQ); 
        
    Hvx_ = Hvpv.get(OPT);
    Hvs_ = Hvpv.get(SLACK);
    Hvy_ = Hvpv.get(EQUAL);
    Hvz_ = Hvpv.get(INEQ);

    // Contributions to optimization part of KKT Hess-vec
    obj_->hessVec(*Hvx_,*vx_,*x_,tol);
    eqcon_->applyAdjointHessian(*dx_,*y_,*vx_,*x_,tol);
    Hvx_->axpy(-1.0,*dx_);
    incon_->applyAdjointHessian(*dx_,*z_,*vx_,*x_,tol);
    Hvx_->axpy(-1.0,*dx_);
    eqcon_->applyAdjointJacobian(*dx_,*vy_,*x_,tol);
    Hvx_->axpy(-1.0,*dx_);
    incon_->applyAdjointJacobian(*dx_,*vz_,*x_,tol);
    Hvx_->axpy(-1.0,*dx_);

    // Contributions of slack part of KKT Hess-vec
    Sigma_->apply(*Hvs_,vs_,tol);
    Hvs_->plus(*vz_); 
  
    // Contribution to equality constraint part of KKT Hess-vec
    eqcon_->applyJacobian(*Hvy_,*vx_,*x_,tol);    

    // Contributions to inequality constraint part of the KKT Hess-vec
    incon_->applyJacobian(*Hvz_,*vx_,*x_,tol);
    Hvz_->axpy(-1.0,*vs_);

  }

}; // class InteriorPointKKT


} // namespace ROL


#endif // ROL_InteriorPointKKT_H


