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


#ifndef ROL_PARAMETRIZEDEQUALITYCONSTRAINT_SIMOPT_H
#define ROL_PARAMETRIZEDEQUALITYCONSTRAINT_SIMOPT_H

#include "ROL_ParametrizedEqualityConstraint.hpp"
#include "ROL_Vector_SimOpt.hpp"
#include "ROL_Types.hpp"
#include <iostream>

namespace ROL {

template <class Real>
class ParametrizedEqualityConstraint_SimOpt : public ParametrizedEqualityConstraint<Real> {
public:
  virtual void update( const Vector<Real> &u, const Vector<Real> &z, bool flag = true, int iter = -1 ) {}

  virtual void value(Vector<Real> &c,
                     const Vector<Real> &u,
                     const Vector<Real> &z,
                     Real &tol) = 0;

  virtual void solve(Vector<Real> &u, 
                     const Vector<Real> &z,
                     Real &tol) {
    u.zero();
  }

  virtual void applyJacobian_1(Vector<Real> &jv,
                               const Vector<Real> &v,
                               const Vector<Real> &u,
                               const Vector<Real> &z,
                               Real &tol) {
    Real ctol = std::sqrt(ROL_EPSILON);
    // Compute step length
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON)) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    // Update state vector to u + hv
    Teuchos::RCP<Vector<Real> > unew = u.clone();
    unew->set(u);
    unew->axpy(h,v);
    // Compute new constraint value
    this->update(*unew,z);
    this->value(jv,*unew,z,ctol);
    // Compute current constraint value
    Teuchos::RCP<Vector<Real> > cold = jv.clone();
    this->update(u,z);
    this->value(*cold,u,z,ctol);
    // Compute Newton quotient
    jv.axpy(-1.0,*cold);
    jv.scale(1.0/h);
  }

  virtual void applyJacobian_2(Vector<Real> &jv,
                               const Vector<Real> &v,
                               const Vector<Real> &u,
                               const Vector<Real> &z,
                               Real &tol) { 
    Real ctol = std::sqrt(ROL_EPSILON);
    // Compute step length
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON)) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    // Update state vector to u + hv
    Teuchos::RCP<Vector<Real> > znew = z.clone();
    znew->set(z);
    znew->axpy(h,v);
    // Compute new constraint value
    this->update(u,*znew);
    this->value(jv,u,*znew,ctol);
    // Compute current constraint value
    Teuchos::RCP<Vector<Real> > cold = jv.clone();
    this->update(u,z);
    this->value(*cold,u,z,ctol);
    // Compute Newton quotient
    jv.axpy(-1.0,*cold);
    jv.scale(1.0/h);
  }

  virtual void applyInverseJacobian_1(Vector<Real> &ijv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) {
    ijv.zero();
  }

  virtual void applyAdjointJacobian_1(Vector<Real> &ajv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) {
    Real ctol = std::sqrt(ROL_EPSILON);
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON)) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    Teuchos::RCP<Vector<Real> > cold = v.clone();
    Teuchos::RCP<Vector<Real> > cnew = v.clone();
    this->update(u,z);
    this->value(*cold,u,z,ctol);
    Teuchos::RCP<Vector<Real> > unew = u.clone();
    ajv.zero();
    for (int i = 0; i < u.dimension(); i++) {
      unew->set(u);
      unew->axpy(h,*(u.basis(i)));
      this->update(*unew,z);
      this->value(*cnew,*unew,z,ctol);
      cnew->axpy(-1.0,*cold);
      cnew->scale(1.0/h);
      ajv.axpy(cnew->dot(v),*(u.basis(i)));
    }
    this->update(u,z);
  }

  virtual void applyAdjointJacobian_2(Vector<Real> &ajv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) {
    Real ctol = std::sqrt(ROL_EPSILON);
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON)) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    Teuchos::RCP<Vector<Real> > cold = v.clone();
    Teuchos::RCP<Vector<Real> > cnew = v.clone();
    this->update(u,z);
    this->value(*cold,u,z,ctol);
    Teuchos::RCP<Vector<Real> > znew = z.clone();
    ajv.zero();
    for (int i = 0; i < z.dimension(); i++) {
      znew->set(z);
      znew->axpy(h,*(z.basis(i)));
      this->update(u,*znew);
      this->value(*cnew,u,*znew,ctol);
      cnew->axpy(-1.0,*cold);
      cnew->scale(1.0/h);
      ajv.axpy(cnew->dot(v),*(z.basis(i)));
    }
    this->update(u,z);
  }

  virtual void applyInverseAdjointJacobian_1(Vector<Real> &iajv,
                                             const Vector<Real> &v,
                                             const Vector<Real> &u,
                                             const Vector<Real> &z,
                                             Real &tol) {
    iajv.zero();
  };

  virtual void applyAdjointHessian_11(Vector<Real> &ahwv,
                                      const Vector<Real> &w,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) {
    Real jtol = std::sqrt(ROL_EPSILON);
    // Compute step size
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON)) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    // Evaluate Jacobian at new state
    Teuchos::RCP<Vector<Real> > unew = u.clone();
    unew->set(u);
    unew->axpy(h,v);
    this->update(*unew,z);
    this->applyAdjointJacobian_1(ahwv,w,*unew,z,jtol);
    // Evaluate Jacobian at old state
    Teuchos::RCP<Vector<Real> > jv = u.clone();
    this->update(u,z);
    this->applyAdjointJacobian_1(*jv,w,u,z,jtol);
    // Compute Newton quotient
    ahwv.axpy(-1.0,*jv);
    ahwv.scale(1.0/h);
  }

  virtual void applyAdjointHessian_12(Vector<Real> &ahwv,
                                      const Vector<Real> &w,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) {
    Real jtol = std::sqrt(ROL_EPSILON);
    // Compute step size
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON)) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    // Evaluate Jacobian at new state
    Teuchos::RCP<Vector<Real> > unew = u.clone();
    unew->set(u);
    unew->axpy(h,v);
    this->update(*unew,z);
    this->applyAdjointJacobian_2(ahwv,w,*unew,z,jtol);
    // Evaluate Jacobian at old state
    Teuchos::RCP<Vector<Real> > jv = z.clone();
    this->update(u,z);
    this->applyAdjointJacobian_2(*jv,w,u,z,jtol);
    // Compute Newton quotient
    ahwv.axpy(-1.0,*jv);
    ahwv.scale(1.0/h);
  }

  virtual void applyAdjointHessian_21(Vector<Real> &ahwv,
                                      const Vector<Real> &w,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) {
    Real jtol = std::sqrt(ROL_EPSILON);
    // Compute step size
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON)) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    // Evaluate Jacobian at new control
    Teuchos::RCP<Vector<Real> > znew = z.clone();
    znew->set(z);
    znew->axpy(h,v);
    this->update(u,*znew,true);
    this->applyAdjointJacobian_1(ahwv,w,u,*znew,jtol);
    // Evaluate Jacobian at old control
    Teuchos::RCP<Vector<Real> > jv = u.clone();
    this->update(u,z,true);
    this->applyAdjointJacobian_1(*jv,w,u,z,jtol);
    // Compute Newton quotient
    ahwv.axpy(-1.0,*jv);
    ahwv.scale(1.0/h);
  }

  virtual void applyAdjointHessian_22(Vector<Real> &ahwv,
                                      const Vector<Real> &w,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) {
    Real jtol = std::sqrt(ROL_EPSILON);
    // Compute step size
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON)) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    // Evaluate Jacobian at new control
    Teuchos::RCP<Vector<Real> > znew = z.clone();
    znew->set(z);
    znew->axpy(h,v);
    this->update(u,*znew);
    this->applyAdjointJacobian_2(ahwv,w,u,*znew,jtol);
    // Evaluate Jacobian at old control
    Teuchos::RCP<Vector<Real> > jv = z.clone();
    this->update(u,z);
    this->applyAdjointJacobian_2(*jv,w,u,z,jtol);
    // Compute Newton quotient
    ahwv.axpy(-1.0,*jv);
    ahwv.scale(1.0/h);
}

  virtual std::vector<Real> solveAugmentedSystem(Vector<Real> &v1,
                                                 Vector<Real> &v2,
                                                 const Vector<Real> &b1,
                                                 const Vector<Real> &b2,
                                                 const Vector<Real> &x,
                                                 Real &tol) {
    std::vector<Real> vec;
    return vec;
  }

  virtual void applyPreconditioner(Vector<Real> &pv,
                                   const Vector<Real> &v,
                                   const Vector<Real> &x,
                                   Real &tol) {
    const Vector_SimOpt<Real> &xs = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    Teuchos::RCP<Vector<Real> > ijv = (xs.get_1())->clone();
    applyInverseJacobian_1(*ijv, v, *(xs.get_1()), *(xs.get_2()), tol);
    applyInverseAdjointJacobian_1(pv, *ijv, *(xs.get_1()), *(xs.get_2()), tol);
  }

  ParametrizedEqualityConstraint_SimOpt(void) : ParametrizedEqualityConstraint<Real>() {}

  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    const Vector_SimOpt<Real> &xs = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    this->update(*(xs.get_1()),*(xs.get_2()),flag,iter);
  }

  virtual bool isFeasible( const Vector<Real> &v ) { return true; }

  virtual void value(Vector<Real> &c,
                     const Vector<Real> &x,
                     Real &tol) {
    const Vector_SimOpt<Real> &xs = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    this->value(c,*(xs.get_1()),*(xs.get_2()),tol);
  }


  virtual void applyJacobian(Vector<Real> &jv,
                             const Vector<Real> &v,
                             const Vector<Real> &x,
                             Real &tol) { 
    const Vector_SimOpt<Real> &xs = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    const Vector_SimOpt<Real> &vs = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(v));
    this->applyJacobian_1(jv,*(vs.get_1()),*(xs.get_1()),*(xs.get_2()),tol);
    Teuchos::RCP<Vector<Real> > jv2 = jv.clone();
    this->applyJacobian_2(*jv2,*(vs.get_2()),*(xs.get_1()),*(xs.get_2()),tol);
    jv.plus(*jv2);
  }


  virtual void applyAdjointJacobian(Vector<Real> &ajv,
                                    const Vector<Real> &v,
                                    const Vector<Real> &x,
                                    Real &tol) { 
    Vector_SimOpt<Real> &ajvs = Teuchos::dyn_cast<Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<Vector<Real> >(ajv));
    const Vector_SimOpt<Real> &xs = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    Teuchos::RCP<Vector<Real> > ajv1 = (ajvs.get_1())->clone();
    this->applyAdjointJacobian_1(*ajv1,v,*(xs.get_1()),*(xs.get_2()),tol);
    ajvs.set_1(*ajv1);
    Teuchos::RCP<Vector<Real> > ajv2 = (ajvs.get_2())->clone();
    this->applyAdjointJacobian_2(*ajv2,v,*(xs.get_1()),*(xs.get_2()),tol);
    ajvs.set_2(*ajv2);
  }


  virtual void applyAdjointHessian(Vector<Real> &ahwv,
                                   const Vector<Real> &w,
                                   const Vector<Real> &v,
                                   const Vector<Real> &x,
                                   Real &tol) {
    Vector_SimOpt<Real> &ahwvs = Teuchos::dyn_cast<Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<Vector<Real> >(ahwv));
    const Vector_SimOpt<Real> &xs = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    const Vector_SimOpt<Real> &vs = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(v));
    // Block-row 1
    Teuchos::RCP<Vector<Real> > C11 = (ahwvs.get_1())->clone();
    Teuchos::RCP<Vector<Real> > C21 = (ahwvs.get_1())->clone();
    this->applyAdjointHessian_11(*C11,w,*(vs.get_1()),*(xs.get_1()),*(xs.get_2()),tol);
    this->applyAdjointHessian_21(*C21,w,*(vs.get_2()),*(xs.get_1()),*(xs.get_2()),tol);
    C11->plus(*C21);
    ahwvs.set_1(*C11); 
    // Block-row 2
    Teuchos::RCP<Vector<Real> > C12 = (ahwvs.get_2())->clone();
    Teuchos::RCP<Vector<Real> > C22 = (ahwvs.get_2())->clone();
    this->applyAdjointHessian_12(*C12,w,*(vs.get_1()),*(xs.get_1()),*(xs.get_2()),tol);
    this->applyAdjointHessian_22(*C22,w,*(vs.get_2()),*(xs.get_1()),*(xs.get_2()),tol);
    C22->plus(*C12);
    ahwvs.set_2(*C22); 
  }

  virtual Real checkSolve(const Vector<Real> &u, 
                          const Vector<Real> &z, 
                          const bool printToScreen = true) {
    // Solve equality constraint for u. 
    Real tol = ROL_EPSILON;
    Teuchos::RCP<Vector<Real> > s = u.clone();
    solve(*s,z,tol);
    // Evaluate equality constraint residual at (u,z).
    Teuchos::RCP<Vector<Real> > c = u.clone();
    value(*c,*s,z,tol);
    // Output norm of residual.
    Real cnorm = c->norm();
    if ( printToScreen ) {
      std::stringstream hist;
      hist << std::scientific << std::setprecision(8);
      hist << "\nTest SimOpt solve at feasible (u,z): \n  ||c(u,z)|| = " << cnorm << "\n";
      std::cout << hist.str();
    }
    return cnorm;
  }

  virtual Real checkJacobian_1(const Vector<Real> &w, 
                               const Vector<Real> &v, 
                               const Vector<Real> &u,
                               const Vector<Real> &z,
                               const bool printToScreen = true) {
    Real tol = ROL_EPSILON;
    Teuchos::RCP<Vector<Real> > Jv = w.clone();
    applyJacobian_1(*Jv,v,u,z,tol);
    Real wJv = w.dot(*Jv);
    Teuchos::RCP<Vector<Real> > Jw = v.clone();
    applyAdjointJacobian_1(*Jw,w,u,z,tol);
    Real vJw = v.dot(*Jw);
    Real diff = std::abs(wJv-vJw);
    if ( printToScreen ) {
      std::stringstream hist;
      hist << std::scientific << std::setprecision(8);
      hist << "\nTest SimOpt consistency of Jacobian_1 and its adjoint: \n  |<w,Jv> - <adj(J)w,v>| = " 
           << diff << "\n";
      hist << "  |<w,Jv>|               = " << std::abs(wJv) << "\n";
      hist << "  Relative Error         = " << diff / (std::abs(wJv)+ROL_UNDERFLOW) << "\n";
      std::cout << hist.str();
    }
    return diff;
  }

  virtual Real checkJacobian_2(const Vector<Real> &w, 
                               const Vector<Real> &v, 
                               const Vector<Real> &u,
                               const Vector<Real> &z,
                               const bool printToScreen = true) {
    Real tol = ROL_EPSILON;
    Teuchos::RCP<Vector<Real> > Jv = w.clone();
    applyJacobian_2(*Jv,v,u,z,tol);
    Real wJv = w.dot(*Jv);
    Teuchos::RCP<Vector<Real> > Jw = v.clone();
    applyAdjointJacobian_2(*Jw,w,u,z,tol);
    Real vJw = v.dot(*Jw);
    Real diff = std::abs(wJv-vJw);
    if ( printToScreen ) {
      std::stringstream hist;
      hist << std::scientific << std::setprecision(8);
      hist << "\nTest SimOpt consistency of Jacobian_2 and its adjoint: \n  |<w,Jv> - <adj(J)w,v>| = "
           << diff << "\n";
      hist << "  |<w,Jv>|               = " << std::abs(wJv) << "\n";
      hist << "  Relative Error         = " << diff / (std::abs(wJv)+ROL_UNDERFLOW) << "\n";
      std::cout << hist.str();
    }
    return diff;
  }

  virtual Real checkInverseJacobian_1(const Vector<Real> &jv, 
                                      const Vector<Real> &v, 
                                      const Vector<Real> &u, 
                                      const Vector<Real> &z, 
                                      const bool printToScreen = true) {
    Real tol = ROL_EPSILON;
    Teuchos::RCP<Vector<Real> > Jv = jv.clone();
    applyJacobian_1(*Jv,v,u,z,tol);
    Teuchos::RCP<Vector<Real> > iJJv = u.clone();
    applyInverseJacobian_1(*iJJv,*Jv,u,z,tol);
    Teuchos::RCP<Vector<Real> > diff = v.clone();
    diff->set(v);
    diff->axpy(-1.0,*iJJv);
    Real dnorm = diff->norm();
    if ( printToScreen ) {
      std::stringstream hist;
      hist << std::scientific << std::setprecision(8);
      hist << "\nTest SimOpt consistency of inverse Jacobian_1: \n  ||v-inv(J)Jv|| = " 
           << dnorm << "\n";
      hist << "  ||v||          = " << v.norm() << "\n";
      hist << "  Relative Error = " << dnorm / (v.norm()+ROL_UNDERFLOW) << "\n";
      std::cout << hist.str();
    }
    return dnorm;
  }

  virtual Real checkInverseAdjointJacobian_1(const Vector<Real> &jv, 
                                             const Vector<Real> &v, 
                                             const Vector<Real> &u, 
                                             const Vector<Real> &z, 
                                             const bool printToScreen = true) {
    Real tol = ROL_EPSILON;
    Teuchos::RCP<Vector<Real> > Jv = jv.clone();
    applyAdjointJacobian_1(*Jv,v,u,z,tol);
    Teuchos::RCP<Vector<Real> > iJJv = v.clone();
    applyInverseAdjointJacobian_1(*iJJv,*Jv,u,z,tol);
    Teuchos::RCP<Vector<Real> > diff = v.clone();
    diff->set(v);
    diff->axpy(-1.0,*iJJv);
    Real dnorm = diff->norm();
    if ( printToScreen ) {
      std::stringstream hist;
      hist << std::scientific << std::setprecision(8);
      hist << "\nTest SimOpt consistency of inverse adjoint Jacobian_1: \n  ||v-inv(adj(J))adj(J)v|| = "
           << dnorm << "\n";
      hist << "  ||v||                   = " << v.norm() << "\n";
      hist << "  Relative Error          = " << dnorm / (v.norm()+ROL_UNDERFLOW) << "\n";
      std::cout << hist.str();
    }
    return dnorm;
  }

}; // class ParametrizedEqualityConstraint_SimOpt

} // namespace ROL

#endif
