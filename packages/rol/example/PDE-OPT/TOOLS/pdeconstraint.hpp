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

/*! \file  constraint.hpp
    \brief Defines the SimOpt constraint for the 'poisson' example.
*/

#ifndef ROL_PDEOPT_PDECONSTRAINT_H
#define ROL_PDEOPT_PDECONSTRAINT_H

#include "ROL_ParametrizedEqualityConstraint_SimOpt.hpp"
#include "pde.hpp"
#include "assembler.hpp"
#include "pdevector.hpp"

template<class Real>
class PDE_Constraint : public ROL::ParametrizedEqualityConstraint_SimOpt<Real> {
private:
  const Teuchos::RCP<PDE<Real> > pde_;
  Teuchos::RCP<Assembler<Real> > assembler_;

  bool computeJ1_, computeJ2_, computeJ3_;
  bool computeH11_, computeH12_, computeH21_, computeH22_;

public:
  PDE_Constraint(const Teuchos::RCP<PDE<Real> > &pde,
                 const Teuchos::RCP<MeshManager<Real> > &meshMgr,
                 const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                 Teuchos::ParameterList &parlist,
                 std::ostream &outStream = std::cout)
    : pde_(pde), computeJ1_(true), computeJ2_(true), computeJ3_(true),
      computeH11_(true), computeH12_(true), computeH21_(true), computeH22_(true) {
    assembler_ = Teuchos::rcp(new Assembler<Real>(pde_->getFields(),meshMgr,comm,parlist,outStream));
    assembler_->setCellNodes(*pde_);
  }

  PDE_Constraint(const Teuchos::RCP<PDE<Real> > &pde,
                 const Teuchos::RCP<Assembler<Real> > &assembler)
    : pde_(pde), assembler_(assembler), computeJ1_(true), computeJ2_(true), computeJ3_(true),
      computeH11_(true), computeH12_(true), computeH21_(true), computeH22_(true) {
    assembler_->setCellNodes(*pde_);
  }

  void setParameter(const std::vector<Real> &param) {
    ROL::ParametrizedEqualityConstraint_SimOpt<Real>::setParameter(param);
    pde_->setParameter(param);
  }

  const Teuchos::RCP<Assembler<Real> > getAssembler(void) const {
    return assembler_;
  }

  using ROL::EqualityConstraint_SimOpt<Real>::update_1;
  void update_1(const ROL::Vector<Real> &u, bool flag = true, int iter = -1) {
    computeJ1_ = (flag ? true : computeJ1_);
  }

  using ROL::EqualityConstraint_SimOpt<Real>::update_2;
  void update_2(const ROL::Vector<Real> &z, bool flag = true, int iter = -1) {
    computeJ2_ = (flag ? true : computeJ2_);
  }

  using ROL::EqualityConstraint_SimOpt<Real>::update;
  void update(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, bool flag = true, int iter = -1) {
    update_1(u,flag,iter);
    update_2(z,flag,iter);
    computeH11_ = (flag ? true : computeH11_);
    computeH12_ = (flag ? true : computeH12_);
    computeH21_ = (flag ? true : computeH21_);
    computeH22_ = (flag ? true : computeH22_);
  }

  using ROL::EqualityConstraint_SimOpt<Real>::value;
  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> >       cf = getField(c);
    Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
    Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
    Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

    Real one(1);
    assembler_->assemblePDEResidual(pde_,uf,zf,zp);
    cf->scale(one,*(assembler_->getPDEResidual()));
  }

  void applyJacobian_1(ROL::Vector<Real> &jv,
                 const ROL::Vector<Real> &v,
                 const ROL::Vector<Real> &u,
                 const ROL::Vector<Real> &z, Real &tol) {
    try {
      if (computeJ1_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
        Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
        Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEJacobian1(pde_,uf,zf,zp);
        computeJ1_ = false;
      }
    }
    catch ( Exception::Zero & ez ) {
     computeJ1_ = true;
    }
    catch ( Exception::NotImplemented & eni ) {
      ROL::EqualityConstraint_SimOpt<Real>::applyJacobian_1(jv,v,u,z,tol);
    }
    if ( computeJ1_ ) {
      jv.zero();
    }
    else {
      Teuchos::RCP<Tpetra::MultiVector<> >      jvf = getField(jv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      assembler_->applyPDEJacobian1(jvf,vf,false);
    }
  }


  void applyJacobian_2(ROL::Vector<Real> &jv,
                 const ROL::Vector<Real> &v,
                 const ROL::Vector<Real> &u,
                 const ROL::Vector<Real> &z, Real &tol) {
    // Apply Jacobian of field controls to vector
    int NotImplemented(0), IsZero(0);
    try {
      if (computeJ2_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
        Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
        Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEJacobian2(pde_,uf,zf,zp);
        computeJ2_ = false;
      }
    }
    catch ( Exception::Zero & ez ) {
      computeJ2_ = true;
      IsZero++;
    }
    catch ( Exception::NotImplemented & eni ) {
      computeJ2_ = true;
      NotImplemented++;
    }
    if ( !computeJ2_ ) {
      Teuchos::RCP<Tpetra::MultiVector<> >      jvf = getField(jv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      assembler_->applyPDEJacobian2(jvf,vf,false);
    }
    // Apply Jacobian of parametric controls to vector
    try {
      if (computeJ3_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
        Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
        Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEJacobian3(pde_,uf,zf,zp);
        computeJ3_ = false;
      }
    }
    catch ( Exception::Zero & ez ) {
      computeJ3_ = true;
      IsZero++;
    }
    catch ( Exception::NotImplemented & eni ) {
      computeJ3_ = true;
      NotImplemented++;
    }
    if ( !computeJ3_ ) {
      Teuchos::RCP<Tpetra::MultiVector<> >  jvf = getField(jv);
      Teuchos::RCP<const std::vector<Real> > vp = getConstParameter(v);
      assembler_->applyPDEJacobian3(jvf,vp);
    }
    // Zero Jacobian if all routines return Exception::Zero
    if ( IsZero == 2 || (IsZero == 1 && NotImplemented == 1) ) {
      jv.zero();
    }
    // Default to finite differences if all routines return Exception::NotImplemented
    if ( NotImplemented==2 ) {
      ROL::EqualityConstraint_SimOpt<Real>::applyJacobian_2(jv,v,u,z,tol);
    }    
  }


  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    try {
      if (computeJ1_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
        Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
        Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEJacobian1(pde_,uf,zf,zp);
        computeJ1_ = false;
      }
    }
    catch ( Exception::Zero & ez ) {
     computeJ1_ = true;
    }
    catch ( Exception::NotImplemented & eni ) {
      ROL::EqualityConstraint_SimOpt<Real>::applyAdjointJacobian_1(ajv,v,u,z,tol);
    }
    if ( computeJ1_ ) {
      ajv.zero();
    }
    else {
      Teuchos::RCP<Tpetra::MultiVector<> >     ajvf = getField(ajv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      assembler_->applyPDEJacobian1(ajvf,vf,true);
    }
  }


  void applyAdjointJacobian_2(ROL::Vector<Real> &ajv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    // Apply Jacobian of field controls to vector
    int NotImplemented(0), IsZero(0);
    try {
      if (computeJ2_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
        Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
        Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEJacobian2(pde_,uf,zf,zp);
        computeJ2_ = false;
      }
    }
    catch ( Exception::Zero & ez ) {
      computeJ2_ = true;
      IsZero++;
    }
    catch ( Exception::NotImplemented & eni ) {
      computeJ2_ = true;
      NotImplemented++;
    }
    if ( !computeJ2_ ) {
      Teuchos::RCP<Tpetra::MultiVector<> >     ajvf = getField(ajv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      assembler_->applyPDEJacobian2(ajvf,vf,false);
    }
    // Apply Jacobian of parametric controls to vector
    try {
      if (computeJ3_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
        Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
        Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEJacobian3(pde_,uf,zf,zp);
        computeJ3_ = false;
      }
    }
    catch ( Exception::Zero & ez ) {
      computeJ3_ = true;
      IsZero++;
    }
    catch ( Exception::NotImplemented & eni ) {
      computeJ3_ = true;
      NotImplemented++;
    }
    if ( !computeJ3_ ) {
      Teuchos::RCP<std::vector<Real> >         ajvp = getParameter(ajv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      assembler_->applyPDEAdjointJacobian3(ajvp,vf);
    }
    // Zero Jacobian if all routines return Exception::Zero
    if ( IsZero == 2 || (IsZero == 1 && NotImplemented == 1) ) {
      ajv.zero();
    }
    // Default to finite differences if all routines return Exception::NotImplemented
    if ( NotImplemented==2 ) {
      ROL::EqualityConstraint_SimOpt<Real>::applyJacobian_2(ajv,v,u,z,tol);
    }    
  }


  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    try {
      if (computeH11_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > wf = getConstField(w);
        Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
        Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
        Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEHessian11(pde_,wf,uf,zf,zp);
        computeH11_ = false;
      }
    }
    catch (Exception::Zero &ez) {
      computeH11_ = true;
    }
    catch (Exception::NotImplemented &eni) {
      ROL::EqualityConstraint_SimOpt<Real>::applyAdjointHessian_11(ahwv,w,v,u,z,tol);
      //throw Exception::NotImplemented(">>> (PDE_Constraint::applyAdjointHessian_11): Hessian not implemented.");
    }
    if ( computeH11_ ) {
      ahwv.zero();
    }
    else {
      Teuchos::RCP<Tpetra::MultiVector<> >    ahwvf = getField(ahwv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      assembler_->applyPDEHessian11(ahwvf,vf);
    }
  }


  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    try {
      if ( computeH12_ ) {
        Teuchos::RCP<const Tpetra::MultiVector<> > wf = getConstField(w);
        Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
        Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
        Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEHessian12(pde_,wf,uf,zf,zp);
        computeH12_ = false;
      }
    }
    catch (Exception::Zero &ez) {
      computeH12_ = true;
    }
    catch (Exception::NotImplemented &eni) {
      ROL::EqualityConstraint_SimOpt<Real>::applyAdjointHessian_12(ahwv,w,v,u,z,tol);
      //throw Exception::NotImplemented(">>> (PDE_Constraint::applyAdjointHessian_12): Hessian not implemented.");
    }
    if ( computeH12_ ) {
      ahwv.zero();
    }
    else {
      Teuchos::RCP<Tpetra::MultiVector<> >    ahwvf = getField(ahwv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      assembler_->applyPDEHessian12(ahwvf,vf);
    }
  }


  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    try {
      if (computeH21_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > wf = getConstField(w);
        Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
        Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
        Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEHessian21(pde_,wf,uf,zf,zp);
        computeH21_ = false;
      }
    }
    catch (Exception::Zero &ez) {
      computeH21_ = true;
    }
    catch (Exception::NotImplemented &eni) {
      ROL::EqualityConstraint_SimOpt<Real>::applyAdjointHessian_21(ahwv,w,v,u,z,tol);
    }
    if ( computeH21_ ) {
      ahwv.zero();
    }
    else {
      Teuchos::RCP<Tpetra::MultiVector<> >    ahwvf = getField(ahwv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      assembler_->applyPDEHessian21(ahwvf,vf);
    }
  }


  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    try {
      if (computeH22_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > wf = getConstField(w);
        Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
        Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
        Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEHessian22(pde_,wf,uf,zf,zp);
        computeH22_ = false;
      }
    }
    catch (Exception::Zero &ez) {
      computeH22_ = true;
    }
    catch (Exception::NotImplemented &eni) {
      ROL::EqualityConstraint_SimOpt<Real>::applyAdjointHessian_22(ahwv,w,v,u,z,tol);
    }
    if ( computeH22_ ) {
      ahwv.zero();
    }
    else {
      Teuchos::RCP<Tpetra::MultiVector<> >    ahwvf = getField(ahwv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      assembler_->applyPDEHessian22(ahwvf,vf);
    }
  }


  void applyInverseJacobian_1(ROL::Vector<Real> &ijv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    try {
      if (computeJ1_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
        Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
        Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEJacobian1(pde_,uf,zf,zp);
        computeJ1_ = false;
      }
    }
    catch ( Exception::NotImplemented & eni ) {
      ROL::EqualityConstraint_SimOpt<Real>::applyInverseJacobian_1(ijv,v,u,z,tol);
    }
    Teuchos::RCP<Tpetra::MultiVector<> >     ijvf = getField(ijv);
    Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
    assembler_->applyInverseJacobian1(ijvf,vf,false);
  }


  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &iajv,
                               const ROL::Vector<Real> &v,
                               const ROL::Vector<Real> &u,
                               const ROL::Vector<Real> &z, Real &tol) {
    try {
      if (computeJ1_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
        Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
        Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEJacobian1(pde_,uf,zf,zp);
        computeJ1_ = false;
      }
    }
    catch ( Exception::NotImplemented & eni ) {
      ROL::EqualityConstraint_SimOpt<Real>::applyInverseAdjointJacobian_1(iajv,v,u,z,tol);
    }
    Teuchos::RCP<Tpetra::MultiVector<> >    iajvf = getField(iajv);
    Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
    assembler_->applyInverseJacobian1(iajvf,vf,true);
  }

private: // Vector accessor functions

  Teuchos::RCP<const Tpetra::MultiVector<> > getConstField(const ROL::Vector<Real> &x) const {
    Teuchos::RCP<const Tpetra::MultiVector<> > xp;
    try {
      xp = Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(x).getVector();
    }
    catch (std::exception &e) {
      Teuchos::RCP<const ROL::TpetraMultiVector<Real> > xvec
        = Teuchos::dyn_cast<const PDE_OptVector<Real> >(x).getField();
      if (xvec == Teuchos::null) {
        xp = Teuchos::null;
      }
      else {
        xp = xvec->getVector();
      }
    }
    return xp;
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getField(ROL::Vector<Real> &x) const {
    Teuchos::RCP<Tpetra::MultiVector<> > xp;
    try {
      xp = Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(x).getVector();
    }
    catch (std::exception &e) {
      Teuchos::RCP<ROL::TpetraMultiVector<Real> > xvec
        = Teuchos::dyn_cast<PDE_OptVector<Real> >(x).getField();
      if ( xvec == Teuchos::null ) {
        xp = Teuchos::null;
      }
      else {
        xp = xvec->getVector();
      }
    }
    return xp;
  }

  Teuchos::RCP<const std::vector<Real> > getConstParameter(const ROL::Vector<Real> &x) const {
    Teuchos::RCP<const std::vector<Real> > xp;
    try {
      Teuchos::RCP<const ROL::StdVector<Real> > xvec
        = Teuchos::dyn_cast<const PDE_OptVector<Real> >(x).getParameter();
      if ( xvec == Teuchos::null ) {
        xp = Teuchos::null;
      }
      else {
        xp = xvec->getVector();
      }
    }
    catch (std::exception &e) {
      xp = Teuchos::null;
    }
    return xp;
  }

  Teuchos::RCP<std::vector<Real> > getParameter(ROL::Vector<Real> &x) const {
    Teuchos::RCP<std::vector<Real> > xp;
    try {
      Teuchos::RCP<ROL::StdVector<Real> > xvec
        = Teuchos::dyn_cast<PDE_OptVector<Real> >(x).getParameter();
      if ( xvec == Teuchos::null ) {
        xp = Teuchos::null;
      }
      else {
        xp = xvec->getVector();
      }
    }
    catch (std::exception &e) {
      xp = Teuchos::null;
    }
    return xp;
  }
};

#endif
