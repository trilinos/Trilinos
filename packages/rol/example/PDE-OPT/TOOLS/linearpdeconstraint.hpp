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

#ifndef ROL_PDEOPT_LINEARPDECONSTRAINT_H
#define ROL_PDEOPT_LINEARPDECONSTRAINT_H

#include "ROL_ParametrizedEqualityConstraint_SimOpt.hpp"
#include "pde.hpp"
#include "assembler.hpp"
#include "pdevector.hpp"

template<class Real>
class Linear_PDE_Constraint : public ROL::ParametrizedEqualityConstraint_SimOpt<Real> {
private:
  const Teuchos::RCP<PDE<Real> > pde_;
  Teuchos::RCP<Assembler<Real> > assembler_;
  bool initZvec_, initZpar_;
  bool assembleRHS_, assembleJ1_, assembleJ2_, assembleJ3_;

  Teuchos::RCP<Tpetra::MultiVector<> > cvec_;
  Teuchos::RCP<Tpetra::MultiVector<> > uvec_;
  Teuchos::RCP<Tpetra::MultiVector<> > zvec_;
  Teuchos::RCP<std::vector<Real> >     zpar_;

  void assemble(const Teuchos::RCP<const Tpetra::MultiVector<> > &zf,
                const Teuchos::RCP<const std::vector<Real> > &zp) {
    // Assemble affine term.
    if (assembleRHS_) {
      assembler_->assemblePDEResidual(pde_,uvec_,zvec_,zpar_);
    }
    assembleRHS_ = false;
    // Assemble jacobian_1.
    if (assembleJ1_) {
      assembler_->assemblePDEJacobian1(pde_,uvec_,zvec_,zpar_);
    }
    assembleJ1_ = false;
    // Assemble jacobian_2.
    if (assembleJ2_ && zf != Teuchos::null) {
      assembler_->assemblePDEJacobian2(pde_,uvec_,zvec_,zpar_);
    }
    assembleJ2_ = false;
    // Assemble jacobian_3.
    if (assembleJ3_ && zp != Teuchos::null) {
      assembler_->assemblePDEJacobian3(pde_,uvec_,zvec_,zpar_);
    }
    assembleJ3_ = false;
  }

public:
  Linear_PDE_Constraint(const Teuchos::RCP<PDE<Real> > &pde,
                        const Teuchos::RCP<MeshManager<Real> > &meshMgr,
                        const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                        Teuchos::ParameterList &parlist,
                        std::ostream &outStream = std::cout)
    : ROL::ParametrizedEqualityConstraint_SimOpt<Real>(), pde_(pde),
      initZvec_(true), initZpar_(true),
      assembleRHS_(true), assembleJ1_(true), assembleJ2_(true), assembleJ3_(true) {
    // Construct assembler.
    assembler_ = Teuchos::rcp(new Assembler<Real>(pde_->getFields(),meshMgr,comm,parlist,outStream));
    assembler_->setCellNodes(*pde_);
    // Initialize zero vectors.
    cvec_ = assembler_->createResidualVector();
    uvec_ = assembler_->createStateVector();
    uvec_->putScalar(static_cast<Real>(0));
  }

  void setParameter(const std::vector<Real> &param) {
    ROL::ParametrizedEqualityConstraint_SimOpt<Real>::setParameter(param);
    pde_->setParameter(param);
    assembleRHS_  = true;
    assembleJ1_   = true;
    assembleJ2_   = true;
    assembleJ3_   = true;
  }

  const Teuchos::RCP<Assembler<Real> > getAssembler(void) const {
    return assembler_;
  }

  using ROL::EqualityConstraint_SimOpt<Real>::update_1;
  void update_1(const ROL::Vector<Real> &u, bool flag = true, int iter = -1) {}

  using ROL::EqualityConstraint_SimOpt<Real>::update_2;
  void update_2(const ROL::Vector<Real> &z, bool flag = true, int iter = -1) {
    // Initialize field component of z.
    Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
    if (initZvec_ && zf != Teuchos::null) {
      zvec_ = assembler_->createControlVector();
      zvec_->putScalar(static_cast<Real>(0));
    }
    initZvec_ = false;
    // Initialize parameter component of z.
    Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);
    if (initZpar_ && zp != Teuchos::null) {
      zpar_ = Teuchos::rcp(new std::vector<Real>(zp->size(),static_cast<Real>(0)));
    }
    initZpar_ = false;
  }

  using ROL::EqualityConstraint_SimOpt<Real>::update;
  void update(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, bool flag = true, int iter = -1) {
    update_1(u,flag,iter);
    update_2(z,flag,iter);
  }

  using ROL::EqualityConstraint_SimOpt<Real>::value;
  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> >       cf = getField(c);
    Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
    Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
    Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

    assemble(zf,zp);

    const Real one(1);
    cf->scale(one,*(assembler_->getPDEResidual()));

    assembler_->applyPDEJacobian1(cvec_,uf,false);
    cf->update(one,*cvec_,one);
    if (zf != Teuchos::null) {
      assembler_->applyPDEJacobian2(cvec_,zf,false);
      cf->update(one,*cvec_,one);
    }
    if (zp != Teuchos::null) {
      assembler_->applyPDEJacobian3(cvec_,zp,false);
      cf->update(one,*cvec_,one);
    }
  }

  void solve(ROL::Vector<Real> &c, ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> >       cf = getField(c);
    Teuchos::RCP<Tpetra::MultiVector<> >       uf = getField(u);
    Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
    Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

    assemble(zf,zp);

    const Real one(1);
    cf->scale(one,*(assembler_->getPDEResidual()));
    if (zf != Teuchos::null) {
      assembler_->applyPDEJacobian2(cvec_,zf,false);
      cf->update(one,*cvec_,one);
    }
    if (zp != Teuchos::null) {
      assembler_->applyPDEJacobian3(cvec_,zp,false);
      cf->update(one,*cvec_,one);
    }

    assembler_->applyInverseJacobian1(uf,cf,false);
    uf->scale(-one);

    assembler_->applyPDEJacobian1(cvec_,uf,false);
    cf->update(one,*cvec_,one);
  }

  void applyJacobian_1(ROL::Vector<Real> &jv,
                 const ROL::Vector<Real> &v,
                 const ROL::Vector<Real> &u,
                 const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> >      jvf = getField(jv);
    Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
    Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
    Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

    assemble(zf,zp);

    assembler_->applyPDEJacobian1(jvf,vf,false);
  }


  void applyJacobian_2(ROL::Vector<Real> &jv,
                 const ROL::Vector<Real> &v,
                 const ROL::Vector<Real> &u,
                 const ROL::Vector<Real> &z, Real &tol) {
    jv.zero();
    Teuchos::RCP<Tpetra::MultiVector<> >      jvf = getField(jv);
    Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
    Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

    assemble(zf,zp);

    const Real one(1);
    Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
    if (vf != Teuchos::null) {
      assembler_->applyPDEJacobian2(cvec_,vf,false);
      jvf->update(one,*cvec_,one);
    }
    Teuchos::RCP<const std::vector<Real> >     vp = getConstParameter(v);
    if (vp != Teuchos::null) {
      assembler_->applyPDEJacobian3(cvec_,vp,false);
      jvf->update(one,*cvec_,one);
    }
  }


  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> >     ajvf = getField(ajv);
    Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
    Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
    Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

    assemble(zf,zp);

    assembler_->applyPDEJacobian1(ajvf,vf,true);
  }


  void applyAdjointJacobian_2(ROL::Vector<Real> &ajv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> >     ajvf = getField(ajv);
    Teuchos::RCP<std::vector<Real> >         ajvp = getParameter(ajv);
    Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
    Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
    Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

    assemble(zf,zp);

    if (zf != Teuchos::null) {
      assembler_->applyPDEJacobian2(ajvf,vf,true);
    }
    if (zp != Teuchos::null) {
      assembler_->applyPDEAdjointJacobian3(ajvp,vf);
    }
  }


  void applyInverseJacobian_1(ROL::Vector<Real> &ijv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> >     ijvf = getField(ijv);
    Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
    Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
    Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

    assemble(zf,zp);

    assembler_->applyInverseJacobian1(ijvf,vf,false);
  }


  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &iajv,
                               const ROL::Vector<Real> &v,
                               const ROL::Vector<Real> &u,
                               const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> >    iajvf = getField(iajv);
    Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
    Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
    Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

    assemble(zf,zp);

    assembler_->applyInverseJacobian1(iajvf,vf,true);
  }


  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }


  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }


  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }


  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
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
