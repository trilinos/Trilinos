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

#ifndef ROL_FRACTIONALCONSTRAINT_H
#define ROL_FRACTIONALCONSTRAINT_H

#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_KrylovFactory.hpp"
#include "fractional_operator.hpp"
#include "fractional_vector.hpp"

template <class Real>
class FractionalConstraint : public ROL::Constraint_SimOpt<Real> {
private:
  const ROL::Ptr<PDE<Real> > pde_local_;
  const ROL::Ptr<PDE<Real> > pde_cylinder_;

  ROL::Ptr<Assembler<Real> > assembler_local_;
  ROL::Ptr<Assembler<Real> > assembler_cylinder_;

  ROL::Ptr<Tpetra::CrsMatrix<> > Klocal_, Blocal_, Mlocal_;
  ROL::Ptr<Tpetra::CrsMatrix<> > Kcylinder_, Mcylinder_;
  ROL::Ptr<Tpetra::MultiVector<> > Flocal_;

  ROL::Ptr<ROL::Krylov<Real> > krylov_;

  ROL::Ptr<FractionalOperator<Real> >        A_;
  ROL::Ptr<FractionalControlOperator<Real> > B_;
  ROL::Ptr<FractionalPreconditioner<Real> >  M_;
  ROL::Ptr<FractionalVector<Real> >          Vec_;

  ROL::Ptr<Tpetra::MultiVector<> > ulocal_;
  ROL::Ptr<Tpetra::MultiVector<> > zlocal_;
  ROL::Ptr<Tpetra::MultiVector<> > ucylinder_;
  ROL::Ptr<Tpetra::MultiVector<> > zcylinder_;

  Teuchos::ParameterList parlist_;

  Real fracPower_;

  bool isAssembled_;

  void assemble(const ROL::Vector<Real> &z) {
    if ( !isAssembled_ ) {
      // Assemble local components
      assembler_local_->assemblePDEJacobian1(Klocal_,pde_local_,ulocal_,zlocal_);
      assembler_local_->assemblePDEJacobian2(Blocal_,pde_local_,ulocal_,zlocal_);
      assembler_local_->assemblePDERieszMap1(Mlocal_,pde_local_);
      // Assemble cylinder components
      assembler_cylinder_->assemblePDEJacobian1(Kcylinder_,pde_cylinder_,ucylinder_,zcylinder_);
      assembler_cylinder_->assemblePDERieszMap1(Mcylinder_,pde_cylinder_);
      // Create fractional operator and vector
      A_   = ROL::makePtr<FractionalOperator<Real>>(Klocal_,Mlocal_,Kcylinder_,Mcylinder_);
      B_   = ROL::makePtr<FractionalControlOperator<Real>>(Blocal_,Mcylinder_->getGlobalNumCols());
      M_   = ROL::makePtr<FractionalPreconditioner<Real>>(Klocal_,Mlocal_,Kcylinder_,Mcylinder_,parlist_);

      isAssembled_ = true;
    }
    ROL::Ptr<const Tpetra::MultiVector<> > zf
      = dynamic_cast<const ROL::TpetraMultiVector<Real>&>(z).getVector(); 
    assembler_local_->assemblePDEResidual(Flocal_,pde_local_,ulocal_,zf);
    Vec_ = ROL::makePtr<FractionalVector<Real>>(Flocal_,Klocal_->getRowMap(),Mcylinder_->getGlobalNumCols(),fracPower_);
  }

public:
  FractionalConstraint(const ROL::Ptr<PDE<Real> >                & pde_local,
                       const ROL::Ptr<MeshManager<Real> >        & mesh_local,
                       const ROL::Ptr<const Teuchos::Comm<int> > & comm_local,
                       const ROL::Ptr<PDE<Real> >                & pde_cylinder,
                       const ROL::Ptr<MeshManager<Real> >        & mesh_cylinder,
                       const ROL::Ptr<const Teuchos::Comm<int> > & comm_cylinder,
                       Teuchos::ParameterList                        & parlist,
                       std::ostream                                  & outStream = std::cout)
    : pde_local_(pde_local), pde_cylinder_(pde_cylinder), parlist_(parlist), isAssembled_(false) {
    assembler_local_ = ROL::makePtr<Assembler<Real>>(pde_local_->getFields(),
                                                        mesh_local,
                                                        comm_local,
                                                        parlist,
                                                        outStream);
    assembler_local_->setCellNodes(*pde_local_);
    ulocal_  = assembler_local_->createStateVector();
    zlocal_  = assembler_local_->createControlVector();
    assembler_cylinder_ = ROL::makePtr<Assembler<Real>>(pde_cylinder_->getFields(),
                                                           mesh_cylinder,
                                                           comm_cylinder,
                                                           parlist,
                                                           outStream);
    assembler_cylinder_->setCellNodes(*pde_cylinder_);
    ucylinder_ = assembler_cylinder_->createStateVector();
    zcylinder_ = assembler_cylinder_->createControlVector();
    krylov_ = ROL::KrylovFactory<Real>(parlist);

    fracPower_ = parlist.sublist("Problem").get("Fractional Power",0.5);
  }

  void setParameter(const std::vector<Real> &param) {
    ROL::Constraint_SimOpt<Real>::setParameter(param);
    pde_local_->setParameter(param);
    pde_cylinder_->setParameter(param);
    isAssembled_ = false;
  }

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assemble(z);
    A_->apply(c,u,tol);
    c.axpy(static_cast<Real>(-1),*(Vec_->get()));
  }

  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assemble(z);
    A_->apply(jv,v,tol);
  }

  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assemble(z);
    B_->setTranspose(false);
    B_->apply(jv,v,tol);
  }

  void applyInverseJacobian_1(ROL::Vector<Real> &ijv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assemble(z);
    int iter(0), flag(0);
    krylov_->run(ijv,*A_,v,*M_,iter,flag);
  }

  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assemble(z);
    A_->apply(ajv,v,tol);
  }
 
  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &iajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assemble(z);
    int iter(0), flag(0);
    krylov_->run(iajv,*A_,v,*M_,iter,flag);
  }

  void applyAdjointJacobian_2(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assemble(z);
    B_->setTranspose(true);
    B_->apply(ajv,v,tol);
  }

  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }

  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }

  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }

  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }

  ROL::Ptr<Assembler<Real> > getLocalAssembler(void) const {
    return assembler_local_;
  }

  ROL::Ptr<Assembler<Real> > getCylinderAssembler(void) const {
    return assembler_cylinder_;
  }
};

#endif
