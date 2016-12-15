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

#include "ROL_EqualityConstraint_SimOpt.hpp"
#include "fractional_operator.hpp"
#include "fractional_vector.hpp"

template <class Real>
class FractionalConstraint : public ROL::EqualityConstraint_SimOpt<Real> {
private:
  const Teuchos::RCP<PDE<Real> > pde_local_;
  const Teuchos::RCP<PDE<Real> > pde_cylinder_;

  Teuchos::RCP<Assembler<Real> > assembler_local_;
  Teuchos::RCP<Assembler<Real> > assembler_cylinder_;

  Teuchos::RCP<Tpetra::CrsMatrix<> > Klocal_, Mlocal_;
  Teuchos::RCP<Tpetra::CrsMatrix<> > Kcylinder_, Mcylinder_;
  Teuchos::RCP<Tpetra::MultiVector<> > Flocal_;

  Teuchos::RCP<Tpetra::MultiVector<> > ulocal_;
  Teuchos::RCP<Tpetra::MultiVector<> > zlocal_;
  Teuchos::RCP<Tpetra::MultiVector<> > ucylinder_;
  Teuchos::RCP<Tpetra::MultiVector<> > zcylinder_;

  void assemble(void) {
    // Assemble local components
    assembler_local_->assemblePDEJacobian1(Klocal_,pde_local_,ulocal_,zlocal_);
    assembler_local_->assemblePDERieszMap1(Mlocal_,pde_local_);
    assembler_local_->assemblePDEResidual(Flocal_,pde_local_,ulocal_,zlocal_);
    // Assemble cylinder components
    assembler_cylinder_->assemblePDEJacobian1(Kcylinder_,pde_cylinder_,ucylinder_,zcylinder_);
    assembler_cylinder_->assemblePDERieszMap1(Mcylinder_,pde_cylinder_);
  }

public:
  LinearFractionalConstraint(const Teuchos::RCP<PDE<Real> >                & pde_local,
                             const Teuchos::RCP<MeshManager<Real> >        & mesh_local,
                             const Teuchos::RCP<const Teuchos::Comm<int> > & comm_local,
                             const Teuchos::RCP<PDE<Real> >                & pde_cylinder,
                             const Teuchos::RCP<MeshManager<Real> >        & mesh_cylinder,
                             const Teuchos::RCP<const Teuchos::Comm<int> > & comm_cylinder,
                             Teuchos::ParameterList                        & parlist,
                             std::ostream                                  & outStream = std::cout)
    : pde_local_(pde_local), pde_cylinder_(pde_cylinder) {
    assembler_local_ = Teuchos::rcp(new Assembler<Real>(pde_local_->getFields(),
                                                        mesh_local,
                                                        comm_local,
                                                        parlist,
                                                        outStream));
    assembler_local_->setCellNodes(*pde_local_);
    ulocal_ = assembler_local_->createStateVector();
    zlocal_ = assembler_local_->createControlVector();
    assembler_cylinder_ = Teuchos::rcp(new Assembler<Real>(pde_cylinder_->getFields(),
                                                           mesh_cylinder,
                                                           comm_cylinder,
                                                           parlist,
                                                           outStream));
    assembler_cylinder_->setCellNodes(*pde_cylinder_);
    ucylinder_ = assembler_cylinder_->createStateVector();
    zcylinder_ = assembler_cylinder_->createControlVector();
  }

  void setParameter(const std::vector<Real> &param) {
    ROL::EqualityConstraint_SimOpt<Real>::setParameter(param);
    pde_local_->setParameter(param);
    pde_cylinder_->setParameter(param);
  }

  void apply( ROL::Vector<Real> &Hv, const ROL::Vector<Real> &v, Real &tol ) const {
    Teuchos::RCP<Tpetra::MultiVector<> >       Hvf = getField(Hv);
    Teuchos::RCP<const Tpetra::MultiVector<> >  vf = getConstField(v);

    Teuchos::ArrayView<const int> indices;
    Teuchos::ArrayView<const Real> values;
    for (size_t r = 0; r < Mcylinder_->getGlobalNumRows(); ++r) {
      Mcylinder_->getGlobalRowView(r,indices,values);
      for (int c = 0; c < indices.size(); ++c) {
        Klocal_->apply(*(vf->getVector(indices[c])),*(Hvf->getVectorNonConst(r)),Teuchos::NO_TRANS,values[c],static_cast<Real>(0));
      }
    }
    for (size_t r = 0; r < Kcylinder_->getGlobalNumRows(); ++r) {
      Kcylinder_->getGlobalRowView(r,indices,values);
      for (int c = 0; c < indices.size(); ++c) {
        Mlocal_->apply(*(vf->getVector(indices[c])),*(Hvf->getVectorNonConst(r)),Teuchos::NO_TRANS,values[c],static_cast<Real>(1));
      }
    }
  }

private: // Vector accessor functions

  Teuchos::RCP<const Tpetra::MultiVector<> > getConstField(const ROL::Vector<Real> &x) const {
    return Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(x).getVector();
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getField(ROL::Vector<Real> &x) const {
    return Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(x).getVector();
  }
};

template <class Real>
class FractionalPreconditioner : public ROL::LinearOperator<Real> {
private:
  const Teuchos::RCP<PDE<Real> > pde_local_;
  const Teuchos::RCP<PDE<Real> > pde_cylinder_;
  Teuchos::RCP<Assembler<Real> > assembler_local_;
  Teuchos::RCP<Assembler<Real> > assembler_cylinder_;

  Teuchos::RCP<Tpetra::CrsMatrix<> > Klocal_, Mlocal_;
  Teuchos::RCP<Tpetra::CrsMatrix<> > Kcylinder_, Mcylinder_;

public:
  FractionalPreconditioner(const Teuchos::RCP<PDE<Real> > &pde_local,
                           const Teuchos::RCP<MeshManager<Real> > &mesh_local,
                           const Teuchos::RCP<PDE<Real> > &pde_cylinder,
                           const Teuchos::RCP<MeshManager<Real> > &mesh_cylinder,
                           const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                           Teuchos::ParameterList &parlist,
                           std::ostream &outStream = std::cout)
    : pde_local_(pde_local), pde_cylinder_(pde_cylinder) {
    Teuchos::RCP<Tpetra::MultiVector<> > uvec;
    Teuchos::RCP<Tpetra::MultiVector<> > zvec;
    // Assemble local components
    assembler_local_ = Teuchos::rcp(new Assembler<Real>(pde_local_->getFields(),mesh_local,comm,parlist,outStream));
    assembler_local_->setCellNodes(*pde_local_);
    uvec = assembler_local_->createStateVector();
    zvec = assembler_local_->createControlVector();
    assembler_local_->assemblePDEJacobian1(Klocal_,pde_local_,uvec,zvec);
    assembler_local_->assemblePDERieszMap1(Mlocal_,pde_local_);
    // Assemble cylinder components
    assembler_cylinder_ = Teuchos::rcp(new Assembler<Real>(pde_cylinder_->getFields(),mesh_cylinder,comm,parlist,outStream));
    assembler_cylinder_->setCellNodes(*pde_cylinder_);
    uvec = assembler_cylinder_->createStateVector();
    zvec = assembler_cylinder_->createControlVector();
    assembler_cylinder_->assemblePDEJacobian1(Kcylinder_,pde_cylinder_,uvec,zvec);
    assembler_cylinder_->assemblePDERieszMap1(Mcylinder_,pde_cylinder_);
  }

  void apply( ROL::Vector<Real> &Hv, const ROL::Vector<Real> &v, Real &tol ) const {
    Hv.set(v);
//    Teuchos::RCP<Tpetra::MultiVector<> >       Hvf = getField(v);
//    Teuchos::RCP<const Tpetra::MultiVector<> >  vf = getConstField(v);
//
//    Teuchos::ArrayView<const int> indices;
//    Teuchos::ArrayView<const Real> values;
//    for (size_t r = 0; r < Mcylinder_->getGlobalNumRows(); ++r) {
//      Mcylinder_->getGlobalRowView(r,indices,values);
//      for (int c = 0; c < indices.size(); ++c) {
//        Klocal_->apply(*(vf->getVector(indices[c])),*(Hvf->getVectorNonConst(r)),Teuchos::NO_TRANS,values[c],static_cast<Real>(0));
//      }
//    }
//    for (size_t r = 0; r < Kcylinder_->getGlobalNumRows(); ++r) {
//      Kcylinder_->getGlobalRowView(r,indices,values);
//      for (int c = 0; c < indices.size(); ++c) {
//        Mlocal_->apply(*(vf->getVector(indices[c])),*(Hvf->getVectorNonConst(r)),Teuchos::NO_TRANS,values[c],static_cast<Real>(1));
//      }
//    }
  }

private: // Vector accessor functions

  Teuchos::RCP<const Tpetra::MultiVector<> > getConstField(const ROL::Vector<Real> &x) const {
    return Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(x).getVector();
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getField(ROL::Vector<Real> &x) const {
    return Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(x).getVector();
  }
};

#endif
