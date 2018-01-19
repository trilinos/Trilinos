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

#ifndef ROL_FRACTIONALOPERATOR_H
#define ROL_FRACTIONALOPERATOR_H

#include "ROL_LinearOperator.hpp"
#include "../../TOOLS/assembler.hpp"

// Do not instantiate the template in this translation unit.
extern template class Assembler<double>;

template <class Real>
class FractionalOperator : public ROL::LinearOperator<Real> {
private:
  ROL::Ptr<Tpetra::CrsMatrix<> > Klocal_, Mlocal_;
  ROL::Ptr<Tpetra::CrsMatrix<> > Kcylinder_, Mcylinder_;
  ROL::Ptr<Tpetra::CrsMatrix<> > KcylinderT_, McylinderT_;

  void transposeMats(void) {
    Tpetra::RowMatrixTransposer<> transposerK(Kcylinder_);
    KcylinderT_ = transposerK.createTranspose();
    Tpetra::RowMatrixTransposer<> transposerM(Mcylinder_);
    McylinderT_ = transposerM.createTranspose();
  }

  bool transpose_;

public:
  FractionalOperator(const ROL::Ptr<PDE<Real> > &pde_local,
                     const ROL::Ptr<MeshManager<Real> > &mesh_local,
                     const ROL::Ptr<PDE<Real> > &pde_cylinder,
                     const ROL::Ptr<MeshManager<Real> > &mesh_cylinder,
                     const ROL::Ptr<const Teuchos::Comm<int> > &comm,
                     Teuchos::ParameterList &parlist,
                     std::ostream &outStream = std::cout)
    : transpose_(false) {
    ROL::Ptr<Tpetra::MultiVector<> > uvec;
    ROL::Ptr<Tpetra::MultiVector<> > zvec;
    // Assemble local components
    ROL::Ptr<Assembler<Real> > assembler_local
      = ROL::makePtr<Assembler<Real>>(pde_local->getFields(),mesh_local,comm,parlist,outStream);
    assembler_local->setCellNodes(*pde_local);
    uvec = assembler_local->createStateVector();
    zvec = assembler_local->createControlVector();
    assembler_local->assemblePDEJacobian1(Klocal_,pde_local,uvec,zvec);
    assembler_local->assemblePDERieszMap1(Mlocal_,pde_local);
    // Assemble cylinder components
    ROL::Ptr<Assembler<Real> > assembler_cylinder
      = ROL::makePtr<Assembler<Real>>(pde_cylinder->getFields(),mesh_cylinder,comm,parlist,outStream);
    assembler_cylinder->setCellNodes(*pde_cylinder);
    uvec = assembler_cylinder->createStateVector();
    zvec = assembler_cylinder->createControlVector();
    assembler_cylinder->assemblePDEJacobian1(Kcylinder_,pde_cylinder,uvec,zvec);
    assembler_cylinder->assemblePDERieszMap1(Mcylinder_,pde_cylinder);
    // Transpose cylinder components
    transposeMats();
  }

  FractionalOperator(const ROL::Ptr<Tpetra::CrsMatrix<> > &Klocal,
                     const ROL::Ptr<Tpetra::CrsMatrix<> > &Mlocal,
                     const ROL::Ptr<Tpetra::CrsMatrix<> > &Kcylinder,
                     const ROL::Ptr<Tpetra::CrsMatrix<> > &Mcylinder)
    : Klocal_(Klocal), Mlocal_(Mlocal), Kcylinder_(Kcylinder), Mcylinder_(Mcylinder),
      transpose_(false) {
    // Transpose cylinder components
    transposeMats();
  }

  void setTranspose(const bool trans = true) {
    transpose_ = trans;
  }

  void apply( ROL::Vector<Real> &Hv, const ROL::Vector<Real> &v, Real &tol ) const {
    ROL::Ptr<Tpetra::MultiVector<> >       Hvf = getField(Hv);
    ROL::Ptr<const Tpetra::MultiVector<> >  vf = getConstField(v);

    if ( !transpose_ ) {
      size_t numRowEntries(0);
      Teuchos::Array<int> indices;
      Teuchos::Array<Real> values;
      Teuchos::Array<size_t> col(1), row(1);
      for (size_t r = 0; r < Mcylinder_->getGlobalNumRows(); ++r) {
        row[0] = r;
        numRowEntries = Mcylinder_->getNumEntriesInGlobalRow(r);
        indices.resize(numRowEntries); values.resize(numRowEntries);
        Mcylinder_->getGlobalRowCopy(r,indices(),values(),numRowEntries);
        for (int c = 0; c < indices.size(); ++c) {
          col[0] = static_cast<size_t>(indices[c]);
          Klocal_->apply(*(vf->subView(col())),*(Hvf->subViewNonConst(row())),Teuchos::NO_TRANS,values[c],static_cast<Real>(0));
          //Klocal_->apply(*(vf->getVector(indices[c])),*(Hvf->getVectorNonConst(r)),Teuchos::NO_TRANS,values[c],static_cast<Real>(0));
        }
      }
      for (size_t r = 0; r < Kcylinder_->getGlobalNumRows(); ++r) {
        row[0] = r;
        numRowEntries = Kcylinder_->getNumEntriesInGlobalRow(r);
        indices.resize(numRowEntries); values.resize(numRowEntries);
        Kcylinder_->getGlobalRowCopy(r,indices(),values(),numRowEntries);
        for (int c = 0; c < indices.size(); ++c) {
          col[0] = static_cast<size_t>(indices[c]);
          Mlocal_->apply(*(vf->subView(col())),*(Hvf->subViewNonConst(row())),Teuchos::NO_TRANS,values[c],static_cast<Real>(1));
          //Mlocal_->apply(*(vf->getVector(indices[c])),*(Hvf->getVectorNonConst(r)),Teuchos::NO_TRANS,values[c],static_cast<Real>(1));
        }
      }
    }
  }

private: // Vector accessor functions

  ROL::Ptr<const Tpetra::MultiVector<> > getConstField(const ROL::Vector<Real> &x) const {
    return dynamic_cast<const ROL::TpetraMultiVector<Real>&>(x).getVector();
  }

  ROL::Ptr<Tpetra::MultiVector<> > getField(ROL::Vector<Real> &x) const {
    return dynamic_cast<ROL::TpetraMultiVector<Real>&>(x).getVector();
  }
};

template <class Real>
class FractionalPreconditioner : public ROL::LinearOperator<Real> {
private:
  ROL::Ptr<Tpetra::CrsMatrix<> > Klocal_, Mlocal_;
  ROL::Ptr<Tpetra::CrsMatrix<> > Kcylinder_, Mcylinder_;
  mutable ROL::Ptr<Tpetra::CrsMatrix<> > M_;
  mutable ROL::Ptr<Solver<Real> > solver_;

  mutable Teuchos::ParameterList parlist_;

  bool transpose_;

public:
  FractionalPreconditioner(const ROL::Ptr<PDE<Real> > &pde_local,
                           const ROL::Ptr<MeshManager<Real> > &mesh_local,
                           const ROL::Ptr<PDE<Real> > &pde_cylinder,
                           const ROL::Ptr<MeshManager<Real> > &mesh_cylinder,
                           const ROL::Ptr<const Teuchos::Comm<int> > &comm,
                           Teuchos::ParameterList &parlist,
                           std::ostream &outStream = std::cout)
    : parlist_(parlist), transpose_(false) {
    ROL::Ptr<Tpetra::MultiVector<> > uvec;
    ROL::Ptr<Tpetra::MultiVector<> > zvec;
    // Assemble local components
    ROL::Ptr<Assembler<Real> > assembler_local
      = ROL::makePtr<Assembler<Real>>(pde_local->getFields(),mesh_local,comm,parlist,outStream);
    assembler_local->setCellNodes(*pde_local);
    uvec = assembler_local->createStateVector();
    zvec = assembler_local->createControlVector();
    assembler_local->assemblePDEJacobian1(Klocal_,pde_local,uvec,zvec);
    assembler_local->assemblePDERieszMap1(Mlocal_,pde_local);
    // Assemble cylinder components
    ROL::Ptr<Assembler<Real> > assembler_cylinder
      = ROL::makePtr<Assembler<Real>>(pde_cylinder->getFields(),mesh_cylinder,comm,parlist,outStream);
    assembler_cylinder->setCellNodes(*pde_cylinder);
    uvec = assembler_cylinder->createStateVector();
    zvec = assembler_cylinder->createControlVector();
    assembler_cylinder->assemblePDEJacobian1(Kcylinder_,pde_cylinder,uvec,zvec);
    assembler_cylinder->assemblePDERieszMap1(Mcylinder_,pde_cylinder);
    // Create linear solver object
    solver_ = ROL::makePtr<Solver<Real>>(parlist);
  }

  FractionalPreconditioner(const ROL::Ptr<Tpetra::CrsMatrix<> > &Klocal,
                           const ROL::Ptr<Tpetra::CrsMatrix<> > &Mlocal,
                           const ROL::Ptr<Tpetra::CrsMatrix<> > &Kcylinder,
                           const ROL::Ptr<Tpetra::CrsMatrix<> > &Mcylinder,
                           Teuchos::ParameterList                   &parlist)
    : Klocal_(Klocal), Mlocal_(Mlocal), Kcylinder_(Kcylinder), Mcylinder_(Mcylinder),
      parlist_(parlist), transpose_(false) {
    // Create linear solver object
    solver_ = ROL::makePtr<Solver<Real>>(parlist);
  }

  void setTranspose(const bool trans = true) {
    transpose_ = trans;
  }

  void apply( ROL::Vector<Real> &Hv, const ROL::Vector<Real> &v, Real &tol ) const {
    ROL::Ptr<Tpetra::MultiVector<> >       Hvf = getField(Hv);
    ROL::Ptr<const Tpetra::MultiVector<> >  vf = getConstField(v);

    size_t numRowEntries(0);
    Real massVal(0), stiffVal(0);
    Teuchos::Array<int> indices;
    Teuchos::Array<Real> values;
    Teuchos::Array<size_t> row(1);
    for (size_t r = 0; r < Mcylinder_->getGlobalNumRows(); ++r) {
      row[0] = r;
      numRowEntries = Mcylinder_->getNumEntriesInGlobalRow(r);
      indices.resize(numRowEntries); values.resize(numRowEntries);
      Mcylinder_->getGlobalRowCopy(r,indices(),values(),numRowEntries);
      massVal = static_cast<Real>(0);
      for (int c = 0; c < indices.size(); ++c) {
        if ( indices[c] == static_cast<int>(r) ) {
          massVal = values[c];
          break;
        }
      }
      numRowEntries = Kcylinder_->getNumEntriesInGlobalRow(r);
      indices.resize(numRowEntries); values.resize(numRowEntries);
      Kcylinder_->getGlobalRowCopy(r,indices(),values(),numRowEntries);
      stiffVal = static_cast<Real>(0);
      for (int c = 0; c < indices.size(); ++c) {
        if ( indices[c] == static_cast<int>(r) ) {
          stiffVal = values[c];
          break;
        }
      }
      M_ = Tpetra::MatrixMatrix::add(massVal, transpose_, *Klocal_, stiffVal, transpose_, *Mlocal_);
      M_->apply(*(vf->subView(row())),*(Hvf->subViewNonConst(row())));
      //M_->apply(*(vf->getVector(r)),*(Hvf->getVectorNonConst(r)));
    }
  }

  void applyInverse( ROL::Vector<Real> &Hv, const ROL::Vector<Real> &v, Real &tol ) const {
    ROL::Ptr<Tpetra::MultiVector<> >       Hvf = getField(Hv);
    ROL::Ptr<const Tpetra::MultiVector<> >  vf = getConstField(v);

    size_t numRowEntries(0);
    Real massVal(0), stiffVal(0);
    Teuchos::Array<int> indices;
    Teuchos::Array<Real> values;
    Teuchos::Array<size_t> row(1);
    for (size_t r = 0; r < Mcylinder_->getGlobalNumRows(); ++r) {
      row[0] = r;
      numRowEntries = Mcylinder_->getNumEntriesInGlobalRow(r);
      indices.resize(numRowEntries); values.resize(numRowEntries);
      Mcylinder_->getGlobalRowCopy(r,indices(),values(),numRowEntries);
      massVal = static_cast<Real>(0);
      for (int c = 0; c < indices.size(); ++c) {
        if ( indices[c] == static_cast<int>(r) ) {
          massVal = values[c];
          break;
        }
      }
      numRowEntries = Kcylinder_->getNumEntriesInGlobalRow(r);
      indices.resize(numRowEntries); values.resize(numRowEntries);
      Kcylinder_->getGlobalRowCopy(r,indices(),values(),numRowEntries);
      stiffVal = static_cast<Real>(0);
      for (int c = 0; c < indices.size(); ++c) {
        if ( indices[c] == static_cast<int>(r) ) {
          stiffVal = values[c];
          break;
        }
      }
      M_ = Tpetra::MatrixMatrix::add(massVal, false, *Klocal_, stiffVal, false, *Mlocal_);
      solver_ = ROL::makePtr<Solver<Real>>(parlist_);
      solver_->setA(M_);
      solver_->solve(Hvf->subViewNonConst(row()),vf->subView(row()),transpose_);
      //solver_->solve(Hvf->getVectorNonConst(r),vf->getVector(r),transpose_);
    }
  }

private: // Vector accessor functions

  ROL::Ptr<const Tpetra::MultiVector<> > getConstField(const ROL::Vector<Real> &x) const {
    return dynamic_cast<const ROL::TpetraMultiVector<Real>&>(x).getVector();
  }

  ROL::Ptr<Tpetra::MultiVector<> > getField(ROL::Vector<Real> &x) const {
    return dynamic_cast<ROL::TpetraMultiVector<Real>&>(x).getVector();
  }
};

template <class Real>
class FractionalControlOperator : public ROL::LinearOperator<Real> {
private:
  const ROL::Ptr<Tpetra::CrsMatrix<> > Blocal_;
  const int numCylinder_;

  bool transpose_;

public:
  FractionalControlOperator(const ROL::Ptr<Tpetra::CrsMatrix<> > &Blocal,
                            const int numCylinder)
    : Blocal_(Blocal), numCylinder_(numCylinder), transpose_(false) {}

  void setTranspose(const bool trans = true) {
    transpose_ = trans;
  }

  void apply( ROL::Vector<Real> &Hv, const ROL::Vector<Real> &v, Real &tol ) const {
    ROL::Ptr<Tpetra::MultiVector<> >       Hvf = getField(Hv);
    ROL::Ptr<const Tpetra::MultiVector<> >  vf = getConstField(v);

    Teuchos::Array<size_t> col(1,0);
    if ( !transpose_ ) {
      Hv.zero();
      Blocal_->apply(*vf, *(Hvf->subViewNonConst(col())));
      //Blocal_->apply(*vf, *(Hvf->getVectorNonConst(0)));
    }
    else {
      Blocal_->apply(*(vf->subView(col())), *Hvf, Teuchos::TRANS);
      //Blocal_->apply(*(vf->getVector(0)), *Hvf, Teuchos::TRANS);
    }
  }

private: // Vector accessor functions

  ROL::Ptr<const Tpetra::MultiVector<> > getConstField(const ROL::Vector<Real> &x) const {
    return dynamic_cast<const ROL::TpetraMultiVector<Real>&>(x).getVector();
  }

  ROL::Ptr<Tpetra::MultiVector<> > getField(ROL::Vector<Real> &x) const {
    return dynamic_cast<ROL::TpetraMultiVector<Real>&>(x).getVector();
  }
};

#endif
