// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_FRACTIONALVECTOR_H
#define ROL_FRACTIONALVECTOR_H

#include "ROL_TpetraMultiVector.hpp"
#include "../../TOOLS/assembler.hpp"
#include <math.h>

// Do not instantiate the template in this translation unit.
extern template class Assembler<double>;

template <class Real>
class FractionalVector {
private:
  const ROL::Ptr<PDE<Real> > pde_local_;
  const ROL::Ptr<PDE<Real> > pde_cylinder_;
  ROL::Ptr<Assembler<Real> > assembler_local_;
  ROL::Ptr<Assembler<Real> > assembler_cylinder_;

  ROL::Ptr<Tpetra::MultiVector<> > Flocal_;
  ROL::Ptr<Tpetra::CrsMatrix<> > Klocal_;
  ROL::Ptr<Tpetra::CrsMatrix<> > Mcylinder_;

  ROL::Ptr<Tpetra::MultiVector<> > Fptr_;
  ROL::Ptr<ROL::Vector<Real> > F_;

public:
  FractionalVector(const ROL::Ptr<PDE<Real> > &pde_local,
                   const ROL::Ptr<MeshManager<Real> > &mesh_local,
                   const ROL::Ptr<PDE<Real> > &pde_cylinder,
                   const ROL::Ptr<MeshManager<Real> > &mesh_cylinder,
                   const ROL::Ptr<const Teuchos::Comm<int> > &comm,
                   Teuchos::ParameterList &parlist,
                   std::ostream &outStream = std::cout)
    : pde_local_(pde_local), pde_cylinder_(pde_cylinder) {
    ROL::Ptr<Tpetra::MultiVector<> > uvec;
    ROL::Ptr<Tpetra::MultiVector<> > zvec;
    // Assemble local components
    assembler_local_ = ROL::makePtr<Assembler<Real>>(pde_local_->getFields(),mesh_local,comm,parlist,outStream);
    assembler_local_->setCellNodes(*pde_local_);
    uvec = assembler_local_->createStateVector();   uvec->putScalar(static_cast<Real>(0));
    zvec = assembler_local_->createControlVector(); zvec->putScalar(static_cast<Real>(0));
    assembler_local_->assemblePDEJacobian1(Klocal_,pde_local_,uvec,zvec);
    assembler_local_->assemblePDEResidual(Flocal_,pde_local_,uvec,zvec);
    // Assemble cylinder components
    assembler_cylinder_ = ROL::makePtr<Assembler<Real>>(pde_cylinder_->getFields(),mesh_cylinder,comm,parlist,outStream);
    assembler_cylinder_->setCellNodes(*pde_cylinder_);
    assembler_cylinder_->assemblePDERieszMap1(Mcylinder_,pde_cylinder_);
    // Build fractional vector
    Real s     = parlist.sublist("Problem").get("Fractional Power",0.5);
    Real alpha = static_cast<Real>(1) - static_cast<Real>(2)*s;
    Real ds    = std::pow(static_cast<Real>(2), alpha) * tgamma(static_cast<Real>(1)-s)/tgamma(s);
    Fptr_ = ROL::makePtr<Tpetra::MultiVector<>>(Klocal_->getRowMap(),Mcylinder_->getGlobalNumCols());
    Fptr_->getVectorNonConst(0)->scale(-ds,*Flocal_);
    F_ = ROL::makePtr<ROL::TpetraMultiVector<Real>>(Fptr_);
  }

  FractionalVector(const ROL::Ptr<const Tpetra::MultiVector<> > &F,
                   const ROL::Ptr<const Tpetra::Map<> > &map,
                   const int numCylinder,
                   const Real s) {
    // Build fractional vector
    Real alpha = static_cast<Real>(1) - static_cast<Real>(2)*s;
    Real ds    = std::pow(static_cast<Real>(2), alpha) * tgamma(static_cast<Real>(1)-s)/tgamma(s);
    Fptr_ = ROL::makePtr<Tpetra::MultiVector<>>(map,numCylinder);
    Fptr_->getVectorNonConst(0)->scale(-ds,*F);
    F_ = ROL::makePtr<ROL::TpetraMultiVector<Real>>(Fptr_);
  }

  const ROL::Ptr<const ROL::Vector<Real> > get(void) const {
    return F_;
  }
};

#endif
