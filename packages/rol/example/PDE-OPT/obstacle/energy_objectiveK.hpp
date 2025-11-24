// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PDEOPT_ENERGY_OBJECTIVEK_HPP
#define ROL_PDEOPT_ENERGY_OBJECTIVEK_HPP

#include "ROL_Objective.hpp"
#include "../TOOLS/assemblerK.hpp"

// Do not instantiate the template in this translation unit.
extern template class Assembler<double,Kokkos::HostSpace>;

template<class Real, class DeviceType>
class EnergyObjective : public ROL::Objective<Real> {
private:
  const ROL::Ptr<PDE<Real,DeviceType>> pde_;
  ROL::Ptr<Assembler<Real,DeviceType>> assembler_;
  bool assembleRHS_, assembleJ1_;

  ROL::Ptr<Tpetra::MultiVector<>> cvec_;
  ROL::Ptr<Tpetra::MultiVector<>> uvec_;

  ROL::Ptr<Tpetra::MultiVector<>> res_;
  ROL::Ptr<Tpetra::CrsMatrix<>>   jac_;

  void assemble(void) {
    // Assemble affine term.
    if (assembleRHS_)
      assembler_->assemblePDEResidual(res_,pde_,uvec_);
    assembleRHS_ = false;
    // Assemble jacobian_1.
    if (assembleJ1_)
      assembler_->assemblePDEJacobian1(jac_,pde_,uvec_);
    assembleJ1_ = false;
  }

public:
  EnergyObjective(const ROL::Ptr<PDE<Real,DeviceType>> &pde,
                  const ROL::Ptr<MeshManager<Real,DeviceType>> &meshMgr,
                  const ROL::Ptr<const Teuchos::Comm<int>> &comm,
                  ROL::ParameterList &parlist,
                  std::ostream &outStream = std::cout)
    : pde_(pde), assembleRHS_(true), assembleJ1_(true) {
    // Construct assembler.
    assembler_ = ROL::makePtr<Assembler<Real,DeviceType>>(pde_->getFields(),meshMgr,comm,parlist,outStream);
    assembler_->setCellNodes(*pde_);
    // Initialize zero vectors.
    cvec_ = assembler_->createResidualVector();
    uvec_ = assembler_->createStateVector();
    uvec_->putScalar(static_cast<Real>(0));
    assemble();
  }

  const ROL::Ptr<Assembler<Real,DeviceType>> getAssembler(void) const {
    return assembler_;
  }

  Real value(const ROL::Vector<Real> &u, Real &tol) {
    auto uf = getConstField(u);
    const Real half(0.5), one(1);
    jac_->apply(*uf,*cvec_);
    cvec_->update(one,*res_,half);
    Teuchos::Array<Real> val(1,0);
    cvec_->dot(*uf, val.view(0,1));
    return val[0];
  }

  void gradient(ROL::Vector<Real> &g, const ROL::Vector<Real> &u, Real &tol) {
    auto gf = getField(g);
    auto uf = getConstField(u);
    const Real one(1);
    gf->scale(one,*res_);
    jac_->apply(*uf,*cvec_);
    gf->update(one,*cvec_,one);
  }

  void hessVec(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, Real &tol) {
    auto hvf = getField(hv);
    auto vf = getConstField(v);
    auto uf = getConstField(u);
    jac_->apply(*vf,*hvf);
  }

  void precond(ROL::Vector<Real> &Pv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, Real &tol) {
    Pv.set(v.dual());
  }

private: // Vector accessor functions

  ROL::Ptr<const Tpetra::MultiVector<>> getConstField(const ROL::Vector<Real> &x) const {
    ROL::Ptr<const Tpetra::MultiVector<>> xp;
    try {
      xp = dynamic_cast<const ROL::TpetraMultiVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      ROL::Ptr<const ROL::TpetraMultiVector<Real>> xvec
        = dynamic_cast<const PDE_OptVector<Real>&>(x).getField();
      if (xvec == ROL::nullPtr) {
        xp = ROL::nullPtr;
      }
      else {
        xp = xvec->getVector();
      }
    }
    return xp;
  }

  ROL::Ptr<Tpetra::MultiVector<> > getField(ROL::Vector<Real> &x) const {
    ROL::Ptr<Tpetra::MultiVector<>> xp;
    try {
      xp = dynamic_cast<ROL::TpetraMultiVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      ROL::Ptr<ROL::TpetraMultiVector<Real>> xvec
        = dynamic_cast<PDE_OptVector<Real>&>(x).getField();
      if ( xvec == ROL::nullPtr ) {
        xp = ROL::nullPtr;
      }
      else {
        xp = xvec->getVector();
      }
    }
    return xp;
  }

};

#endif
