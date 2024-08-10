// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde.hpp
    \brief Implements the local PDE interface for the Poisson-Boltzmann control problem.
*/

#ifndef PDE_POISSON_BOLTZMANN_HPP
#define PDE_POISSON_BOLTZMANN_HPP

#include "../TOOLS/pde.hpp"
#include "../TOOLS/fe.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real>
class PDE_Poisson_Boltzmann : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real> > cellCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real> > volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;
  std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real> > fe_vol_;

public:
  PDE_Poisson_Boltzmann(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("Problem").get("Order of FE discretization",1);
    if (basisOrder == 1) {
      basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >>();
    }
    else if (basisOrder == 2) {
      basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >>();
    }
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();                  // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                                 // create cubature factory
    int cubDegree = parlist.sublist("PDE Poisson Boltzmann").get("Cubature Degree",2); // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                                 // create default cubature
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real> > & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_vol_->N()->dimension(0);
    int f = fe_vol_->N()->dimension(1);
    int p = fe_vol_->N()->dimension(2);
    int d = cellCub_->getDimension();
    // INITIALIZE RESIDUAL
    res = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // COMPUTE STIFFNESS TERM
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *gradU_eval,
                                                  *(fe_vol_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // ADD NONLINEAR TERM
    ROL::Ptr<Intrepid::FieldContainer<Real> > valU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valU_eval)(i,j) = std::exp(-(*valU_eval)(i,j));
      }
    }
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *valU_eval,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // ADD CONTROL TERM
    if ( z_coeff != ROL::nullPtr ) {
      ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval =
        ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      fe_vol_->evaluateValue(valZ_eval, z_coeff);
      Intrepid::RealSpaceTools<Real>::scale(*valZ_eval,static_cast<Real>(-1));
      Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                    *valZ_eval,
                                                    *(fe_vol_->NdetJ()),
                                                    Intrepid::COMP_CPP, true);
    }
  }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_vol_->N()->dimension(0);
    int f = fe_vol_->N()->dimension(1);
    int p = fe_vol_->N()->dimension(2);
    // INITIALIZE JACOBIAN
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // COMPUTE STIFFNESS TERM
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *(fe_vol_->gradN()),
                                                  *(fe_vol_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // ADD NONLINEAR TERM
    ROL::Ptr<Intrepid::FieldContainer<Real> > valU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valU_eval)(i,j) = -std::exp(-(*valU_eval)(i,j));
      }
    }
    Intrepid::FieldContainer<Real> NexpU(c,f,p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(NexpU,
                                                                *valU_eval,
                                                                *(fe_vol_->N()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  NexpU,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
  }

  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    if ( z_coeff != ROL::nullPtr ) {
      // GET DIMENSIONS
      int c = fe_vol_->N()->dimension(0);
      int f = fe_vol_->N()->dimension(1);
      // INITIALIZE JACOBIAN
      jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
      // ADD CONTROL TERM
      Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                    *(fe_vol_->N()),
                                                    *(fe_vol_->NdetJ()),
                                                    Intrepid::COMP_CPP, false);
      Intrepid::RealSpaceTools<Real>::scale(*jac,static_cast<Real>(-1));
    }
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_vol_->N()->dimension(0);
    int f = fe_vol_->N()->dimension(1);
    int p = fe_vol_->N()->dimension(2);
    // INITIALIZE HESSIAN
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // COMPUTE NONLINEAR TERM
    ROL::Ptr<Intrepid::FieldContainer<Real> > valU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valL_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(valL_eval, l_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valU_eval)(i,j) = (*valL_eval)(i,j)*std::exp(-(*valU_eval)(i,j));
      }
    }
    Intrepid::FieldContainer<Real> NLexpU(c,f,p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(NLexpU,
                                                                *valU_eval,
                                                                *(fe_vol_->N()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  NLexpU,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Poisson_Boltzmann:Hessian_12: Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Poisson_Boltzmann:Hessian_21: Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Poisson_Boltzmann:Hessian_22: Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    // GET DIMENSIONS
    int c = fe_vol_->N()->dimension(0);
    int f = fe_vol_->N()->dimension(1);
    // INITIALIZE RIESZ
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    *riesz = *fe_vol_->stiffMat();
    Intrepid::RealSpaceTools<Real>::add(*riesz,*(fe_vol_->massMat()));
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    // GET DIMENSIONS
    int c = fe_vol_->N()->dimension(0);
    int f = fe_vol_->N()->dimension(1);
    // INITIALIZE RIESZ
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    *riesz = *fe_vol_->massMat();
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > getFields() {
    return basisPtrs_;
  }

  void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real> > &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds) {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fe_vol_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtr_,cellCub_);
  }

  const ROL::Ptr<FE<Real> > getFE(void) const {
    return fe_vol_;
  }

}; // PDE_Poisson_Boltzmann

#endif
