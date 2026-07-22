// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde.hpp
    \brief Implements the local PDE interface for a semilinear control problem.
*/

#ifndef PDE_SEMILINEARK_HPP
#define PDE_SEMILINEARK_HPP

#include "../TOOLS/pdeK.hpp"
#include "../TOOLS/feK.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real, class DeviceType>
class PDE_Semilinear : public PDE<Real, DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type = FE<Real, DeviceType>;
  using fst = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct = Intrepid2::CellTools<DeviceType>;
  using rst = Intrepid2::RealSpaceTools<DeviceType>;

private:
  // Finite element basis information
  basis_ptr basisPtr_;
  std::vector<basis_ptr> basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid2::Cubature<DeviceType, Real, Real>> cellCub_;
  // Cell node information
  scalar_view volCellNodes_;
  std::vector<std::vector<scalar_view>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<fe_type> fe_vol_;
  std::vector<std::vector<int>> fidx_;
  std::vector<std::vector<scalar_view>> bdryCellDofValues_;

  Real DirichletFunc(const std::vector<Real> & coords, int sideset, int locSideId) const {
    return 0;
  }

public:
  PDE_Semilinear(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("Problem").get("Order of FE discretization",1);
    if (basisOrder == 1) {
      basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType, Real, Real>>();
    }
    else if (basisOrder == 2) {
      basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<DeviceType, Real, Real>>();
    }
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();                  // get the cell type from any basis
    Intrepid2::DefaultCubatureFactory cubFactory;                                      // create cubature factory
    int cubDegree = parlist.sublist("PDE Semilinear").get("Cubature Degree",2); // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create<DeviceType, Real, Real>(cellType, cubDegree);         // create default cubature
  }

  void residual(scalar_view & res,
                const scalar_view u_coeff,
                const scalar_view z_coeff = scalar_view(),
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // GET DIMENSIONS
    int c = fe_vol_->N().extent_int(0);
    int f = fe_vol_->N().extent_int(1);
    int p = fe_vol_->N().extent_int(2);
    int d = cellCub_->getDimension();
    // INITIALIZE RESIDUAL
    res = scalar_view("res", c, f);
    // COMPUTE STIFFNESS TERM
    scalar_view gradU_eval("gradU_eval", c, p, d);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    fst::integrate(res, gradU_eval, fe_vol_->gradNdetJ(), false);
    // ADD NONLINEAR TERM
    scalar_view valU_eval("valU_eval", c, p);
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        // (valU_eval)(i,j) = std::exp(-(valU_eval)(i,j));
        (valU_eval)(i,j) = std::pow((valU_eval)(i,j),3);
      }
    }
    fst::integrate(res, valU_eval, fe_vol_->NdetJ(), true);
    // ADD CONTROL TERM
    if ( z_coeff != scalar_view() ) {
      scalar_view valZ_eval("valZ_eval", c, p);
      fe_vol_->evaluateValue(valZ_eval, z_coeff);
      rst::scale(valZ_eval, static_cast<Real>(-1));
      fst::integrate(res, valZ_eval, fe_vol_->NdetJ(), true);
    }
    // =====================================
    // DIRICHLET CONDITIONS
    // =====================================
    int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      int numLocalSideIds = bdryCellLocIds_[i].size();
      for (int j = 0; j < numLocalSideIds; ++j) {
        int numCellsSide = bdryCellLocIds_[i][j].size();
        int numVBdryDofs = fidx_[j].size();
        for (int k = 0; k < numCellsSide; ++k) {
          int cidx = bdryCellLocIds_[i][j][k];
          for (int l = 0; l < numVBdryDofs; ++l) {
            for (int m = 0; m < d; ++m) {
              res(cidx,fidx_[j][l]) = u_coeff(cidx,fidx_[j][l]) - (bdryCellDofValues_[i][j])(k,fidx_[j][l]);
            }
          }
        }
      }
    }
  }

  void Jacobian_1(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    // GET DIMENSIONS
    int c = fe_vol_->N().extent_int(0);
    int f = fe_vol_->N().extent_int(1);
    int p = fe_vol_->N().extent_int(2);
    // INITIALIZE JACOBIAN
    jac = scalar_view("jac1", c, f, f);
    // COMPUTE STIFFNESS TERM
    fst::integrate(jac, fe_vol_->gradN(), fe_vol_->gradNdetJ(), false);
    // ADD NONLINEAR TERM
    scalar_view valU_eval("valU_eval", c, p);
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (valU_eval)(i,j) = 3*std::pow(valU_eval(i,j),2);
      }
    }
    scalar_view NexpU("NexpU", c, f, p);
    fst::scalarMultiplyDataField(NexpU, valU_eval, fe_vol_->N());
    fst::integrate(jac, NexpU, fe_vol_->NdetJ(), true);
    // =====================================
    // DIRICHLET CONDITIONS
    // =====================================
    int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      int numLocalSideIds = bdryCellLocIds_[i].size();
      for (int j = 0; j < numLocalSideIds; ++j) {
        int numCellsSide = bdryCellLocIds_[i][j].size();
        int numBdryDofs = fidx_[j].size();
        for (int k = 0; k < numCellsSide; ++k) {
          int cidx = bdryCellLocIds_[i][j][k];
          for (int l = 0; l < numBdryDofs; ++l) {
            for (int m = 0; m < f; ++m) {
              jac(cidx,fidx_[j][l],m) = static_cast<Real>(0);
            }
            jac(cidx,fidx_[j][l],fidx_[j][l]) = static_cast<Real>(1);
          }
        }
      }
    }
  }

  void Jacobian_2(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    if ( z_coeff != scalar_view() ) {
      // GET DIMENSIONS
      int c = fe_vol_->N().extent_int(0);
      int f = fe_vol_->N().extent_int(1);
      // INITIALIZE JACOBIAN
      jac = scalar_view("jac2", c, f, f);
      // ADD CONTROL TERM
      fst::integrate(jac, fe_vol_->N(), fe_vol_->NdetJ(), false);
      rst::scale(jac, static_cast<Real>(-1));
      // =====================================
      // DIRICHLET CONDITIONS
      // =====================================
      int numSideSets = bdryCellLocIds_.size();
      if (numSideSets > 0) {
        for (int i = 0; i < numSideSets; ++i) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m = 0; m < f; ++m) {
                  jac(cidx,fidx_[j][l],m) = static_cast<Real>(0);
                }
              }
            }
          }
        }
      }
    }
  }

  void Hessian_11(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    // GET DIMENSIONS
    int c = fe_vol_->N().extent_int(0);
    int f = fe_vol_->N().extent_int(1);
    int p = fe_vol_->N().extent_int(2);
    // INITIALIZE HESSIAN
    hess = scalar_view("hess11", c, f, f);

    // =====================================
    // DIRICHLET CONDITIONS
    // =====================================
    int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      int numLocalSideIds = bdryCellLocIds_[i].size();
      for (int j = 0; j < numLocalSideIds; ++j) {
        int numCellsSide = bdryCellLocIds_[i][j].size();
        int numVBdryDofs = fidx_[j].size();
        for (int k = 0; k < numCellsSide; ++k) {
          int cidx = bdryCellLocIds_[i][j][k];
          for (int l = 0; l < numVBdryDofs; ++l) {
            l_coeff(cidx,fidx_[j][l]) = static_cast<Real>(0);
          }
        }
      }
    }
    // =====================================

    // COMPUTE NONLINEAR TERM
    scalar_view valU_eval("valU_eval", c, p);
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    scalar_view valL_eval("valL_eval", c, p);
    fe_vol_->evaluateValue(valL_eval, l_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        // (valU_eval)(i,j) = (valL_eval)(i,j)*std::exp(-(valU_eval)(i,j));
        (valU_eval)(i,j) = (valL_eval)(i,j)*6*(valU_eval)(i,j);
      }
    }
    scalar_view NLexpU("NLexpU", c, f, p);
    fst::scalarMultiplyDataField(NLexpU, valU_eval, fe_vol_->N());
    fst::integrate(hess, NLexpU, fe_vol_->NdetJ(), false);
  }

  void Hessian_12(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Semilinear:Hessian_12: Hessian is zero.");
  }

  void Hessian_21(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Semilinear:Hessian_21: Hessian is zero.");
  }

  void Hessian_22(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Semilinear:Hessian_22: Hessian is zero.");
  }

  void RieszMap_1(scalar_view & riesz) override {
    // GET DIMENSIONS
    int c = fe_vol_->N().extent_int(0);
    int f = fe_vol_->N().extent_int(1);
    // INITIALIZE RIESZ
    riesz = scalar_view("riesz1", c, f, f);
    Kokkos::deep_copy(riesz, fe_vol_->stiffMat());
    rst::add(riesz,fe_vol_->massMat());
  }

  void RieszMap_2(scalar_view & riesz) override {
    // GET DIMENSIONS
    int c = fe_vol_->N().extent_int(0);
    int f = fe_vol_->N().extent_int(1);
    // INITIALIZE RIESZ
    riesz = scalar_view("riesz2", c, f, f);
    Kokkos::deep_copy(riesz, fe_vol_->massMat());
  }

  std::vector<basis_ptr> getFields() override {
    return basisPtrs_;
  }

  void setCellNodes(const scalar_view & volCellNodes,
                    const std::vector<std::vector<scalar_view>> & bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> & bdryCellLocIds) override {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fe_vol_ = ROL::makePtr<fe_type>(volCellNodes_, basisPtr_, cellCub_);
    fidx_ = fe_vol_->getBoundaryDofs();
    // Compute Dirichlet values at DOFs.
    int d = basisPtr_->getBaseCellTopology().getDimension();
    int numSidesets = bdryCellLocIds_.size();
    bdryCellDofValues_.resize(numSidesets);
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellDofValues_[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
        int f = basisPtr_->getCardinality();
        bdryCellDofValues_[i][j] = scalar_view("bdryVals", c, f);
        scalar_view coords("coords", c, f, d);
        if (c > 0) {
          fe_vol_->computeDofCoords(coords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<f; ++l) {
            std::vector<Real> dofpoint(d);
            for (int m=0; m<d; ++m) {
              dofpoint[m] = coords(k, l, m);
            }
            bdryCellDofValues_[i][j](k, l) = DirichletFunc(dofpoint, i, j);
          }
        }
      }
    }
  }

  const ROL::Ptr<fe_type> getFE(void) const {
    return fe_vol_;
  }

}; // PDE_Semilinear

#endif
