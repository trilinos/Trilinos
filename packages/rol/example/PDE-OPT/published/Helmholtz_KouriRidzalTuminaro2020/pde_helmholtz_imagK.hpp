// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde_helmholtz.hpp
    \brief Implements the local PDE interface for the optimal control of
           Helmholtz.
*/

#ifndef PDE_HELMHOLTZ_IMAGK_HPP
#define PDE_HELMHOLTZ_IMAGK_HPP

#include "../../TOOLS/pdeK.hpp"
#include "../../TOOLS/feK.hpp"
#include "../../TOOLS/fieldhelperK.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real, class DeviceType>
class PDE_Helmholtz_Imag : public PDE<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  // Finite element basis information
  basis_ptr basisPtr_;
  std::vector<basis_ptr> basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid2::Cubature<DeviceType,Real,Real>> cellCub_, bdryCub_;
  // Cell node information
  scalar_view volCellNodes_;
  std::vector<std::vector<scalar_view>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<fe_type> fe_;
  std::vector<ROL::Ptr<fe_type>> feBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<scalar_view>> bdryCellDofValues_;

  Real waveNumber_;
  Real dampFactor_;
  Real impFactor_;
  
  scalar_view getBoundaryCoeff(const scalar_view cell_coeff, int sideSet, int cell) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][cell];
    const int numCellsSide = bdryCellLocId.size();
    const int f = basisPtr_->getCardinality();
    
    scalar_view bdry_coeff("bdry_coeff", numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j)
        bdry_coeff(i, j) = cell_coeff(bdryCellLocId[i], j);
    }
    return bdry_coeff;
  }

public:
  PDE_Helmholtz_Imag(ROL::ParameterList &parlist)
    : waveNumber_(parlist.sublist("Problem").get("Wave Number",10.0)),
      dampFactor_(parlist.sublist("Problem").get("Damping Factor", 1.0)),
      impFactor_ (parlist.sublist("Problem").get("Impedance Factor", 1.0)) {
    // Finite element fields.
    basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType,Real,Real>>();
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();        // get the cell type from the basis
    Intrepid2::DefaultCubatureFactory cubFactory;                            // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 2);    // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create<DeviceType,Real,Real>(cellType, cubDegree); // create default cubature
    int d = cellType.getDimension();
    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2); // set cubature degree, e.g., 2
    bdryCub_ = cubFactory.create<DeviceType,Real,Real>(bdryCellType, bdryCubDegree);
  }

  void residual(scalar_view & res,
                const scalar_view u_coeff,
                const scalar_view z_coeff = scalar_view(),
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);
    // Initialize residuals.
    res = scalar_view("res",c,f);
    // Evaluate/interpolate finite element fields on cells.
    scalar_view valu_eval("valu_eval", c, p);
    fe_->evaluateValue(valu_eval, u_coeff);
    // Integrate PDE term
    fst::integrate(res,valu_eval,fe_->NdetJ(),false);
    rst::scale(res, waveNumber_*dampFactor_);
    // APPLY ROBIN CONTROLS: Sideset 0
    int sideset = 0;
    int numLocalSideIds = bdryCellLocIds_[sideset].size();
    const int numCubPerSide = bdryCub_->getNumPoints();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[sideset][j].size();
      if (numCellsSide) {
        scalar_view robinRes("robinRes", numCellsSide, f);
        // Get U coefficients on Robin boundary
        scalar_view u_coeff_bdry = getBoundaryCoeff(u_coeff, sideset, j);
        // Evaluate U on FE basis
        scalar_view valu_eval_bdry("valu_eval_bdry", numCellsSide, numCubPerSide);
        feBdry_[j]->evaluateValue(valu_eval_bdry, u_coeff_bdry);
        // Compute Robin residual
        fst::integrate(robinRes,valu_eval_bdry,feBdry_[j]->NdetJ(),false);
        // Add Robin residual to volume residual
        for (int k = 0; k < numCellsSide; ++k) {
          int cidx = bdryCellLocIds_[sideset][j][k];
          for (int l = 0; l < f; ++l)
            res(cidx,l) += waveNumber_*impFactor_*robinRes(k,l);
        }
      }
    }
  }

  void Jacobian_1(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    // Initialize Jacobian.
    jac = scalar_view("jac",c,f,f);
    // Add PDE term
    Kokkos::deep_copy(jac,fe_->massMat());
    rst::scale(jac, dampFactor_*waveNumber_);
    // APPLY ROBIN CONTROL: Sideset 0
    int sideset = 0;
    int numLocalSideIds = bdryCellLocIds_[sideset].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[sideset][j].size();
      if (numCellsSide) {
        // Add Neumann control Jacobian to volume residual
        for (int k = 0; k < numCellsSide; ++k) {
          int cidx = bdryCellLocIds_[sideset][j][k];
          for (int l = 0; l < f; ++l) { 
            for (int m = 0; m < f; ++m) { 
              jac(cidx,l,m) += impFactor_*waveNumber_*(feBdry_[j]->massMat())(k,l,m);
            }
          }
        }
      }
    }
  }


  void Jacobian_2(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Helmholtz::Jacobian_2): Jacobian is zero.");
  }

  void Hessian_11(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Helmholtz::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Helmholtz::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Helmholtz::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Helmholtz::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(scalar_view & riesz) override {
    // Retrieve dimensions.
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    // Initialize Jacobian.
    riesz = scalar_view("riesz1",c,f,f);
    Kokkos::deep_copy(riesz,fe_->stiffMat());
    rst::add(riesz, fe_->massMat());
  }

  void RieszMap_2(scalar_view & riesz) override {
    // Retrieve dimensions.
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    // Initialize Jacobian.
    riesz = scalar_view("riesz2",c,f,f);
    Kokkos::deep_copy(riesz,fe_->massMat());
  }

  std::vector<basis_ptr> getFields() override {
    return basisPtrs_;
  }

  void setCellNodes(const scalar_view &volCellNodes,
                    const std::vector<std::vector<scalar_view>> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds) override {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fe_ = ROL::makePtr<fe_type>(volCellNodes_,basisPtr_,cellCub_);
    fidx_ = fe_->getBoundaryDofs();
    // Construct boundary FE
    int sideset = 0;
    int numLocSides = bdryCellNodes[sideset].size();
    feBdry_.resize(numLocSides);
    for (int j = 0; j < numLocSides; ++j) {
      if (bdryCellNodes[sideset][j] != scalar_view()) {
        feBdry_[j] = ROL::makePtr<fe_type>(bdryCellNodes[sideset][j],basisPtr_,bdryCub_,j);
      }
    }
  }

  const ROL::Ptr<fe_type> getFE(void) const {
    return fe_;
  }

  const std::vector<ROL::Ptr<fe_type>> getBdryFE(void) const {
    return feBdry_;
  }

  const std::vector<std::vector<int>> getBdryCellLocIds(const int sideset = 0) const {
    return bdryCellLocIds_[sideset];
  }

}; // PDE_Helmholtz


#endif
