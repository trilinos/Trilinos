// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde.hpp
    \brief Implements the local PDE interface for the Poisson control problem.
*/

#ifndef PDE_POISSONK_HPP
#define PDE_POISSONK_HPP

#include "../TOOLS/pdeK.hpp"
#include "../TOOLS/feK.hpp"

#include "Intrepid2_HGRAD_LINE_Cn_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real, class DeviceType>
class PDE_Poisson : public PDE<Real, DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type = FE<Real,DeviceType>;
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
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<scalar_view>> bdryCellDofValues_;

  bool useStateRiesz_;
  bool useControlRiesz_;
  Real alpha_;

  Real DirichletFunc(const std::vector<Real> & coords, int sideset, int locSideId) const {
    return 0;
  }

  Real evaluateRHS(const std::vector<Real> &x) const {
    const Real pi(M_PI), eight(8);
    Real s1(1), s2(1);
    int dim = x.size();
    for (int i=0; i<dim; ++i) {
      s1 *= std::sin(eight*pi*x[i]);
      s2 *= std::sin(pi*x[i]);
    }
    Real coeff1(64), coeff2(dim);
    return s1/(alpha_*coeff1*coeff2*pi*pi) + coeff2*pi*pi*s2;
  }

  void computeRHS(scalar_view rhs) const {
    // GET DIMENSIONS
    int c = fe_vol_->gradN().extent_int(0);
    int p = fe_vol_->gradN().extent_int(2);
    int d = fe_vol_->gradN().extent_int(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = fe_vol_->cubPts()(i,j,k);
        }
        // Compute forcing term f
        rhs(i,j) = -evaluateRHS(pt);
      }
    }
  }

public:
  PDE_Poisson(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("Problem").get("Basis Order",1);
    int cubDegree  = parlist.sublist("Problem").get("Cubature Degree",4);
    int probDim    = parlist.sublist("Problem").get("Problem Dimension",2);
    if (probDim > 3 || probDim < 1) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> PDE-OPT/poisson/pde_poisson.hpp: Problem dimension is not 1, 2 or 3!");
    }
    if (basisOrder > 2 || basisOrder < 1) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> PDE-OPT/poisson/pde_poisson.hpp: Basis order is not 1 or 2!");
    }
    if (probDim == 1) {
      basisPtr_ = Teuchos::rcp(new Intrepid2::Basis_HGRAD_LINE_Cn_FEM<DeviceType, Real, Real>(basisOrder,Intrepid2::POINTTYPE_EQUISPACED));
    } else if (probDim == 2) {
      if (basisOrder == 1) {
        basisPtr_ = Teuchos::rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType, Real, Real>);
      }
      else if (basisOrder == 2) {
        basisPtr_ = Teuchos::rcp(new Intrepid2::Basis_HGRAD_QUAD_C2_FEM<DeviceType, Real, Real>);
      }
    } else if (probDim == 3) {
      if (basisOrder == 1) {
        basisPtr_ = Teuchos::rcp(new Intrepid2::Basis_HGRAD_HEX_C1_FEM<DeviceType, Real, Real>);
      }
      else if (basisOrder == 2) {
        basisPtr_ = Teuchos::rcp(new Intrepid2::Basis_HGRAD_HEX_C2_FEM<DeviceType, Real, Real>);
      }
    }
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();
    Intrepid2::DefaultCubatureFactory cubFactory;
    cellCub_ = cubFactory.create<DeviceType, Real, Real>(cellType, cubDegree);
    // Problem data.
    useStateRiesz_   = parlist.sublist("Problem").get("Use State Riesz Map", true);      // use Riesz map for state variables?
    useControlRiesz_ = parlist.sublist("Problem").get("Use Control Riesz Map", true);    // use Riesz map for control variables?
    alpha_           = parlist.sublist("Problem").get("Control penalty parameter",1e-2);
  }

  void residual(scalar_view & res,
                const scalar_view u_coeff,
                const scalar_view z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // GET DIMENSIONS
    int c = u_coeff.extent_int(0);
    int p = cellCub_->getNumPoints();
    int f = basisPtr_->getCardinality();
    int d = cellCub_->getDimension();
    // INITIALIZE RESIDUAL
    res = scalar_view("res", c, f);
    // COMPUTE STIFFNESS TERM
    scalar_view gradU_eval("gradU_eval", c, p, d);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    fst::integrate(res, gradU_eval, fe_vol_->gradNdetJ(), false);
    // COMPUTE RHS
    scalar_view rhs("rhs", c, p);
    computeRHS(rhs);
    fst::integrate(res, rhs, fe_vol_->NdetJ(), true);
    // ADD CONTROL TERM TO RESIDUAL
    if ( z_coeff.is_allocated()) {
      scalar_view valZ_eval("valZ_eval", c, p);
      fe_vol_->evaluateValue(valZ_eval, z_coeff);
      rst::scale(valZ_eval,static_cast<Real>(-1));
      fst::integrate(res, valZ_eval, fe_vol_->NdetJ(), true);
    }
    // APPLY DIRICHLET CONDITIONS
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
              //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx_[j][l];
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
    int c = u_coeff.extent_int(0);
    int f = basisPtr_->getCardinality();
    // INITILAIZE JACOBIAN
    jac = scalar_view("jac", c, f, f);
    // COMPUTE STIFFNESS TERM
    fst::integrate(jac, fe_vol_->gradN(), fe_vol_->gradNdetJ(), false);
    // APPLY DIRICHLET CONDITIONS
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
              //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
              for (int m = 0; m < f; ++m) {
                jac(cidx,fidx_[j][l],m) = static_cast<Real>(0);
              }
              jac(cidx,fidx_[j][l],fidx_[j][l]) = static_cast<Real>(1);
            }
          }
        }
      }
    }
  }

  void Jacobian_2(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    if ( z_coeff.is_allocated() ) {
      // GET DIMENSIONS
      int c = u_coeff.extent_int(0);
      int f = basisPtr_->getCardinality();
      // INITIALIZE JACOBIAN
      jac = scalar_view("jac", c, f, f);
      // ADD CONTROL TERM
      fst::integrate(jac, fe_vol_->N(), fe_vol_->NdetJ(), false);
      rst::scale(jac, static_cast<Real>(-1));
      // APPLY DIRICHLET CONDITIONS
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
                //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
                for (int m = 0; m < f; ++m) {
                  jac(cidx,fidx_[j][l],m) = static_cast<Real>(0);
                }
              }
            }
          }
        }
      }
    }
    else {
      throw Exception::Zero(">>> (PDE_Poisson::Jacobian_2): Jacobian is zero.");
    }
  }

  void Hessian_11(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Poisson::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Poisson::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Poisson::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Poisson::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(scalar_view & riesz) override {
    // Optionally disable Riesz map ...
    if (!useStateRiesz_) {
      throw Exception::NotImplemented(">>> (PDE_Poisson::RieszMap_1): Not implemented.");
    }

    // ...otherwise ...

    int c = fe_vol_->N().extent_int(0);
    int f = fe_vol_->N().extent_int(1);
    // INITILAIZE JACOBIAN
    riesz = scalar_view("riesz", c, f, f);
    Kokkos::deep_copy(riesz,fe_vol_->stiffMat());
    rst::add(riesz,fe_vol_->massMat());
  }

  void RieszMap_2(scalar_view & riesz) override {
    // Optionally disable Riesz map ...
    if (!useControlRiesz_) {
      throw Exception::NotImplemented(">>> (PDE_Poisson::RieszMap_2): Not implemented.");
    }

    // ...otherwise ...

    int c = fe_vol_->N().extent_int(0);
    int f = fe_vol_->N().extent_int(1);
    // INITILAIZE JACOBIAN
    riesz = scalar_view("riesz", c, f, f);
    Kokkos::deep_copy(riesz, fe_vol_->massMat());
  }

  std::vector<basis_ptr> getFields() override   {
    return basisPtrs_;
  }

  void setCellNodes(const scalar_view &volCellNodes,
                    const std::vector<std::vector<scalar_view>> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds) override {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fe_vol_ = ROL::makePtr<fe_type>(volCellNodes_,basisPtr_,cellCub_);
    // Set local boundary DOFs.
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

}; // PDE_Poisson

#endif
