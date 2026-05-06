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

#ifndef PDE_ADV_DIFFK_HPP
#define PDE_ADV_DIFFK_HPP

#include "../../TOOLS/pdeK.hpp"
#include "../../TOOLS/feK.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real, class DeviceType>
class PDE_adv_diff : public PDE<Real,DeviceType> {
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
  ROL::Ptr<Intrepid2::Cubature<DeviceType,Real,Real>> cellCub_;
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

  std::vector<Real> mx_, my_;
  Real sx_, sy_;

public:
  PDE_adv_diff(ROL::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("PDE Poisson").get("Basis Order",1);
    if (basisOrder == 1)
      basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType,Real,Real>>();
    else if (basisOrder == 2)
      basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<DeviceType,Real,Real>>();
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();        // get the cell type from any basis
    Intrepid2::DefaultCubatureFactory cubFactory;                            // create cubature factory
    int cubDegree = parlist.sublist("PDE Poisson").get("Cubature Degree",2); // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create<DeviceType,Real,Real>(cellType, cubDegree); // create default cubature

    int nx = parlist.sublist("Problem").get("Number Controls - X", 3);
    int ny = parlist.sublist("Problem").get("Number Controls - Y", 3);
    mx_.resize(nx*ny); my_.resize(ny*ny);
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        int index = i + j*nx;
        mx_[index] = static_cast<Real>(i+1)/static_cast<Real>(nx+1);
        my_[index] = static_cast<Real>(j+1)/static_cast<Real>(ny+1);
      }
    }
    sx_ = static_cast<Real>(1)/static_cast<Real>(6*(nx+1));
    sy_ = static_cast<Real>(1)/static_cast<Real>(6*(ny+1));
  }

  void residual(scalar_view & res,
                const scalar_view u_coeff,
                const scalar_view z_coeff = scalar_view(),
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_vol_->gradN().extent_int(0);
    int f = fe_vol_->gradN().extent_int(1);
    int p = fe_vol_->gradN().extent_int(2);
    int d = fe_vol_->gradN().extent_int(3);
    // INITIALIZE RESIDUAL
    res = scalar_view("res", c, f);
    // COMPUTE PDE COEFFICIENTS
    scalar_view kappa("kappa", c, p);
    scalar_view V("V", c, p, d);
    scalar_view rhs("rhs", c, p);
    computeCoefficients(kappa,V,rhs);
    // COMPUTE DIFFUSION TERM
    // Compute grad(U)
    scalar_view gradU_eval("gradU_eval", c, p, d);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    // Multiply kappa * grad(U)
    scalar_view kappa_gradU("kappa_gradU", c, p, d);
    fst::tensorMultiplyDataData(kappa_gradU,kappa,gradU_eval);
    // Integrate (kappa * grad(U)) . grad(N)
    fst::integrate(res,kappa_gradU,fe_vol_->gradNdetJ(),false);
    // ADD ADVECTION TERM TO RESIDUAL
    // Multiply V . grad(U)
    scalar_view V_gradU("V_gradU", c, p);
    fst::dotMultiplyDataData(V_gradU,V,gradU_eval);
    // Integrate (V . grad(U)) * N
    fst::integrate(res,V_gradU,fe_vol_->NdetJ(),true);
    // ADD RHS TO RESIDUAL
    fst::integrate(res,rhs,fe_vol_->NdetJ(),true);

    // ADD CONTROL TERM TO RESIDUAL
    int size = z_param->size();
    scalar_view ctrl("ctrl", c, p);
    for (int i = 0; i < size; ++i) {
      computeControlOperator(ctrl,(*z_param)[i],i);
      fst::integrate(res,ctrl,fe_vol_->NdetJ(),true);
    }
    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
        int cidx = bdryCellLocIds_[0][j][k];
        for (int l = 0; l < numBdryDofs; ++l) {
          res(cidx,fidx_[j][l])
            = u_coeff(cidx,fidx_[j][l]) - (bdryCellDofValues_[0][j])(k,fidx_[j][l]);
        }
      }
    }
  }

  void Jacobian_1(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_vol_->gradN().extent_int(0);
    int f = fe_vol_->gradN().extent_int(1);
    int p = fe_vol_->gradN().extent_int(2);
    int d = fe_vol_->gradN().extent_int(3);
    // INITILAIZE JACOBIAN
    jac = scalar_view("jac", c, f, f);
    // COMPUTE PDE COEFFICIENTS
    scalar_view kappa("kappa", c, p);
    scalar_view V("V", c, p, d);
    scalar_view rhs("rhs", c, p);
    computeCoefficients(kappa,V,rhs);
    // COMPUTE DIFFUSION TERM
    // Multiply kappa * grad(N)
    scalar_view kappa_gradN("kappa_gradN", c, f, p, d);
    fst::tensorMultiplyDataField(kappa_gradN,kappa,fe_vol_->gradN());
    // Integrate (kappa * grad(N)) . grad(N)
    fst::integrate(jac,kappa_gradN,fe_vol_->gradNdetJ(),false);
    // ADD ADVECTION TERM TO JACOBIAN
    // Multiply V . grad(N)
    scalar_view V_gradN("V_gradN", c, f, p);
    fst::dotMultiplyDataField(V_gradN,V,fe_vol_->gradN());
    // Integrate (V . grad(U)) * N
    fst::integrate(jac,V_gradN,fe_vol_->NdetJ(),true);
    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
        int cidx = bdryCellLocIds_[0][j][k];
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

  void Jacobian_2(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Jacobian_2): Jacobian is zero.");
  }

  void Jacobian_3(std::vector<scalar_view> & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_vol_->gradN().extent_int(0);
    int f = fe_vol_->gradN().extent_int(1);
    int p = fe_vol_->gradN().extent_int(2);
    // ADD CONTROL TERM TO RESIDUAL
    int size = z_param->size();
    scalar_view ctrl("ctrl", c, p);
    for (int i = 0; i < size; ++i) {
      jac[i] = scalar_view("jac", c, f);
      computeControlOperator(ctrl,static_cast<Real>(1),i);
      fst::integrate(jac[i],ctrl,fe_vol_->NdetJ(),false);
      // APPLY DIRICHLET CONDITIONS
      int numLocalSideIds = bdryCellLocIds_[0].size();
      for (int j = 0; j < numLocalSideIds; ++j) {
        int numCellsSide = bdryCellLocIds_[0][j].size();
        int numBdryDofs = fidx_[j].size();
        for (int k = 0; k < numCellsSide; ++k) {
          int cidx = bdryCellLocIds_[0][j][k];
          for (int l = 0; l < numBdryDofs; ++l) {
            (jac[i])(cidx,fidx_[j][l]) = static_cast<Real>(0);
          }
        }
      }
    }
  }

  void Hessian_11(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_12): Hessian is zero.");
  }

  void Hessian_13(std::vector<scalar_view> & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_13): Hessian is zero.");
  }

  void Hessian_21(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_22): Hessian is zero.");
  }

  void Hessian_23(std::vector<scalar_view> & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_23): Hessian is zero.");
  }

  void Hessian_31(std::vector<scalar_view> & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_31): Hessian is zero.");
  }

  void Hessian_32(std::vector<scalar_view> & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_32): Hessian is zero.");
  }

  void Hessian_33(std::vector<std::vector<scalar_view>> & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_33): Hessian is zero.");
  }

  void RieszMap_1(scalar_view & riesz) {
    // GET DIMENSIONS
    int c = fe_vol_->N().extent_int(0);
    int f = fe_vol_->N().extent_int(1);
    // INITIALIZE RIESZ
    riesz = scalar_view("riesz", c, f, f);
    Kokkos::deep_copy(riesz,fe_vol_->stiffMat());
    rst::add(riesz,fe_vol_->massMat());
  }

  void RieszMap_2(scalar_view & riesz) {
    // GET DIMENSIONS
    int c = fe_vol_->N().extent_int(0);
    int f = fe_vol_->N().extent_int(1);
    // INITIALIZE RIESZ
    riesz = scalar_view("riesz", c, f, f);
    Kokkos::deep_copy(riesz,fe_vol_->massMat());
  }

  std::vector<basis_ptr> getFields() {
    return basisPtrs_;
  }

  void setCellNodes(const scalar_view &volCellNodes,
                    const std::vector<std::vector<scalar_view>> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds) {
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
        bdryCellDofValues_[i][j] = scalar_view("bdryCellDofValues", c, f);
        scalar_view coords("coords", c, f, d);
        if (c > 0)
          fe_vol_->computeDofCoords(coords, bdryCellNodes_[i][j]);
        for (int k=0; k<c; ++k) {
          for (int l=0; l<f; ++l) {
            std::vector<Real> dofpoint(d);
            for (int m=0; m<d; ++m)
              dofpoint[m] = coords(k, l, m);
            (bdryCellDofValues_[i][j])(k, l) = evaluateDirichlet(dofpoint, i, j);
          }
        }
      }
    }
  }

  const ROL::Ptr<fe_type> getFE(void) const {
    return fe_vol_;
  }

private:

  Real evaluateDirichlet(const std::vector<Real> & coords, int sideset, int locSideId) const {
    return static_cast<Real>(0);
  }

  Real evaluateDiffusivity(const std::vector<Real> &x) const {
    return static_cast<Real>(2.5);
  }

  void evaluateVelocity(std::vector<Real> &adv, const std::vector<Real> &x) const {
    const Real a(2.5), b(7.5);
    adv[0] = b - a*x[0];
    adv[1] =     a*x[1];
  }

  Real evaluateRHS(const std::vector<Real> &x) const {
    const int ns = 5;             
    const Real half(0.5), two(2);
    Real source(0), arg1(0), arg2(0), mag(0), x0(0), y0(0), sx(0), sy(0);
    // Upper and lower bounds on source magintudes
    const std::vector<Real> ml = {1.5, 1.2, 1.5, 1.2, 1.1};
    const std::vector<Real> mu = {2.5, 1.8, 1.9, 2.6, 1.5};
    // Upper and lower bounds on source locations
    const std::vector<Real> xl = {0.45, 0.75, 0.40, 0.05, 0.85};
    const std::vector<Real> xu = {0.55, 0.85, 0.60, 0.35, 0.95};
    const std::vector<Real> yl = {0.25, 0.55, 0.50, 0.45, 0.45};
    const std::vector<Real> yu = {0.35, 0.65, 0.70, 0.65, 0.55};
    // Upper and lower bounds on source widths
    const std::vector<Real> sxl = {0.03, 0.02, 0.01, 0.02, 0.015};
    const std::vector<Real> sxu = {0.07, 0.04, 0.05, 0.04, 0.025};
    const std::vector<Real> syl = {0.04, 0.01, 0.02, 0.02, 0.01};
    const std::vector<Real> syu = {0.12, 0.05, 0.04, 0.04, 0.03};
    for (int i=0; i<ns; ++i) {
      mag  = half*(ml[i]+mu[i]);
      x0   = half*(xl[i]+xu[i]);
      y0   = half*(yl[i]+yu[i]);
      sx   = half*(sxl[i]+sxu[i]);
      sy   = half*(syl[i]+syu[i]);
      arg1 = std::pow((x[0]-x0)/sx, two);
      arg2 = std::pow((x[1]-y0)/sy, two);
      source += mag*std::exp(-half*(arg1+arg2));
    }
    return source;
  }

  void computeCoefficients(scalar_view &kappa,
                           scalar_view &V,
                           scalar_view &rhs) const {
    // GET DIMENSIONS
    int c = fe_vol_->gradN().extent_int(0);
    int p = fe_vol_->gradN().extent_int(2);
    int d = fe_vol_->gradN().extent_int(3);
    std::vector<Real> pt(d), adv(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k)
          pt[k] = (fe_vol_->cubPts())(i,j,k);
        // Compute diffusivity kappa
        kappa(i,j) = evaluateDiffusivity(pt);
        // Compute advection velocity field V
        evaluateVelocity(adv,pt);
        for (int k = 0; k < d; ++k)
          V(i,j,k) = adv[k];
        // Compute forcing term f
        rhs(i,j) = -evaluateRHS(pt);
      }
    }
  }

  Real evaluateControlOperator(const std::vector<Real> &x, int i) const {
    const Real half(0.5);
    return -std::exp(-half*(x[0]-mx_[i])*(x[0]-mx_[i]) / (sx_*sx_)
                     -half*(x[1]-my_[i])*(x[1]-my_[i]) / (sy_*sy_));
  }
  
  void computeControlOperator(scalar_view &ctrl, Real z, int I) const {
    // GET DIMENSIONS
    int c = fe_vol_->gradN().extent_int(0);
    int p = fe_vol_->gradN().extent_int(2);
    int d = fe_vol_->gradN().extent_int(3);
    std::vector<Real> pt(d), adv(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k)
          pt[k] = (fe_vol_->cubPts())(i,j,k);
        // Compute control operator
        ctrl(i,j) = -z*evaluateControlOperator(pt,I);
      }
    }
  }

}; // PDE_stoch_adv_diff

#endif
