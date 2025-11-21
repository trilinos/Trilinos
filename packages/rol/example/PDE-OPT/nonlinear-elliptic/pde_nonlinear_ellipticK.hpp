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

#ifndef PDE_NONLINEAR_ELLIPTICK_HPP
#define PDE_NONLINEAR_ELLIPTICK_HPP

#include "../TOOLS/pdeK.hpp"
#include "../TOOLS/feK.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real, class DeviceType>
class PDE_Nonlinear_Elliptic : public PDE<Real,DeviceType> {
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
  ROL::Ptr<fe_type> fe_;

  void computeBetaGradU(scalar_view & BgradU, const scalar_view gradU) const {
    // GET DIMENSIONS
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);
    std::vector<Real> U(d), BU(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k)
          U[k] = gradU(i,j,k);
        beta_value(BU,U);
        for (int k = 0; k < d; ++k)
          BgradU(i,j,k) = BU[k];
      }
    }
  }

  void computeBetaJacobianGradU(scalar_view BgradU, const scalar_view gradU) const {
    // GET DIMENSIONS
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);
    std::vector<Real> U(d);
    std::vector<std::vector<Real>> BU(d, U);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k)
          U[k] = gradU(i,j,k);
        beta_jacobian(BU,U);
        for (int k = 0; k < d; ++k) {
          for (int m = 0; m < d; ++m)
            BgradU(i,j,k,m) = BU[k][m];
        }
      }
    }
  }

  void computeBetaHessianGradUGradL(scalar_view & BgradUgradL,
                                    const scalar_view gradU,
                                    const scalar_view gradL) const {
    // GET DIMENSIONS
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);
    std::vector<Real> U(d), L(d);
    std::vector<std::vector<Real>> BU(d, U);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          U[k] = gradU(i,j,k);
          L[k] = gradL(i,j,k);
        }
        beta_hessian(BU,U,L);
        for (int k = 0; k < d; ++k) {
          for (int m = 0; m < d; ++m) {
            BgradUgradL(i,j,k,m) = BU[k][m];
          }
        }
      }
    }
  }

  void computeDiffusivity(scalar_view &kappa) const {
    // GET DIMENSIONS
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k)
          pt[k] = (fe_->cubPts())(i,j,k);
        // Compute diffusivity kappa
        kappa(i,j) = kappa_value(pt);
      }
    }
  }

  void computeRHS(scalar_view &rhs) const {
    // GET DIMENSIONS
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k)
          pt[k] = (fe_->cubPts())(i,j,k);
        // Compute diffusivity kappa
        rhs(i,j) = -rhs_value(pt);
      }
    }
  }

public:
  PDE_Nonlinear_Elliptic(ROL::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("Problem").get("Basis Order",1);
    if (basisOrder == 1)
      basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType,Real,Real>>();
    else if (basisOrder == 2)
      basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<DeviceType,Real,Real>>();
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();        // get the cell type from any basis
    Intrepid2::DefaultCubatureFactory cubFactory;                            // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree",2);     // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create<DeviceType,Real,Real>(cellType, cubDegree); // create default cubature
  }

  virtual void beta_value(std::vector<Real> &BU,
                          const std::vector<Real> &U) const {
    const int size = U.size();
    Real magU(0);
    for (int i = 0; i < size; ++i)
      magU += U[i] * U[i];
    for (int i = 0; i < size; ++i)
      BU[i] = magU * U[i];
  }

  virtual void beta_jacobian(std::vector<std::vector<Real>> &BU,
                             const std::vector<Real> &U) const {
    const int size = U.size();
    Real magU(0);
    for (int i = 0; i < size; ++i)
      magU += U[i] * U[i];
    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < size; ++j)
        BU[i][j] = static_cast<Real>(2) * U[j] * U[i];
      BU[i][i] += magU;
    }
  }

  virtual void beta_hessian(std::vector<std::vector<Real>> &BU,
                            const std::vector<Real> &U,
                            const std::vector<Real> &L) const {
    const int size = U.size();
    Real UL(0);
    for (int i = 0; i < size; ++i)
      UL += static_cast<Real>(2) * U[i] * L[i];
    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < size; ++j)
        BU[i][j] = static_cast<Real>(2) * (L[j] * U[i] + L[i] * U[j]);
      BU[i][i] += UL;
    }
  }

  Real kappa_value(const std::vector<Real> &x) const {
    const std::vector<Real> param = PDE<Real,DeviceType>::getParameter();
    const int ns = param.size(), d = x.size();
    if (ns > d+1) { 
      // random diffusion coefficient from i. babuska, f. nobile, r. tempone 2010.
      // simplified model for random stratified media.
      const Real one(1), two(2), three(3), eight(8), sixteen(16), half(0.5);
      const Real lc = one/sixteen, sqrtpi = std::sqrt(M_PI);
      const Real xi = std::sqrt(sqrtpi*lc), sqrt3 = std::sqrt(three);
      Real f(0), phi(0);                                          
      Real val = one + sqrt3*param[0]*std::sqrt(sqrtpi*lc*half);
      Real arg = one + sqrt3*std::sqrt(sqrtpi*lc*half);
      for (int i = 1; i < ns-(d+1); ++i) {
        f = floor(half*static_cast<Real>(i+1));     
        phi = ((i+1)%2 ? std::sin(f*M_PI*x[0]) : std::cos(f*M_PI*x[0]));
        val += xi*std::exp(-std::pow(f*M_PI*lc,two)/eight)*phi*sqrt3*param[i];
        arg += xi*std::exp(-std::pow(f*M_PI*lc,two)/eight)*std::abs(phi)*sqrt3;
      }
      return half + two*std::exp(val)/std::exp(arg);
    }
    return static_cast<Real>(1);
  }

  Real rhs_value(const std::vector<Real> &x) const {
    const std::vector<Real> param = PDE<Real,DeviceType>::getParameter();
    const int ns = param.size(), d = x.size();
    if (ns) {
      const int nk = ns - (d+1);
      const Real ten(10), two(2), pi(M_PI);
      Real val = std::pow(ten,two*param[nk]-two);
      for (int i = 0; i < d; ++i)
        val *= std::sin(two*param[nk+i+1]*pi*x[i]);
      return val;
    }
    return static_cast<Real>(0);
  }

  void residual(scalar_view & res,
                const scalar_view u_coeff,
                const scalar_view z_coeff = scalar_view(),
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // GET DIMENSIONS
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);
    // INITIALIZE RESIDUAL
    res = scalar_view("res", c, f);
    // COMPUTE STIFFNESS TERM
    scalar_view gradU_eval("gradU_eval", c, p, d);
    fe_->evaluateGradient(gradU_eval, u_coeff);
    scalar_view kappa("kappa", c, p);
    computeDiffusivity(kappa);
    scalar_view KgradU("KgradU", c, p, d);
    fst::scalarMultiplyDataData(KgradU,kappa,gradU_eval);
    fst::integrate(res,KgradU,fe_->gradNdetJ(),false);
    // ADD MASS TERM
    scalar_view valU_eval("valU_eval", c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    fst::integrate(res,valU_eval,fe_->NdetJ(),true);
    // ADD NONLINEAR TERM
    scalar_view BgradU("BgradU", c, p, d);
    computeBetaGradU(BgradU,gradU_eval);
    fst::integrate(res,BgradU,fe_->gradNdetJ(),true);
    // ADD RHS
    scalar_view rhs("rhs", c, p);
    computeRHS(rhs);
    fst::integrate(res,rhs,fe_->NdetJ(),true);
    // ADD CONTROL TERM
    if ( z_coeff != scalar_view() ) {
      scalar_view valZ_eval("valZ_eval", c, p);
      fe_->evaluateValue(valZ_eval, z_coeff);
      rst::scale(valZ_eval,static_cast<Real>(-1));
      fst::integrate(res,valZ_eval,fe_->NdetJ(),true);
    }
  }

  void Jacobian_1(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // GET DIMENSIONS
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);
    // INITIALIZE JACOBIAN
    jac = scalar_view("jac", c, f, f);
    // COMPUTE STIFFNESS TERM
    scalar_view kappa("kappa", c, p);
    computeDiffusivity(kappa);
    scalar_view KgradN("KgradN", c, f, p, d);
    fst::scalarMultiplyDataField(KgradN,kappa,fe_->gradN());
    fst::integrate(jac,KgradN,fe_->gradNdetJ(),false);
    // ADD MASS TERM
    fst::integrate(jac,fe_->N(),fe_->NdetJ(),true);
    // ADD NONLINEAR TERM
    scalar_view gradU_eval("gradU_eval", c, p, d);
    fe_->evaluateGradient(gradU_eval, u_coeff);
    scalar_view BgradU("BgradU", c, p, d, d);
    computeBetaJacobianGradU(BgradU,gradU_eval);
    scalar_view BgradUgradN("BgradUgradN", c, f, p, d);
    fst::tensorMultiplyDataField(BgradUgradN,BgradU,fe_->gradN());
    fst::integrate(jac,BgradUgradN,fe_->gradNdetJ(),true);
  }

  void Jacobian_2(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if ( z_coeff != scalar_view() ) {
      // GET DIMENSIONS
      int c = fe_->N().extent_int(0);
      int f = fe_->N().extent_int(1);
      // INITIALIZE JACOBIAN
      jac = scalar_view("jac", c, f, f);
      // ADD CONTROL TERM
      fst::integrate(jac,fe_->N(),fe_->NdetJ(),false);
      rst::scale(jac,static_cast<Real>(-1));
    }
  }

  void Hessian_11(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // GET DIMENSIONS
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);
    // INITIALIZE JACOBIAN
    hess = scalar_view("hess", c, f, f);
    // COMPUTE NONLINEAR TERM
    scalar_view gradU_eval("gradU_eval", c, p, d);
    fe_->evaluateGradient(gradU_eval, u_coeff);
    scalar_view gradL_eval("gradL_eval", c, p, d);
    fe_->evaluateGradient(gradL_eval, l_coeff);
    scalar_view BgradUgradL("BgradUgradL", c, p, d, d);
    computeBetaHessianGradUGradL(BgradUgradL,gradU_eval,gradL_eval);
    scalar_view BgradUgradLgradN("BgradUgradLgradN", c, f, p, d);
    fst::tensorMultiplyDataField(BgradUgradLgradN,BgradUgradL,fe_->gradN());
    fst::integrate(hess,BgradUgradLgradN,fe_->gradNdetJ(),false);
  }

  void Hessian_12(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Nonlinear_Elliptic:Hessian_12: Hessian is zero.");
  }

  void Hessian_21(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Nonlinear_Elliptic:Hessian_21: Hessian is zero.");
  }

  void Hessian_22(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Nonlinear_Elliptic:Hessian_22: Hessian is zero.");
  }

  void RieszMap_1(scalar_view & riesz) override {
    // GET DIMENSIONS
    int c = fe_->N().extent_int(0);
    int f = fe_->N().extent_int(1);
    // INITIALIZE RIESZ
    riesz = scalar_view("riesz2", c, f, f);
    Kokkos::deep_copy(riesz,fe_->stiffMat());
    rst::add(riesz,fe_->massMat());
  }

  void RieszMap_2(scalar_view & riesz) override {
    // GET DIMENSIONS
    int c = fe_->N().extent_int(0);
    int f = fe_->N().extent_int(1);
    // INITIALIZE RIESZ
    riesz = scalar_view("riesz2", c, f, f);
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
  }

  const ROL::Ptr<fe_type> getFE(void) const {
    return fe_;
  }

}; // PDE_Nonlinear_Elliptic

#endif
