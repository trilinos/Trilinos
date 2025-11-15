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

#ifndef PDE_ALLEN_CAHNK_HPP
#define PDE_ALLEN_CAHNK_HPP

#include "../TOOLS/pdeK.hpp"
#include "../TOOLS/feK.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real, class DeviceType>
class PDE_Allen_Cahn : public PDE<Real,DeviceType> {
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
  ROL::Ptr<fe_type> fe_vol_;
  std::vector<std::vector<ROL::Ptr<fe_type>>> fe_bdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fidx_;

  Real uScale_, vScale_;
  Real robinCoeff_;

public:
  PDE_Allen_Cahn(ROL::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("Problem").get("Order of FE discretization",1);
    if (basisOrder == 1)
      basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType,Real,Real>>();
    else if (basisOrder == 2)
      basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<DeviceType,Real,Real>>();
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();                  // get the cell type from any basis
    Intrepid2::DefaultCubatureFactory cubFactory;                                      // create cubature factory
    int cubDegree = parlist.sublist("PDE Poisson Boltzmann").get("Cubature Degree",2); // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create<DeviceType,Real,Real>(cellType, cubDegree);           // create default cubature

    int d = cellType.getDimension();
    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2); // set cubature degree, e.g., 2
    bdryCub_ = cubFactory.create<DeviceType,Real,Real>(bdryCellType, bdryCubDegree);

    uScale_     = parlist.sublist("Problem").get("Scale Third Order Term",1.0);
    vScale_     = parlist.sublist("Problem").get("Scale First Order Term",1.0);
    robinCoeff_ = parlist.sublist("Problem").get("Robin Coefficient",1e0);
  }

  void residual(scalar_view & res,
                const scalar_view u_coeff,
                const scalar_view z_coeff = scalar_view(),
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN().extent_int(0);
    const int f = fe_vol_->gradN().extent_int(1);
    const int p = fe_vol_->gradN().extent_int(2);
    const int d = fe_vol_->gradN().extent_int(3);
    // INITIALIZE RESIDUAL
    res = scalar_view("res", c, f);
    // STORAGE
    scalar_view valU_eval("valU_eval", c, p);
    scalar_view gradU_eval("gradU_eval", c, p, d);
    scalar_view lambda("lambda", c, p);
    scalar_view lambda_gradU_eval("lambda_gradU_eval", c, p, d);
    scalar_view phi_valU_eval("phi_valU_eval", c, p);
    // COMPUTE PDE COEFFICIENTS
    computeCoefficients(lambda);
    // EVALUATE STATE AT QUADRATURE POINTS
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    // ADD STIFFNESS TERM TO RESIDUAL
    fst::tensorMultiplyDataData(lambda_gradU_eval,lambda,gradU_eval);
    fst::integrate(res,lambda_gradU_eval,fe_vol_->gradNdetJ(),false);
    // ADD NONLINEAR TERM TO RESIDUAL
    computeNonlinearity(phi_valU_eval, valU_eval,0);
    fst::integrate(res,phi_valU_eval,fe_vol_->NdetJ(),true);
    // APPLY ROBIN CONDITIONS
    std::vector<std::vector<scalar_view>> bdryCellDofValues;
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      int numLocalSideIds = bdryCellLocIds_[i].size();
      for (int j = 0; j < numLocalSideIds; ++j) {
        int numCellsSide = bdryCellLocIds_[i][j].size();
        const int numCubPerSide = bdryCub_->getNumPoints();
        if (numCellsSide) {
          // Get U and Z coefficients on Robin boundary
          scalar_view u_coeff_bdry = getBoundaryCoeff(u_coeff, i, j);
          scalar_view z_coeff_bdry = getBoundaryCoeff(z_coeff, i, j);
          // Evaluate U and Z on FE basis
          scalar_view valU_eval_bdry("valU_eval_bdry", numCellsSide, numCubPerSide);
          scalar_view valZ_eval_bdry("valZ_eval_bdry", numCellsSide, numCubPerSide);
          fe_bdry_[i][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
          fe_bdry_[i][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
          // Compute Robin residual
          scalar_view robinRes("robinRes", numCellsSide, f);
          scalar_view robinVal("robinVal", numCellsSide, numCubPerSide);
          computeRobin(robinVal,valU_eval_bdry,valZ_eval_bdry,i,j,0);
          fst::integrate(robinRes,robinVal,fe_bdry_[i][j]->NdetJ(),false);
          // Add Robin residual to volume residual
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
            for (int l = 0; l < f; ++l)
              res(cidx,l) += robinRes(k,l);
          }
        }
      }
    }
  }

  void Jacobian_1(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN().extent_int(0);
    const int f = fe_vol_->gradN().extent_int(1);
    const int p = fe_vol_->gradN().extent_int(2);
    const int d = fe_vol_->gradN().extent_int(3);
    // INITIALIZE JACOBIAN
    jac = scalar_view("jac", c, f, f);
    // STORAGE
    scalar_view valU_eval("valU_eval", c, p);
    scalar_view lambda("lambda", c, p);
    scalar_view lambda_gradN_eval("lambda_gradN_eval", c, f, p, d);
    scalar_view dphi_valU_eval("dphi_valU_eval", c, p);
    scalar_view NdphiU("NdphiU", c, f, p);
    // COMPUTE PDE COEFFICIENTS
    computeCoefficients(lambda);
    // EVALUATE STATE AT QUADRATURE POINTS
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    // ADD STIFFNESS TERM TO JACOBIAN
    fst::tensorMultiplyDataField(lambda_gradN_eval,lambda,fe_vol_->gradN());
    fst::integrate(jac,lambda_gradN_eval,fe_vol_->gradNdetJ(),false);
    // ADD NONLINEAR TERM
    computeNonlinearity(dphi_valU_eval, valU_eval, 1);
    fst::scalarMultiplyDataField(NdphiU,dphi_valU_eval,fe_vol_->N());
    fst::integrate(jac,NdphiU,fe_vol_->NdetJ(),true);
    // APPLY ROBIN CONDITIONS
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      int numLocalSideIds = bdryCellLocIds_[i].size();
      for (int j = 0; j < numLocalSideIds; ++j) {
        int numCellsSide = bdryCellLocIds_[i][j].size();
        const int numCubPerSide = bdryCub_->getNumPoints();
        if (numCellsSide) {
          // Get U and Z coefficients on Robin boundary
          scalar_view u_coeff_bdry = getBoundaryCoeff(u_coeff, i, j);
          scalar_view z_coeff_bdry = getBoundaryCoeff(z_coeff, i, j);
          // Evaluate U and Z on FE basis
          scalar_view valU_eval_bdry("valU_eval_bdry", numCellsSide, numCubPerSide);
          scalar_view valZ_eval_bdry("valZ_eval_bdry", numCellsSide, numCubPerSide);
          fe_bdry_[i][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
          fe_bdry_[i][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
          // Compute Robin residual
          scalar_view robinVal_N("robinVal_N", numCellsSide, f, numCubPerSide);
          scalar_view robinJac("robinJac", numCellsSide, f, f);
          scalar_view robinVal("robinVal", numCellsSide, numCubPerSide);
          computeRobin(robinVal,valU_eval_bdry,valZ_eval_bdry,i,j,1,1);
          fst::scalarMultiplyDataField(robinVal_N,robinVal,fe_bdry_[i][j]->N());
          fst::integrate(robinJac,robinVal_N,fe_bdry_[i][j]->NdetJ(),false);
          // Add Robin residual to volume residual
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
            for (int l = 0; l < f; ++l) { 
              for (int m = 0; m < f; ++m)
                jac(cidx,l,m) += robinJac(k,l,m);
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
    // GET DIMENSIONS
    const int c = fe_vol_->gradN().extent_int(0);
    const int f = fe_vol_->gradN().extent_int(1);
    // INITIALIZE JACOBIAN
    jac = scalar_view("jac", c, f, f);
    // APPLY ROBIN CONDITIONS
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      int numLocalSideIds = bdryCellLocIds_[i].size();
      for (int j = 0; j < numLocalSideIds; ++j) {
        int numCellsSide = bdryCellLocIds_[i][j].size();
        const int numCubPerSide = bdryCub_->getNumPoints();
        if (numCellsSide) {
          // Get U and Z coefficients on Robin boundary
          scalar_view u_coeff_bdry = getBoundaryCoeff(u_coeff, i, j);
          scalar_view z_coeff_bdry = getBoundaryCoeff(z_coeff, i, j);
          // Evaluate U and Z on FE basis
          scalar_view valU_eval_bdry("valU_eval_bdry", numCellsSide, numCubPerSide);
          scalar_view valZ_eval_bdry("valZ_eval_bdry", numCellsSide, numCubPerSide);
          fe_bdry_[i][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
          fe_bdry_[i][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
          // Compute Robin residual
          scalar_view robinVal_N("robinVal_N", numCellsSide, f, numCubPerSide);
          scalar_view robinJac("robinJac", numCellsSide, f, f);
          scalar_view robinVal("robinVal", numCellsSide, numCubPerSide);
          computeRobin(robinVal,valU_eval_bdry,valZ_eval_bdry,i,j,1,2);
          fst::scalarMultiplyDataField(robinVal_N,robinVal,fe_bdry_[i][j]->N());
          fst::integrate(robinJac,robinVal_N,fe_bdry_[i][j]->NdetJ(),false);
          // Add Robin residual to volume residual
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
            for (int l = 0; l < f; ++l) { 
              for (int m = 0; m < f; ++m)
                jac(cidx,l,m) += robinJac(k,l,m);
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
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN().extent_int(0);
    const int f = fe_vol_->gradN().extent_int(1);
    const int p = fe_vol_->gradN().extent_int(2);
    // INITIALIZE HESSIAN
    hess = scalar_view("hess", c, f, f);
    // STORAGE
    scalar_view valU_eval("valU_eval", c, p);
    scalar_view valL_eval("valL_eval", c, p);
    scalar_view d2phi_valU_eval("d2phi_valU_eval", c, p);
    scalar_view Ld2phiU("Ld2phiU", c, p);
    scalar_view NLd2phiU("NLd2phiU", c, f, p);
    // EVALUATE STATE AT QUADRATURE POINTS
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    fe_vol_->evaluateValue(valL_eval, l_coeff);
    // COMPUTE NONLINEAR TERM
    computeNonlinearity(d2phi_valU_eval, valU_eval, 2);
    fst::scalarMultiplyDataData(Ld2phiU,d2phi_valU_eval,valL_eval);
    fst::scalarMultiplyDataField(NLd2phiU,Ld2phiU,fe_vol_->N());
    fst::integrate(hess,NLd2phiU,fe_vol_->NdetJ(),false);
  }

  void Hessian_12(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Allen_Cahn:Hessian_12: Hessian is zero.");
  }

  void Hessian_21(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Allen_Cahn:Hessian_21: Hessian is zero.");
  }

  void Hessian_22(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Allen_Cahn:Hessian_22: Hessian is zero.");
  }

  void RieszMap_1(scalar_view & riesz) override {
    // GET DIMENSIONS
    const int c = fe_vol_->N().extent_int(0);
    const int f = fe_vol_->N().extent_int(1);
    // INITIALIZE RIESZ
    riesz = scalar_view("riesz1", c, f, f);
    Kokkos::deep_copy(riesz,fe_vol_->stiffMat());
    rst::add(riesz,fe_vol_->massMat());
  }

  void RieszMap_2(scalar_view & riesz) override {
    // GET DIMENSIONS
    const int c = fe_vol_->N().extent_int(0);
    const int f = fe_vol_->N().extent_int(1);
    // INITIALIZE RIESZ
    riesz = scalar_view("riesz2", c, f, f);
    Kokkos::deep_copy(riesz,fe_vol_->massMat());
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
    fe_vol_ = ROL::makePtr<fe_type>(volCellNodes_,basisPtr_,cellCub_);
    // Set local boundary DOFs.
    fidx_ = fe_vol_->getBoundaryDofs();
    // Construct boundary FEs
    const int numSidesets = bdryCellNodes.size();
    fe_bdry_.resize(numSidesets);
    for(int i = 0; i < numSidesets; ++i) {
      int numLocSides = bdryCellNodes[i].size();
      fe_bdry_[i].resize(numLocSides);
      for (int j = 0; j < numLocSides; ++j) {
        if (bdryCellNodes_[i][j] != scalar_view()) {
          fe_bdry_[i][j] = ROL::makePtr<fe_type>(bdryCellNodes_[i][j],basisPtr_,bdryCub_,j);
        }
      }
    }
  }

  const ROL::Ptr<fe_type> getFE(void) const {
    return fe_vol_;
  }

  const std::vector<std::vector<ROL::Ptr<fe_type>>> getBdryFE(void) const {
    return fe_bdry_;
  }

  const scalar_view getCellNodes(void) const {
    return volCellNodes_;
  }

  const std::vector<std::vector<std::vector<int>>> getBdryCellLocIds(void) const {
    return bdryCellLocIds_;
  }

private:

  Real evaluateLambda(const std::vector<Real> &x) const {
    const std::vector<Real> param = PDE<Real,DeviceType>::getParameter();
    Real val(1);
    if (param.size()) {
      val += param[0];
    }
    return val;
  }

  void computeCoefficients(scalar_view &lambda) const {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN().extent_int(0);
    const int p = fe_vol_->gradN().extent_int(2);
    const int d = fe_vol_->gradN().extent_int(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k)
          pt[k] = (fe_vol_->cubPts())(i,j,k);
        lambda(i,j) = evaluateLambda(pt);
      }
    }
  }

  Real evaluateDelta(const std::vector<Real> &x) const {
    const std::vector<Real> param = PDE<Real,DeviceType>::getParameter();
    Real val(1);
    if (param.size())
      val += param[1];
    return val;
  }

  Real evaluateNonlinearity(const std::vector<Real> &x, const Real u, const int deriv) const {
    Real val(0);
    if (deriv == 0)
      val = uScale_*std::pow(u,3) - vScale_*u;
    if (deriv == 1)
      val = uScale_*static_cast<Real>(3)*std::pow(u,2) - vScale_;
    if (deriv == 2)
      val = uScale_*static_cast<Real>(6)*u;
    return evaluateDelta(x) * val;
  }

  void computeNonlinearity(scalar_view &val,
                           const scalar_view &u,
                           const int deriv = 0) const {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN().extent_int(0);
    const int p = fe_vol_->gradN().extent_int(2);
    const int d = fe_vol_->gradN().extent_int(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k)
          pt[k] = (fe_vol_->cubPts())(i,j,k);
        val(i,j) = evaluateNonlinearity(pt, u(i,j), deriv);
      }
    }
  }

  Real evaluateRobin(const Real u, const Real z, const std::vector<Real> &x,
                     const int sideset, const int locSideId,
                     const int deriv = 0, const int component = 1) const {
    const std::vector<Real> param = PDE<Real,DeviceType>::getParameter();
    Real h(robinCoeff_);
    if (param.size())
      h += param[2];
    if ( deriv == 1 )
      return (component==1) ? h : -h;
    if ( deriv > 1 )
      return static_cast<Real>(0);
    return h * (u - z);
  }

  void computeRobin(scalar_view &robin,
                    const scalar_view &u,
                    const scalar_view &z,
                    const int sideset,
                    const int locSideId,
                    const int deriv = 0,
                    const int component = 1) const {
    const int c = u.extent_int(0);
    const int p = u.extent_int(1);
    const int d = fe_vol_->gradN().extent_int(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k)
          pt[k] = (fe_bdry_[sideset][locSideId]->cubPts())(i,j,k);
        robin(i,j) = evaluateRobin(u(i,j),z(i,j),pt,sideset,locSideId,deriv,component);
      }
    }
  }

  scalar_view getBoundaryCoeff(const scalar_view cell_coeff, const int sideset, const int locSideId) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideset][locSideId];
    const int numCellsSide = bdryCellLocId.size();
    const int f = fe_vol_->N().extent_int(1);
    
    scalar_view bdry_coeff("bdry_coeff", numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j)
        bdry_coeff(i, j) = cell_coeff(bdryCellLocId[i], j);
    }
    return bdry_coeff;
  }

}; // PDE_Allen_Cahn

#endif
