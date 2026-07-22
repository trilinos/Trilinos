// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  obj.hpp
    \brief Provides the interface for local (cell-based) objective function computations.
*/

#ifndef PDEOPT_QOI_L2TRACKING_ALLEN_CAHNK_HPP
#define PDEOPT_QOI_L2TRACKING_ALLEN_CAHNK_HPP

#include "../TOOLS/qoiK.hpp"
#include "pde_allen_cahnK.hpp"

template <class Real, class DeviceType>
class QoI_State_Cost_Allen_Cahn : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  ROL::Ptr<fe_type> fe_;
  scalar_view target_;

  Real targetFunc(const std::vector<Real> & x) const {
    Real val(0);
    for (auto xi : x) val += xi*xi;
    return val;
  }

public:
  QoI_State_Cost_Allen_Cahn(const ROL::Ptr<fe_type> &fe) : fe_(fe) {
    int c = fe_->cubPts().extent_int(0);
    int p = fe_->cubPts().extent_int(1);
    int d = fe_->cubPts().extent_int(2);
    std::vector<Real> pt(d);
    target_ = scalar_view("target_",c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k)
          pt[k] = (fe_->cubPts())(i,j,k);
        target_(i,j) = targetFunc(pt);
      }
    }
  }

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = u_coeff.extent_int(0);
    int p = fe_->cubPts().extent_int(1);
    // Initialize output val
    val = scalar_view("val", c);
    // Evaluate state on FE basis
    scalar_view valU_eval("valU_eval", c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute difference between state and target
    rst::subtract(valU_eval,target_);
    // Compute squared L2-norm of diff
    fe_->computeIntegral(val,valU_eval,valU_eval);
    // Scale by one half
    rst::scale(val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = u_coeff.extent_int(0);
    int p = fe_->cubPts().extent_int(1);
    int f = fe_->N().extent_int(1);
    // Initialize output grad
    grad = scalar_view("grad", c, f);
    // Evaluate state on FE basis
    scalar_view valU_eval("valU_eval", c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute difference between state and target
    rst::subtract(valU_eval,target_);
    // Compute gradient of squared L2-norm of diff
    fst::integrate(grad,valU_eval,fe_->NdetJ(),false);
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_State_Cost_Allen_Cahn::gradient_2 is zero.");
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    int c = v_coeff.extent_int(0);
    int p = fe_->cubPts().extent_int(1);
    int f = fe_->N().extent_int(1);
    hess = scalar_view("hess", c, f);
    scalar_view valV_eval("valV_eval", c, p);
    fe_->evaluateValue(valV_eval, v_coeff);
    fst::integrate(hess,valV_eval,fe_->NdetJ(),false);
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_State_Cost_Allen_Cahn::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_State_Cost_Allen_Cahn::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_State_Cost_Allen_Cahn::HessVec_22 is zero.");
  }

}; // QoI_L2Tracking

template <class Real, class DeviceType>
class QoI_Control_Cost_Allen_Cahn : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> fe_vol_;
  const std::vector<std::vector<ROL::Ptr<fe_type>>> fe_bdry_;
  const std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;

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

public:
  QoI_Control_Cost_Allen_Cahn(const ROL::Ptr<fe_type> &fe_vol,
                  const std::vector<std::vector<ROL::Ptr<fe_type>>> &fe_bdry,
                  const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds)
    : fe_vol_(fe_vol), fe_bdry_(fe_bdry), bdryCellLocIds_(bdryCellLocIds) {}

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    const int c = fe_vol_->gradN().extent_int(0);
    // Initialize output val
    val = scalar_view("val", c);
    // Compute cost integral
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      const int numLocSides = bdryCellLocIds_[i].size();
      for (int j = 0; j < numLocSides; ++j) {
        const int numCellsSide  = bdryCellLocIds_[i][j].size();
        if ( numCellsSide ) {
          const int numCubPerSide = fe_bdry_[i][j]->cubPts().extent_int(1);
          // Evaluate control on FE basis
          scalar_view z_coeff_bdry = getBoundaryCoeff(z_coeff, i, j);
          scalar_view valZ_eval("valZ_eval", numCellsSide, numCubPerSide);
          fe_bdry_[i][j]->evaluateValue(valZ_eval, z_coeff_bdry);
          // Integrate cell cost
          scalar_view intVal("intVal", numCellsSide);
          fe_bdry_[i][j]->computeIntegral(intVal,valZ_eval,valZ_eval);
          // Add to integral value
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
            val(cidx) += static_cast<Real>(0.5)*intVal(k);
          }
        }
      }
    }
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Control_Cost_Allen_Cahn::gradient_1 is zero.");
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    const int c = fe_vol_->gradN().extent_int(0);
    const int f = fe_vol_->gradN().extent_int(1);
    // Initialize output val
    grad = scalar_view("grad", c, f);
    // Compute cost integral
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      const int numLocSides = bdryCellLocIds_[i].size();
      for (int j = 0; j < numLocSides; ++j) {
        const int numCellsSide  = bdryCellLocIds_[i][j].size();
        if ( numCellsSide ) {
          const int numCubPerSide = fe_bdry_[i][j]->cubPts().extent_int(1);
          // Evaluate control on FE basis
          scalar_view z_coeff_bdry = getBoundaryCoeff(z_coeff, i, j);
          scalar_view valZ_eval("valZ_eval", numCellsSide, numCubPerSide);
          fe_bdry_[i][j]->evaluateValue(valZ_eval, z_coeff_bdry);
          // Compute gradient of squared L2-norm of diff
          scalar_view intGrad("intGrad", numCellsSide, f);
          fst::integrate(intGrad,valZ_eval,fe_bdry_[i][j]->NdetJ(),false);
          // Add to integral value
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
            for (int l = 0; l < f; ++l)
              grad(cidx,l) += intGrad(k,l);
          }
        }
      }
    }
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Control_Cost_Allen_Cahn::HessVec_11 is zero.");
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Control_Cost_Allen_Cahn::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Control_Cost_Allen_Cahn::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    const int c = fe_vol_->gradN().extent_int(0);
    const int f = fe_vol_->gradN().extent_int(1);
    // Initialize output val
    hess = scalar_view("hess", c, f);
    // Compute cost integral
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      const int numLocSides = bdryCellLocIds_[i].size();
      for (int j = 0; j < numLocSides; ++j) {
        const int numCellsSide  = bdryCellLocIds_[i][j].size();
        if ( numCellsSide ) {
          const int numCubPerSide = fe_bdry_[i][j]->cubPts().extent_int(1);
          // Evaluate direction on FE basis
          scalar_view v_coeff_bdry = getBoundaryCoeff(v_coeff, i, j);
          scalar_view valV_eval("valV_eval", numCellsSide, numCubPerSide);
          fe_bdry_[i][j]->evaluateValue(valV_eval, v_coeff_bdry);
          // Compute hessian times a vector of cost
          scalar_view intHess("intHess", numCellsSide, f);
          fst::integrate(intHess,valV_eval,fe_bdry_[i][j]->NdetJ(),false);
          // Add to integral value
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
            for (int l = 0; l < f; ++l)
              hess(cidx,l) += intHess(k,l);
          }
        }
      }
    }
  }

}; // QoI_L2Penalty

#endif
