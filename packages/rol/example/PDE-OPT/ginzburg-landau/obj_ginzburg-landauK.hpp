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

#ifndef PDEOPT_QOI_GINZBURGLANDAUK_HPP
#define PDEOPT_QOI_GINZBURGLANDAUK_HPP

#include "../TOOLS/qoiK.hpp"
#include "pde_ginzburg-landauK.hpp"

template <class Real, class DeviceType>
class QoI_GinzburgLandau_StateTracking : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> fe_;
  const ROL::Ptr<FieldHelper<Real,DeviceType>> fieldHelper_;

  std::vector<scalar_view> target_;

  Real epsilon0_;
  Real lambda_;

protected:
  void computeTarget(void) {
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);

    target_.clear(); target_.resize(2);
    target_[0] = scalar_view("target",c,p);
    target_[1] = scalar_view("target",c,p);
   
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (fe_->cubPts())(i,j,k);
        }
        (target_[0])(i,j) = evaluateRealTarget(x);
        (target_[1])(i,j) = evaluateImagTarget(x);
      }
    } 
  }
  
public:
  QoI_GinzburgLandau_StateTracking(const ROL::Ptr<fe_type> &fe,
                                   const ROL::Ptr<FieldHelper<Real,DeviceType>> &fieldHelper,
                                   ROL::ParameterList &parlist)
    : fe_(fe), fieldHelper_(fieldHelper) {
    lambda_   = parlist.sublist("Problem").get("Current Loading",1.0);
    epsilon0_ = parlist.sublist("Problem").get("State Scaling",1.0);
  }

  virtual Real evaluateRealTarget(const std::vector<Real> &x) const = 0;

  virtual Real evaluateImagTarget(const std::vector<Real> &x) const = 0;

  Real value(scalar_view &val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c = fe_->gradN().extent_int(0);
    const int p = fe_->gradN().extent_int(2);
    // Initialize output val
    val = scalar_view("val", c);
    // Get components of the control
    std::vector<scalar_view> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    scalar_view valU_eval("valU_eval", c, p);
    for (int i=0; i<2; ++i) {
      Kokkos::deep_copy(valU_eval, static_cast<Real>(0));
      fe_->evaluateValue(valU_eval, U[i]);
      rst::subtract(valU_eval,target_[i]);
      fe_->computeIntegral(val,valU_eval,valU_eval,true);
    }
    rst::scale(val,static_cast<Real>(0.5)*lambda_/epsilon0_);
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view &grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);
    const int p = fe_->gradN().extent_int(2);
    // Initialize output grad
    std::vector<scalar_view> G(2);
    // Get components of the control
    std::vector<scalar_view> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    scalar_view valU_eval("valU_eval", c, p);
    for (int i=0; i<2; ++i) {
      Kokkos::deep_copy(valU_eval, static_cast<Real>(0));
      fe_->evaluateValue(valU_eval, U[i]);
      rst::subtract(valU_eval,target_[i]);
      G[i] = scalar_view("grad", c, f);
      fst::integrate(G[i],valU_eval,fe_->NdetJ(),false);
      rst::scale(G[i],lambda_/epsilon0_);
    }

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(scalar_view &grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_GinzburgLandau_StateTracking::gradient_2 is zero.");
  }

  void HessVec_11(scalar_view &hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);
    const int p = fe_->gradN().extent_int(2);
    // Initialize output hessvec
    std::vector<scalar_view> H(2);
    // Get components of the control
    std::vector<scalar_view> V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Evaluate on FE basis
    scalar_view valV_eval("valV_eval", c, p);
    for (int i=0; i<2; ++i) {
      Kokkos::deep_copy(valV_eval, static_cast<Real>(0));
      fe_->evaluateValue(valV_eval, V[i]);
      H[i] = scalar_view("hess", c, f);
      fst::integrate(H[i],valV_eval,fe_->NdetJ(),false);
      rst::scale(H[i],lambda_/epsilon0_);
    }

    fieldHelper_->combineFieldCoeff(hess, H);
  }

  void HessVec_12(scalar_view &hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_TopoOpt::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view &hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_TopoOpt::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view &hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_TopoOpt::HessVec_22 is zero.");
  }

}; // QoI_TopoOpt


template <class Real, class DeviceType>
class QoI_GinzburgLandau_ControlPenalty : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> fe_;
  const std::vector<ROL::Ptr<fe_type>> feBdry_;
  const std::vector<std::vector<int>> bdryCellLocIds_;
  const ROL::Ptr<FieldHelper<Real,DeviceType>> fieldHelper_;
  Real delta0_;
  Real lambda_;

  scalar_view getBoundaryCoeff(const scalar_view cell_coeff, int locSideId) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[locSideId];
    const int numCellsSide = bdryCellLocId.size();
    const int f = fe_->N().extent_int(1);
    
    scalar_view bdry_coeff("bdry_coeff", numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        bdry_coeff(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

public:
  QoI_GinzburgLandau_ControlPenalty(const ROL::Ptr<fe_type> &fe,
                                    const std::vector<ROL::Ptr<fe_type>> &feBdry,
                                    const std::vector<std::vector<int>> &bdryCellLocIds,
                                    const ROL::Ptr<FieldHelper<Real,DeviceType>> &fieldHelper,
                                    ROL::ParameterList &parlist)
  : fe_(fe), feBdry_(feBdry), bdryCellLocIds_(bdryCellLocIds), fieldHelper_(fieldHelper) {
    delta0_ = parlist.sublist("Problem").get("Control Scaling",1.0);
    lambda_ = parlist.sublist("Problem").get("Current Loading",1.0);
  }

  Real value(scalar_view &val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    const int c = fe_->gradN().extent_int(0);
    // Initialize output val
    val = scalar_view("val", c);
    // Get components of the control
    std::vector<scalar_view> Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Compute cost integral
    const int numLocSides = bdryCellLocIds_.size();
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < numLocSides; ++j) {
        const int numCellsSide = bdryCellLocIds_[j].size();
        if ( numCellsSide ) {
          const int numCubPerSide = feBdry_[j]->cubPts().extent_int(1);
          // Evaluate control on FE basis
          scalar_view z_coeff_bdry = getBoundaryCoeff(Z[i], j);
          scalar_view valZ_eval("valZ_eval", numCellsSide, numCubPerSide);
          feBdry_[j]->evaluateValue(valZ_eval, z_coeff_bdry);
          // Integrate cell cost
          scalar_view intVal("intVal", numCellsSide);
          feBdry_[j]->computeIntegral(intVal,valZ_eval,valZ_eval,false);
          // Add to integral value
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[j][k];
            val(cidx) += static_cast<Real>(0.5)*delta0_*lambda_*intVal(k);
          }
        }
      }
    }
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view &grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_GinzburgLandau_ControlCost::gradient_1 is zero.");
  }

  void gradient_2(scalar_view &grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);
    // Initialize output grad
    std::vector<scalar_view> G(2);
    // Get components of the control
    std::vector<scalar_view> Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Compute cost integral
    const int numLocSides = bdryCellLocIds_.size();
    for (int i = 0; i < 2; ++i) {
      G[i] = scalar_view("grad", c, f);
      for (int j = 0; j < numLocSides; ++j) {
        const int numCellsSide  = bdryCellLocIds_[j].size();
        if ( numCellsSide ) {
          const int numCubPerSide = feBdry_[j]->cubPts().extent_int(1);
          // Evaluate control on FE basis
          scalar_view z_coeff_bdry = getBoundaryCoeff(Z[i], j);
          scalar_view valZ_eval("valZ_eval", numCellsSide, numCubPerSide);
          feBdry_[j]->evaluateValue(valZ_eval, z_coeff_bdry);
          // Compute gradient of squared L2-norm
          scalar_view intGrad("intGrad", numCellsSide, f);
          fst::integrate(intGrad,valZ_eval,feBdry_[j]->NdetJ(),false);
          // Add to integral value
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[j][k];
            for (int l = 0; l < f; ++l) {
              (G[i])(cidx,l) += lambda_*delta0_*intGrad(k,l);
            }
          }
        }
      }
    }

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void HessVec_11(scalar_view &hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_GinzburgLandau_ControlCost::HessVec_11 is zero.");
  }

  void HessVec_12(scalar_view &hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_GinzburgLandau_ControlCost::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view &hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_GinzburgLandau_ControlCost::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view &hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);
    // Initialize output grad
    std::vector<scalar_view> H(2);
    // Get components of the control
    std::vector<scalar_view> V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Compute cost integral
    const int numLocSides = bdryCellLocIds_.size();
    for (int i = 0; i < 2; ++i) {
      H[i] = scalar_view("hess", c, f);
      for (int j = 0; j < numLocSides; ++j) {
        const int numCellsSide  = bdryCellLocIds_[j].size();
        if ( numCellsSide ) {
          const int numCubPerSide = feBdry_[j]->cubPts().extent_int(1);
          // Evaluate control on FE basis
          scalar_view v_coeff_bdry = getBoundaryCoeff(V[i], j);
          scalar_view valV_eval("valV_eval", numCellsSide, numCubPerSide);
          feBdry_[j]->evaluateValue(valV_eval, v_coeff_bdry);
          // Compute gradient of squared L2-norm of diff
          scalar_view intHess("intHess", numCellsSide, f);
          fst::integrate(intHess,valV_eval,feBdry_[j]->NdetJ(),false);
          // Add to integral value
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[j][k];
            for (int l = 0; l < f; ++l) {
              (H[i])(cidx,l) += lambda_*delta0_*intHess(k,l);
            }
          }
        }
      }
    }

    fieldHelper_->combineFieldCoeff(hess, H);
  }

}; // QoI_GinzburgLandau_ControlPenalty

#endif
