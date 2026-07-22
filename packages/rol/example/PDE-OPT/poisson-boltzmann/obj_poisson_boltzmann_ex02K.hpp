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

#ifndef PDEOPT_QOI_L2TRACKING_POISSON_BOLTZMANN_EX02K_HPP
#define PDEOPT_QOI_L2TRACKING_POISSON_BOLTZMANN_EX02K_HPP

#include "../TOOLS/qoiK.hpp"
#include "dopingK.hpp"
#include "pde_poisson_boltzmann_ex02K.hpp"

template <class Real, class DeviceType>
class QoI_State_Cost_1_Poisson_Boltzmann : public QoI<Real, DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using fe_type = FE<Real, DeviceType>;
  using fst = Intrepid2::FunctionSpaceTools<DeviceType>;
  using rst = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> fe_vol_;
  const std::vector<std::vector<ROL::Ptr<fe_type>>> fe_bdry_;
  const std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;

  scalar_view getBoundaryCoeff(
      const scalar_view cell_coeff,
      const int sideSet, const int locSideId) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][locSideId];
    const int numCellsSide = bdryCellLocId.size();
    const int f = fe_vol_->N().extent_int(1);

    scalar_view bdry_coeff("bdry_coeff", numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        bdry_coeff(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

public:
  QoI_State_Cost_1_Poisson_Boltzmann(const ROL::Ptr<fe_type> & fe_vol,
                                     const std::vector<std::vector<ROL::Ptr<fe_type>>> & fe_bdry,
                                     const std::vector<std::vector<std::vector<int>>> & bdryCellLocIds,
                                     Teuchos::ParameterList &parlist)
    : fe_vol_(fe_vol), fe_bdry_(fe_bdry), bdryCellLocIds_(bdryCellLocIds) {}

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const int c = fe_vol_->gradN().extent_int(0);
    const int d = fe_vol_->gradN().extent_int(3);
    // Initialize output val
    val = scalar_view("val", c);
    // Compute cost integral
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      if ( i == 3 ) {
        const int numLocSides = bdryCellLocIds_[i].size();
        for (int j = 0; j < numLocSides; ++j) {
          const int numCellsSide  = bdryCellLocIds_[i][j].size();
          if ( numCellsSide ) {
            const int numCubPerSide = fe_bdry_[i][j]->cubPts().extent_int(1);
            // Evaluate control on FE basis
            scalar_view u_coeff_bdry = getBoundaryCoeff(u_coeff, i, j);
            scalar_view gradU_eval("gradU_eval", numCellsSide, numCubPerSide, d);
            fe_bdry_[i][j]->evaluateGradient(gradU_eval, u_coeff_bdry);
            // Compute cost
            scalar_view costU("costU", numCellsSide, numCubPerSide);
            scalar_view weightU("weightU", numCellsSide, numCubPerSide);
            for (int k = 0; k < numCellsSide; ++k) {
              for (int l = 0; l < numCubPerSide; ++l) {
                costU(k,l) = gradU_eval(k, l, 1);
                weightU(k,l) = static_cast<Real>(1);
              }
            }
            // Integrate cell cost
            scalar_view intVal("intVal", numCellsSide);
            fe_bdry_[i][j]->computeIntegral(intVal, costU, weightU);
            // Add to integral value
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              val(cidx) += intVal(k);
            }
          }
        }
      }
    }
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const int c = fe_vol_->gradN().extent_int(0);
    const int f = fe_vol_->gradN().extent_int(1);
    // Initialize output val
    grad = scalar_view("grad", c, f);
    // Compute cost integral
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      if ( i == 3 ) {
        const int numLocSides = bdryCellLocIds_[i].size();
        for (int j = 0; j < numLocSides; ++j) {
          const int numCellsSide  = bdryCellLocIds_[i][j].size();
          if ( numCellsSide ) {
            const int numCubPerSide = fe_bdry_[i][j]->cubPts().extent_int(1);
            // Compute cost
            scalar_view weightU("weightU", numCellsSide, numCubPerSide);
            for (int k = 0; k < numCellsSide; ++k) {
              for (int l = 0; l < numCubPerSide; ++l) {
                weightU(k,l) = static_cast<Real>(1);
              }
            }
            // Compute gradient of squared L2-norm of diff
            scalar_view intGrad("intGrad", numCellsSide, f);
            fst::integrate(intGrad, weightU, (fe_bdry_[i][j]->DNDdetJ(1)), false);
            // Add to integral value
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < f; ++l) {
                grad(cidx, l) += intGrad(k,l);
              }
            }
          }
        }
      }
    }
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_ControlCost::gradient_2 is zero.");
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_ControlCost::HessVec_11 is zero.");
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_ControlCost::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_ControlCost::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_ControlCost::HessVec_22 is zero.");
  }

}; // QoI_State_Cost_2_Poisson_Boltzmann

template <class Real, class DeviceType>
class QoI_State_Cost_2_Poisson_Boltzmann : public QoI<Real, DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using fe_type = FE<Real, DeviceType>;
  using fst = Intrepid2::FunctionSpaceTools<DeviceType>;
  using rst = Intrepid2::RealSpaceTools<DeviceType>;

private:
  ROL::Ptr<fe_type> fe_;

  scalar_view target_;

  Real targetFunc(const std::vector<Real> & x) const {
    return static_cast<Real>(1);
//    int size = x.size();
//    Real val(0);
//    for (int i = 0; i < size; ++i) {
//      val += x[i]*x[i];
//    }
//    return val;
  }

public:
  QoI_State_Cost_2_Poisson_Boltzmann(const ROL::Ptr<fe_type> & fe) : fe_(fe) {
    int c = fe_->cubPts().extent_int(0);
    int p = fe_->cubPts().extent_int(1);
    int d = fe_->cubPts().extent_int(2);
    std::vector<Real> pt(d);
    target_ = scalar_view("target_", c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (fe_->cubPts())(i,j,k);
        }
        target_(i,j) = targetFunc(pt);
      }
    }
  }

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = u_coeff.extent_int(0);
    int p = fe_->cubPts().extent_int(1);
    // Initialize output val
    val = scalar_view("val", c);
    // Evaluate state on FE basis
    scalar_view valU_eval("valU_eval", c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute difference between state and target
    rst::subtract(valU_eval, target_);
    // Compute squared L2-norm of diff
    fe_->computeIntegral(val, valU_eval, valU_eval);
    // Scale by one half
    rst::scale(val, static_cast<Real>(0.5));
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
    grad = scalar_view("grad1", c, f);
    // Evaluate state on FE basis
    scalar_view valU_eval("valU_eval", c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute difference between state and target
    rst::subtract(valU_eval, target_);
    // Compute gradient of squared L2-norm of diff
    fst::integrate(grad, valU_eval, fe_->NdetJ(), false);
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_State_Cost_Poisson_Boltzmann::gradient_2 is zero.");
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    int c = v_coeff.extent_int(0);
    int p = fe_->cubPts().extent_int(1);
    int f = fe_->N().extent_int(1);
    scalar_view valV_eval("valV_eval", c, p);
    hess = scalar_view("hess11", c, f);
    fe_->evaluateValue(valV_eval, v_coeff);
    fst::integrate(hess, valV_eval, fe_->NdetJ(), false);
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_State_Cost_Poisson_Boltzmann::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_State_Cost_Poisson_Boltzmann::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_State_Cost_Poisson_Boltzmann::HessVec_22 is zero.");
  }

}; // QoI_State_Cost

template <class Real, class DeviceType>
class QoI_Control_Cost_1_Poisson_Boltzmann : public QoI<Real, DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using fe_type = FE<Real, DeviceType>;
  using fst = Intrepid2::FunctionSpaceTools<DeviceType>;
  using rst = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> fe_;
  scalar_view target_;

public:
  QoI_Control_Cost_1_Poisson_Boltzmann(const ROL::Ptr<fe_type> & fe,
                                       const ROL::Ptr<Doping<Real, DeviceType>> & dope)
    : fe_(fe) {
    // Compute reference doping profile
    const int c = fe_->gradN().extent_int(0);
    const int p = fe_->gradN().extent_int(2);
    const int d = fe_->gradN().extent_int(3);
    std::vector<Real> pt(d);
    target_ = scalar_view("target_", c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (fe_->cubPts())(i,j,k);
        }
        target_(i,j) = dope->evaluate(pt);
      }
    }
  }

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN().extent_int(0);
    const int p = fe_->gradN().extent_int(2);
    // Initialize output val
    val = scalar_view("val", c);
    // Evaluate state on FE basis
    scalar_view valZ_eval("valZ_eval", c, p);
    fe_->evaluateValue(valZ_eval, z_coeff);
    // Compute difference between control and target
    rst::subtract(valZ_eval, target_);
    // Compute squared L2-norm of diff
    fe_->computeIntegral(val, valZ_eval, valZ_eval);
    // Scale by one half
    rst::scale(val, static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_1_Poisson_Boltzmann::gradient_1 is zero.");
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);
    const int p = fe_->gradN().extent_int(2);
    // Initialize output grad
    grad = scalar_view("grad2", c, f);
    // Evaluate state on FE basis
    scalar_view valZ_eval("valZ_eval", c, p);
    fe_->evaluateValue(valZ_eval, z_coeff);
    // Compute difference between control and target
    rst::subtract(valZ_eval, target_);
    // Compute gradient of squared L2-norm of diff
    fst::integrate(grad, valZ_eval, fe_->NdetJ(), false);
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_1_Poisson_Boltzmann::HessVec_11 is zero.");
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_1_Poisson_Boltzmann::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_1_Poisson_Boltzmann::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);
    const int p = fe_->gradN().extent_int(2);
    scalar_view valV_eval("valV_eval", c, p);
    hess = scalar_view("hess22", c, f);
    fe_->evaluateValue(valV_eval, v_coeff);
    fst::integrate(hess, valV_eval, fe_->NdetJ(), false);
  }

}; // QoI_Control_Cost_1

template <class Real, class DeviceType>
class QoI_Control_Cost_2_Poisson_Boltzmann : public QoI<Real, DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using fe_type = FE<Real, DeviceType>;
  using fst = Intrepid2::FunctionSpaceTools<DeviceType>;
  using rst = Intrepid2::RealSpaceTools<DeviceType>;

private:
  ROL::Ptr<fe_type> fe_;

public:
  QoI_Control_Cost_2_Poisson_Boltzmann(const ROL::Ptr<fe_type> & fe) : fe_(fe) {}

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN().extent_int(0);
    const int p = fe_->gradN().extent_int(2);
    const int d = fe_->gradN().extent_int(3);
    // Initialize output val
    val = scalar_view("val", c);
    // Evaluate state on FE basis
    scalar_view gradZ_eval("gradZ_eval", c, p, d);
    fe_->evaluateGradient(gradZ_eval, z_coeff);
    // Compute squared L2-norm of grad
    fe_->computeIntegral(val, gradZ_eval, gradZ_eval);
    // Scale by one half
    rst::scale(val, static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_2_Poisson_Boltzmann::gradient_1 is zero.");
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);
    const int p = fe_->gradN().extent_int(2);
    const int d = fe_->gradN().extent_int(3);
    // Initialize output grad
    grad = scalar_view("grad", c, f);
    // Evaluate state on FE basis
    scalar_view gradZ_eval("gradZ_eval", c, p, d);
    fe_->evaluateGradient(gradZ_eval, z_coeff);
    // Compute gradient of squared L2-norm of gradient
    fst::integrate(grad, gradZ_eval, fe_->gradNdetJ(), false);
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_2_Poisson_Boltzmann::HessVec_11 is zero.");
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_2_Poisson_Boltzmann::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_2_Poisson_Boltzmann::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);
    const int p = fe_->gradN().extent_int(2);
    const int d = fe_->gradN().extent_int(3);
    scalar_view gradV_eval("gradV_eval", c, p, d);
    hess = scalar_view("hess22", c, f);
    fe_->evaluateGradient(gradV_eval, v_coeff);
    fst::integrate(hess, gradV_eval, fe_->gradNdetJ(), false);
  }

}; // QoI_Control_Cost_2

template <class Real>
class StdObjective_Poisson_Boltzmann : public ROL::StdObjective<Real> {
private:
  const Real J_, w1_, w2_, w3_;

public:
  StdObjective_Poisson_Boltzmann(const Real J, const Real w1, const Real w2, const Real w3)
    : J_(J), w1_(w1), w2_(w2), w3_(w3) {}

  Real value(const std::vector<Real> &x, Real &tol) {
    const Real half(0.5);
    return half*w1_*std::pow(x[0]-J_,2) + w2_*x[1] + w3_*x[2];
  }

  void gradient(std::vector<Real> &g, const std::vector<Real> &x, Real &tol) {
    g[0] = w1_*(x[0]-J_);
    g[1] = w2_;
    g[2] = w3_;
  }

  void hessVec(std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    const Real zero(0);
    hv[0] = w1_*v[0];
    hv[1] = zero;
    hv[2] = zero;
  }

}; // OBJ_SCALAR

#endif
