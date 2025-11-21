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

#ifndef OBJ_DARCYK_HPP
#define OBJ_DARCYK_HPP

#include "../../../../TOOLS/qoiK.hpp"
#include "pde_darcyK.hpp"
#include "permeabilityK.hpp"

template <class Real, class DeviceType>
class QoI_Velocity_Darcy : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> fePrs_, feCtrl_;
  const std::vector<ROL::Ptr<fe_type>> fePrsBdry_, feCtrlBdry_;
  const std::vector<std::vector<int>> bdryCellLocIds_;
  const ROL::Ptr<Permeability<Real,DeviceType>> perm_;
  std::vector<Real> target_;
  bool onlyAxial_;

  scalar_view getBoundaryCoeff(const scalar_view cell_coeff, int locSideId, const ROL::Ptr<fe_type> &fe) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[locSideId];
    const int numCellsSide = bdryCellLocId.size();
    const int f = fe->N().extent_int(1);
    
    scalar_view bdry_coeff("bdry_coeff", numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        bdry_coeff(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

public:
  QoI_Velocity_Darcy(ROL::ParameterList &list,
                      const ROL::Ptr<fe_type> &fePrs,
                      const ROL::Ptr<fe_type> &feCtrl,
                      const std::vector<ROL::Ptr<fe_type>> &fePrsBdry,
                      const std::vector<ROL::Ptr<fe_type>> &feCtrlBdry,
                      const std::vector<std::vector<int>> &bdryCellLocIds,
                      const ROL::Ptr<Permeability<Real,DeviceType>> &perm)
    : fePrs_(fePrs), feCtrl_(feCtrl),
      fePrsBdry_(fePrsBdry), feCtrlBdry_(feCtrlBdry),
      bdryCellLocIds_(bdryCellLocIds),
      perm_(perm) {
    target_.clear(); target_.resize(2);
    target_[0] = list.sublist("Problem").get("Target Radial Velocity",0.0);
    target_[1] = list.sublist("Problem").get("Target Axial Velocity",-15.0);
    onlyAxial_ = list.sublist("Problem").get("Only Use Axial Velocity",false);
  }

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c = fePrs_->gradN().extent_int(0);
    const int d = fePrs_->gradN().extent_int(3);
    // Initialize output val
    val = scalar_view("val", c);
    // Compute cost integral
    std::vector<Real> target = (z_param == ROL::nullPtr) ? target_ : *z_param;
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = fePrsBdry_[l]->cubPts().extent_int(1);
        scalar_view valU_eval("valU_eval", numCellsSide, numCubPerSide, d);
        scalar_view valZ_eval("valZ_eval", numCellsSide, numCubPerSide);
        scalar_view alpha("alpha", numCellsSide, numCubPerSide);
        scalar_view intVal("intVal", numCellsSide);
        scalar_view u_coeff_bdry = getBoundaryCoeff(u_coeff, l, fePrs_);
        scalar_view z_coeff_bdry = getBoundaryCoeff(z_coeff, l, feCtrl_);
        fePrsBdry_[l]->evaluateGradient(valU_eval, u_coeff_bdry);
        feCtrlBdry_[l]->evaluateValue(valZ_eval, z_coeff_bdry);
        perm_->compute(alpha, valZ_eval, fePrsBdry_[l]->cubPts(), 0);
        for (int i = 0; i < numCellsSide; ++i) {
          for (int j = 0; j < numCubPerSide; ++j) {
            for (int k = 0; k < d; ++k) {
              valU_eval(i,j,k) *= -alpha(i,j);
              if (k==0 && onlyAxial_) valU_eval(i,j,k) = static_cast<Real>(0);
              else                    valU_eval(i,j,k) -= target[k];
            }
          }
        }
        fePrsBdry_[l]->computeIntegral(intVal,valU_eval,valU_eval,true);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          val(cidx) += static_cast<Real>(0.5)*intVal(i);
        }
      }
    }
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c = fePrs_->gradN().extent_int(0);
    const int f = fePrs_->gradN().extent_int(1);
    const int d = fePrs_->gradN().extent_int(3);
    // Initialize output grad
    grad = scalar_view("grad", c, f);
    // Compute cost integral
    std::vector<Real> target = (z_param == ROL::nullPtr) ? target_ : *z_param;
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = fePrsBdry_[l]->cubPts().extent_int(1);
        scalar_view valU_eval("valU_eval", numCellsSide, numCubPerSide, d);
        scalar_view valZ_eval("valZ_eval", numCellsSide, numCubPerSide);
        scalar_view alpha("alpha", numCellsSide, numCubPerSide);
        scalar_view intGrad("intGrad", numCellsSide, f);
        scalar_view u_coeff_bdry = getBoundaryCoeff(u_coeff, l, fePrs_);
        scalar_view z_coeff_bdry = getBoundaryCoeff(z_coeff, l, feCtrl_);
        fePrsBdry_[l]->evaluateGradient(valU_eval, u_coeff_bdry);
        feCtrlBdry_[l]->evaluateValue(valZ_eval, z_coeff_bdry);
        perm_->compute(alpha, valZ_eval, fePrsBdry_[l]->cubPts(), 0);
        for (int i = 0; i < numCellsSide; ++i) {
          for (int j = 0; j < numCubPerSide; ++j) {
            for (int k = 0; k < d; ++k) {
              valU_eval(i,j,k) *= -alpha(i,j);
              if (k==0 && onlyAxial_) valU_eval(i,j,k) = static_cast<Real>(0);
              else                    valU_eval(i,j,k) -= target[k];
              valU_eval(i,j,k) *= -alpha(i,j);
            }
          }
        }
        // Compute gradient of squared L2-norm of diff
        fst::integrate(intGrad,valU_eval,fePrsBdry_[l]->gradNdetJ(),false);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          for (int j = 0; j < f; ++j) {
            grad(cidx,j) += intGrad(i,j);
          }
        }
      }
    }
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c  = fePrs_->gradN().extent_int(0);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int d  = fePrs_->gradN().extent_int(3);
    // Initialize output grad
    grad = scalar_view("grad", c, fc);
    // Compute cost integral
    std::vector<Real> target = (z_param == ROL::nullPtr) ? target_ : *z_param;
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = fePrsBdry_[l]->cubPts().extent_int(1);
        scalar_view valU_eval("valU_eval", numCellsSide, numCubPerSide, d);
        scalar_view valZ_eval("valZ_eval", numCellsSide, numCubPerSide);
        scalar_view alpha("alpha", numCellsSide, numCubPerSide);
        scalar_view dalpha("dalpha", numCellsSide, numCubPerSide);
        scalar_view integrand("integrad", numCellsSide, numCubPerSide);
        scalar_view intGrad("intGrad", numCellsSide, fc);
        scalar_view u_coeff_bdry = getBoundaryCoeff(u_coeff, l, fePrs_);
        scalar_view z_coeff_bdry = getBoundaryCoeff(z_coeff, l, feCtrl_);
        fePrsBdry_[l]->evaluateGradient(valU_eval, u_coeff_bdry);
        feCtrlBdry_[l]->evaluateValue(valZ_eval, z_coeff_bdry);
        perm_->compute(alpha, valZ_eval, fePrsBdry_[l]->cubPts(), 0);
        perm_->compute(dalpha, valZ_eval, fePrsBdry_[l]->cubPts(), 1);
        Real dalphaU(0), misfit(0);
        for (int i = 0; i < numCellsSide; ++i) {
          for (int j = 0; j < numCubPerSide; ++j) {
            for (int k = 0; k < d; ++k) {
              dalphaU = -dalpha(i,j) * valU_eval(i,j,k);
              if (k==0 && onlyAxial_) misfit = static_cast<Real>(0);
              else                    misfit = -alpha(i,j) * valU_eval(i,j,k) - target[k];
              integrand(i,j) += dalphaU * misfit;
            }
          }
        }
        // Compute gradient of squared L2-norm of diff
        fst::integrate(intGrad,integrand,feCtrlBdry_[l]->NdetJ(),false);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          for (int j = 0; j < fc; ++j) {
            grad(cidx,j) += intGrad(i,j);
          }
        }
      }
    }
  }

  std::vector<Real> gradient_3(std::vector<scalar_view> & grad,
                               const scalar_view u_coeff,
                               const scalar_view z_coeff = ROL::nullPtr,
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN().extent_int(0);
      const int d = fePrs_->gradN().extent_int(3);
      // Initialize output val
      std::vector<Real> g_param(d,static_cast<Real>(0));
      grad.clear(); grad.resize(d);
      for (int k = 0; k < d; ++k)
        grad[k] = scalar_view("grad", c);
      // Compute cost integral
      const int numLocSides = bdryCellLocIds_.size();
      for (int l = 0; l < numLocSides; ++l) {
        const int numCellsSide  = bdryCellLocIds_[l].size();
        if ( numCellsSide ) {
          const int numCubPerSide = fePrsBdry_[l]->cubPts().extent_int(1);
          std::vector<scalar_view> intVal(d);
          scalar_view valU_eval("valU_eval", numCellsSide, numCubPerSide, d);
          scalar_view valZ_eval("valZ_eval", numCellsSide, numCubPerSide);
          scalar_view integrand("integrand", numCellsSide, numCubPerSide);
          scalar_view weight("weight", numCellsSide, numCubPerSide);
          scalar_view alpha("alpha", numCellsSide, numCubPerSide);
          scalar_view u_coeff_bdry = getBoundaryCoeff(u_coeff, l, fePrs_);
          scalar_view z_coeff_bdry = getBoundaryCoeff(z_coeff, l, feCtrl_);
          fePrsBdry_[l]->evaluateGradient(valU_eval, u_coeff_bdry);
          feCtrlBdry_[l]->evaluateValue(valZ_eval, z_coeff_bdry);
          perm_->compute(alpha, valZ_eval, fePrsBdry_[l]->cubPts(), 0);
          for (int k = 0; k < d; ++k) {
            intVal[k] = scalar_view(numCellsSide);
            if ((k==0 && !onlyAxial_) || k==1) {
              integrand->initialize();
              for (int i = 0; i < numCellsSide; ++i) {
                for (int j = 0; j < numCubPerSide; ++j) {
                  weight(i,j) = static_cast<Real>(1);
                  integrand(i,j) = (*z_param)[k] + alpha(i,j) * valU_eval(i,j,k);
                }
              }
              fePrsBdry_[l]->computeIntegral(intVal[k],weight,integrand,true);
            }
          }
          // Add to integral value
          for (int k = 0; k < d; ++k) {
            for (int i = 0; i < numCellsSide; ++i) {
              int cidx = bdryCellLocIds_[l][i];
              (grad[k])(cidx) += (intVal[k])(i);
            }
          }
        }
      }
      return g_param;
    }
    else {
      throw Exception::Zero(">>> QoI_Velocity_Darcy::gradient_3 is zero.");
    }
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c = fePrs_->gradN().extent_int(0);
    const int f = fePrs_->gradN().extent_int(1);
    const int d = fePrs_->gradN().extent_int(3);
    // Initialize output grad
    hess = scalar_view("hess", c, f);
    // Compute cost integral
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = fePrsBdry_[l]->cubPts().extent_int(1);
        scalar_view valV_eval("valV_eval", numCellsSide, numCubPerSide, d);
        scalar_view valZ_eval("valZ_eval", numCellsSide, numCubPerSide);
        scalar_view alpha("alpha", numCellsSide, numCubPerSide);
        scalar_view intHess("intHess", numCellsSide, f);
        scalar_view v_coeff_bdry = getBoundaryCoeff(v_coeff, l, fePrs_);
        scalar_view z_coeff_bdry = getBoundaryCoeff(z_coeff, l, feCtrl_);
        fePrsBdry_[l]->evaluateGradient(valV_eval, v_coeff_bdry);
        feCtrlBdry_[l]->evaluateValue(valZ_eval, z_coeff_bdry);
        perm_->compute(alpha, valZ_eval, fePrsBdry_[l]->cubPts(), 0);
        for (int i = 0; i < numCellsSide; ++i) {
          for (int j = 0; j < numCubPerSide; ++j) {
            for (int k = 0; k < d; ++k) {
              valV_eval(i,j,k) *= (-alpha(i,j)) * (-alpha(i,j));
                if (k==0 && onlyAxial_) valV_eval(i,j,k) = static_cast<Real>(0);
            }
          }
        }
        // Compute hessian of squared L2-norm of diff
        fst::integrate(intHess,valV_eval,fePrsBdry_[l]->gradNdetJ(),false);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          for (int j = 0; j < f; ++j) {
            hess(cidx,j) += intHess(i,j);
          }
        }
      }
    }
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c  = fePrs_->gradN().extent_int(0);
    const int f  = fePrs_->gradN().extent_int(1);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int d  = fePrs_->gradN().extent_int(3);
    // Initialize output grad
    hess = scalar_view("hess", c, f);
    // Compute cost integral
    std::vector<Real> target = (z_param == ROL::nullPtr) ? target_ : *z_param;
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = fePrsBdry_[l]->cubPts().extent_int(1);
        scalar_view valU_eval("valU_eval", numCellsSide, numCubPerSide, d);
        scalar_view valZ_eval("valZ_eval", numCellsSide, numCubPerSide);
        scalar_view valV_eval("valV_eval", numCellsSide, numCubPerSide);
        scalar_view alpha("alpha", numCellsSide, numCubPerSide);
        scalar_view dalpha("dalpha", numCellsSide, numCubPerSide);
        scalar_view integrand("integrand", numCellsSide, numCubPerSide, d);
        scalar_view intHess("intHess", numCellsSide, f);
        scalar_view u_coeff_bdry = getBoundaryCoeff(u_coeff, l, fePrs_);
        scalar_view z_coeff_bdry = getBoundaryCoeff(z_coeff, l, feCtrl_);
        scalar_view v_coeff_bdry = getBoundaryCoeff(v_coeff, l, feCtrl_);
        fePrsBdry_[l]->evaluateGradient(valU_eval, u_coeff_bdry);
        feCtrlBdry_[l]->evaluateValue(valZ_eval, z_coeff_bdry);
        feCtrlBdry_[l]->evaluateValue(valV_eval, v_coeff_bdry);
        perm_->compute(alpha, valZ_eval, fePrsBdry_[l]->cubPts(), 0);
        perm_->compute(dalpha, valZ_eval, fePrsBdry_[l]->cubPts(), 1);
        Real dalphaV(0), misfit(0), dmisfit(0);
        for (int i = 0; i < numCellsSide; ++i) {
          for (int j = 0; j < numCubPerSide; ++j) {
            for (int k = 0; k < d; ++k) {
              dalphaV = -dalpha(i,j) * valV_eval(i,j);
              if (k==0 && onlyAxial_) {
                misfit  = static_cast<Real>(0);
                dmisfit = static_cast<Real>(0);
              }
              else {
                misfit  = -alpha(i,j)*valU_eval(i,j,k) - target[k];
                dmisfit = -alpha(i,j)*valU_eval(i,j,k);
              }
              integrand(i,j,k) = dalphaV * (misfit + dmisfit);
            }
          }
        }
        // Compute gradient of squared L2-norm of diff
        fst::integrate(intHess,integrand,fePrsBdry_[l]->gradNdetJ(),false);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          for (int j = 0; j < f; ++j) {
            hess(cidx,j) += intHess(i,j);
          }
        }
      }
    }
  }

  void HessVec_13(scalar_view & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN().extent_int(0);
      const int f = fePrs_->gradN().extent_int(1);
      const int d = fePrs_->gradN().extent_int(3);
      // Initialize output grad
      hess = scalar_view("hess", c, f);
      // Compute cost integral
      const int numLocSides = bdryCellLocIds_.size();
      for (int l = 0; l < numLocSides; ++l) {
        const int numCellsSide  = bdryCellLocIds_[l].size();
        if ( numCellsSide ) {
          const int numCubPerSide = fePrsBdry_[l]->cubPts().extent_int(1);
          scalar_view weight("weight", numCellsSide, numCubPerSide, d);
          scalar_view valZ_eval("valZ_eval", numCellsSide, numCubPerSide);
          scalar_view alpha("alpha", numCellsSide, numCubPerSide);
          scalar_view intHess("intHess", numCellsSide, f);
          scalar_view z_coeff_bdry = getBoundaryCoeff(z_coeff, l, feCtrl_);
          feCtrlBdry_[l]->evaluateValue(valZ_eval, z_coeff_bdry);
          perm_->compute(alpha, valZ_eval, fePrsBdry_[l]->cubPts(), 0);
          for (int k = 0; k < d; ++k) {
            if ((k==0 && !onlyAxial_) || k==1) {
              for (int i = 0; i < numCellsSide; ++i) {
                for (int j = 0; j < numCubPerSide; ++j) {
                  weight(i,j,k) = (*v_param)[k] * (-alpha(i,j));
                }
              }
              // Compute gradient of squared L2-norm of diff
              fst::integrate(intHess,weight,fePrsBdry_[l]->gradNdetJ(),false);
              // Add to integral value
              for (int i = 0; i < numCellsSide; ++i) {
                int cidx = bdryCellLocIds_[l][i];
                for (int j = 0; j < f; ++j) {
                  hess(cidx,j) += intHess(i,j);
                }
              }
            }
          }
        }
      }
    }
    else {
      throw Exception::NotImplemented(">>> HessVec_13 not implemented.");
    }
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c  = fePrs_->gradN().extent_int(0);
    const int f  = fePrs_->gradN().extent_int(1);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int d  = fePrs_->gradN().extent_int(3);
    // Initialize output grad
    hess = scalar_view("hess", c, fc);
    // Compute cost integral
    std::vector<Real> target = (z_param == ROL::nullPtr) ? target_ : *z_param;
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = fePrsBdry_[l]->cubPts().extent_int(1);
        scalar_view valU_eval("valU_eval", numCellsSide, numCubPerSide, d);
        scalar_view valZ_eval("valZ_eval", numCellsSide, numCubPerSide);
        scalar_view valV_eval("valV_eval", numCellsSide, numCubPerSide, d);
        scalar_view alpha("alpha", numCellsSide, numCubPerSide);
        scalar_view dalpha("dalpha", numCellsSide, numCubPerSide);
        scalar_view integrand("integrand", numCellsSide, numCubPerSide);
        scalar_view intHess("intHess", numCellsSide, fc);
        scalar_view u_coeff_bdry = getBoundaryCoeff(u_coeff, l, fePrs_);
        scalar_view z_coeff_bdry = getBoundaryCoeff(z_coeff, l, feCtrl_);
        scalar_view v_coeff_bdry = getBoundaryCoeff(v_coeff, l, fePrs_);
        fePrsBdry_[l]->evaluateGradient(valU_eval, u_coeff_bdry);
        feCtrlBdry_[l]->evaluateValue(valZ_eval, z_coeff_bdry);
        fePrsBdry_[l]->evaluateGradient(valV_eval, v_coeff_bdry);
        perm_->compute(alpha, valZ_eval, fePrsBdry_[l]->cubPts(), 0);
        perm_->compute(dalpha, valZ_eval, fePrsBdry_[l]->cubPts(), 1);
        Real dalphaV(0), misfit(0), dmisfit(0);
        for (int i = 0; i < numCellsSide; ++i) {
          for (int j = 0; j < numCubPerSide; ++j) {
            for (int k = 0; k < d; ++k) {
              dalphaV = -dalpha(i,j) * valV_eval(i,j,k);
              if (k==0 && onlyAxial_) {
                misfit  = static_cast<Real>(0);
                dmisfit = static_cast<Real>(0);
              }
              else {
                misfit  = -alpha(i,j)*valU_eval(i,j,k) - target[k];
                dmisfit = -alpha(i,j)*valU_eval(i,j,k);
              }
              integrand(i,j) += dalphaV * (misfit + dmisfit);
            }
          }
        }
        // Compute gradient of squared L2-norm of diff
        fst::integrate(intHess,integrand,feCtrlBdry_[l]->NdetJ(),false);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          for (int j = 0; j < fc; ++j) {
            hess(cidx,j) += intHess(i,j);
          }
        }
      }
    }
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c  = fePrs_->gradN().extent_int(0);
    const int f  = fePrs_->gradN().extent_int(1);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int d  = fePrs_->gradN().extent_int(3);
    // Initialize output grad
    hess = scalar_view("hess", c, fc);
    // Compute cost integral
    std::vector<Real> target = (z_param == ROL::nullPtr) ? target_ : *z_param;
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = fePrsBdry_[l]->cubPts().extent_int(1);
        scalar_view valU_eval("valU_eval", numCellsSide, numCubPerSide, d);
        scalar_view valZ_eval("valZ_eval", numCellsSide, numCubPerSide);
        scalar_view valV_eval("valV_eval", numCellsSide, numCubPerSide);
        scalar_view alpha("alpha", numCellsSide, numCubPerSide);
        scalar_view dalpha("dalpha", numCellsSide, numCubPerSide);
        scalar_view ddalpha("ddalpha", numCellsSide, numCubPerSide);
        scalar_view integrand("integrand", numCellsSide, numCubPerSide);
        scalar_view intHess("intHess", numCellsSide, fc);
        scalar_view u_coeff_bdry = getBoundaryCoeff(u_coeff, l, fePrs_);
        scalar_view z_coeff_bdry = getBoundaryCoeff(z_coeff, l, feCtrl_);
        scalar_view v_coeff_bdry = getBoundaryCoeff(v_coeff, l, feCtrl_);
        fePrsBdry_[l]->evaluateGradient(valU_eval, u_coeff_bdry);
        feCtrlBdry_[l]->evaluateValue(valZ_eval, z_coeff_bdry);
        feCtrlBdry_[l]->evaluateValue(valV_eval, v_coeff_bdry);
        perm_->compute(alpha, valZ_eval, fePrsBdry_[l]->cubPts(), 0);
        perm_->compute(dalpha, valZ_eval, fePrsBdry_[l]->cubPts(), 1);
        perm_->compute(ddalpha, valZ_eval, fePrsBdry_[l]->cubPts(), 2);
        Real dalphaV(0), misfit(0), dmisfit(0), ddalphaV(0);
        for (int i = 0; i < numCellsSide; ++i) {
          for (int j = 0; j < numCubPerSide; ++j) {
            for (int k = 0; k < d; ++k) {
              dalphaV  = -dalpha(i,j) * valV_eval(i,j) * valU_eval(i,j,k);
              ddalphaV = -ddalpha(i,j) * valV_eval(i,j) * valU_eval(i,j,k);
              if (k==0 && onlyAxial_) {
                misfit  = static_cast<Real>(0);
                dmisfit = static_cast<Real>(0);
              }
              else {
                misfit  = -alpha(i,j)*valU_eval(i,j,k) - target[k];
                dmisfit = -dalpha(i,j)*valU_eval(i,j,k);
              }
              integrand(i,j) += dalphaV * dmisfit + ddalphaV * misfit;
            }
          }
        }
        // Compute gradient of squared L2-norm of diff
        fst::integrate(intHess,integrand,feCtrlBdry_[l]->NdetJ(),false);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          for (int j = 0; j < fc; ++j) {
            hess(cidx,j) += intHess(i,j);
          }
        }
      }
    }
  }

  void HessVec_23(scalar_view & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Velocity_Darcy::HessVec_23 is zero.");
  }

  std::vector<Real> HessVec_31(std::vector<scalar_view> & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN().extent_int(0);
      const int d = fePrs_->gradN().extent_int(3);
      // Initialize output val
      std::vector<Real> h_param(d,static_cast<Real>(0));
      hess.clear(); hess.resize(d);
      for (int k = 0; k < d; ++k)
        hess[k] = scalar_view("hess", c);
      // Compute cost integral
      const int numLocSides = bdryCellLocIds_.size();
      for (int l = 0; l < numLocSides; ++l) {
        const int numCellsSide  = bdryCellLocIds_[l].size();
        if ( numCellsSide ) {
          const int numCubPerSide = fePrsBdry_[l]->cubPts().extent_int(1);
          std::vector<scalar_view> intVal(d);
          scalar_view valV_eval("valV_eval", numCellsSide, numCubPerSide);
          scalar_view weight("weight", numCellsSide, numCubPerSide);
          for (int k = 0; k < d; ++k) {
            intVal[k] = scalar_view("intVal", numCellsSide);
            if ((k==0 && !onlyAxial_) || k==1) {
              scalar_view v_coeff_bdry = getBoundaryCoeff(v_coeff, l, fePrs_);
              valV_eval->initialize();
              fePrsBdry_[l]->evaluateValue(valV_eval, v_coeff_bdry);
              for (int i = 0; i < numCellsSide; ++i) {
                for (int j = 0; j < numCubPerSide; ++j) {
                  weight(i,j) = static_cast<Real>(1);
                  valV_eval(i,j) *= static_cast<Real>(-1);
                }
              }
              fePrsBdry_[l]->computeIntegral(intVal[k],weight,valV_eval,true);
            }
          }
          // Add to integral value
          for (int k = 0; k < d; ++k) {
            for (int i = 0; i < numCellsSide; ++i) {
              int cidx = bdryCellLocIds_[l][i];
              (hess[k])(cidx) += (intVal[k])(i);
            }
          }
        }
      }
      return h_param;
    }
    else {
      throw Exception::Zero(">>> QoI_Velocity_Darcy::HessVec_31 is zero.");
    }
  }

  std::vector<Real> HessVec_32(std::vector<scalar_view> & hess,
                               const scalar_view v_coeff,
                               const scalar_view u_coeff,
                               const scalar_view z_coeff = ROL::nullPtr,
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> HessVec_32 is zero.");
  }

  std::vector<Real> HessVec_33(std::vector<scalar_view> & hess,
                               const ROL::Ptr<const std::vector<Real>> & v_param,
                               const scalar_view u_coeff,
                               const scalar_view z_coeff = ROL::nullPtr,
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN().extent_int(0);
      const int d = fePrs_->gradN().extent_int(3);
      // Initialize output val
      std::vector<Real> h_param(d,static_cast<Real>(0));
      hess.clear(); hess.resize(d);
      for (int k = 0; k < d; ++k)
        hess[k] = scalar_view("hess", c);
      // Compute cost integral
      const int numLocSides = bdryCellLocIds_.size();
      for (int l = 0; l < numLocSides; ++l) {
        const int numCellsSide  = bdryCellLocIds_[l].size();
        if ( numCellsSide ) {
          const int numCubPerSide = fePrsBdry_[l]->cubPts().extent_int(1);
          std::vector<scalar_view> intVal(d);
          scalar_view valU_eval("valU_eval", numCellsSide, numCubPerSide);
          scalar_view weight("weight", numCellsSide, numCubPerSide);
          for (int k = 0; k < d; ++k) {
            intVal[k] = scalar_view("intVal", numCellsSide);
            if ((k==0 && !onlyAxial_) || k==1) {
              valU_eval->initialize();
              for (int i = 0; i < numCellsSide; ++i) {
                for (int j = 0; j < numCubPerSide; ++j) {
                  weight(i,j) = static_cast<Real>(1);
                  valU_eval(i,j) = (*v_param)[k];
                }
              }
              fePrsBdry_[l]->computeIntegral(intVal[k],weight,valU_eval,true);
            }
          }
          // Add to integral value
          for (int k = 0; k < d; ++k) {
            for (int i = 0; i < numCellsSide; ++i) {
              int cidx = bdryCellLocIds_[l][i];
              (hess[k])(cidx) += (intVal[k])(i);
            }
          }
        }
      }
      return h_param;
    }
    else {
      throw Exception::Zero(">>> QoI_Velocity_Darcy::HessVec_33 is zero.");
    }
  }

}; // QoI_Velocity_Darcy


template <class Real, class DeviceType>
class QoI_VelocityTracking_Darcy : public QoI<Real, DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> fePrs_, feCtrl_;
  const ROL::Ptr<Permeability<Real,DeviceType>> perm_;
  scalar_view target_, weight_;
  Real rad_, yvel_, frac_, twpow_, top_;
  std::vector<Real> rvec_, zvec_;
  int targetType_;
  bool onlyAxial_, optimizeFrac_, invertFrac_, polyWeight_, veloWeight_;

  Real xTarget(const std::vector<Real> &x) const {
    const Real X = x[0], Y = x[1];
    //return xWeight(x) ? -X*Y/(rad_*rad_-Y*Y) : zero;
    //return xWeight(x) ? -X*Y/std::sqrt(rad_*rad_-Y*Y) : zero;
    //return polyWeight(x) * (-X*Y/std::sqrt(rad_*rad_-Y*Y));
    //return -X*Y/std::sqrt(rad_*rad_-Y*Y);
    if (targetType_ == 1)
      return -X*Y/((rad_*rad_-Y*Y)*(rad_*rad_-Y*Y));
    else if (targetType_ == 2) {
      Real slope = 0;
      Real xVal = 0;
      for (int i=0; i<5; ++i) {
        if ((Y >= zvec_[i]) && (Y < zvec_[i+1])) {
          slope = (rvec_[i+1] - rvec_[i]) / (zvec_[i+1] - zvec_[i]);
          xVal  = X*slope/std::pow(rvec_[i] + slope*(Y-zvec_[i]), 3);
        }
      }
      return xVal;
    }
    else
      throw Exception::NotImplemented(">>> Desired target type (not 1 or 2) not implemented.");
  }

  Real yTarget(const std::vector<Real> &x) const {
    const Real one(1), Y = x[1];
    //return yWeight(x) ? one : zero;
    //return yWeight(x) ? std::sqrt(rad_*rad_-Y*Y) : zero;
    //return polyWeight(x) * std::sqrt(rad_*rad_-Y*Y);
    //return std::sqrt(rad_*rad_-Y*Y);
    if (targetType_ == 1)
      return one/(rad_*rad_-Y*Y);
    else if (targetType_ == 2) {
      Real slope = 0;
      Real yVal = 0;
      for (int i=0; i<5; ++i) {
        if ((Y >= zvec_[i]) && (Y < zvec_[i+1])) {
          slope = (rvec_[i+1] - rvec_[i]) / (zvec_[i+1] - zvec_[i]);
          yVal  = 1.0/std::pow(rvec_[i] + slope*(Y-zvec_[i]), 2);
        }
      }
      return yVal;
    }
    else
      throw Exception::NotImplemented(">>> Desired target type (not 1 or 2) not implemented.");
  }

  Real xWeight(const std::vector<Real> &x) const {
    return yWeight(x);
  }

  Real yWeight(const std::vector<Real> &x) const {
    if (optimizeFrac_) {
      const Real zero(0), one(1), Y = x[1];
      if (invertFrac_)
        return (std::abs(Y)  > frac_*top_ ? one : zero);
      else
        return (std::abs(Y) <= frac_*top_ ? one : zero);
    }
    else if (polyWeight_)
      return polyWeight(x);
    else if (veloWeight_)
      return 1.0/(std::pow(xTarget(x), 2) + std::pow(yTarget(x), 2));
    else
      return 1.0;
  }

  Real polyWeight(const std::vector<Real> &x) const {
    const Real zero(0), one(1), Y = x[1], p = twpow_;
    Real yTOP(0), yBOT(0);
    yTOP = top_;
    yBOT = -yTOP;
    Real val = 0, at = 0, bt = 0;
    at = one / std::pow(-yTOP,p);
    bt = one / std::pow(-yBOT,p);
    if (Y > zero) {
      val = at*std::pow(Y-yTOP,p);
    } else {
      val = bt*std::pow(Y-yBOT,p);
    }
    //std::cout << Y << "  " << val << std::endl;
    return val;
  }

public:
  QoI_VelocityTracking_Darcy(ROL::ParameterList &list,
                             const ROL::Ptr<fe_type> &fePrs,
                             const ROL::Ptr<fe_type> &feCtrl,
                             const ROL::Ptr<Permeability<Real,DeviceType>> &perm)
    : fePrs_(fePrs), feCtrl_(feCtrl), perm_(perm) {

    rvec_.resize(6);
    zvec_.resize(6);

    rad_          = list.sublist("Problem").get("Diffuser Radius",10.0);
    top_          = list.sublist("Problem").get("Diffuser Top",9.9763392);
    rvec_[0]      = list.sublist("Problem").get("r0",0.6875);
    zvec_[0]      = list.sublist("Problem").get("z0",-14.001);
    rvec_[1]      = list.sublist("Problem").get("r1",3.0);
    zvec_[1]      = list.sublist("Problem").get("z1",-9.0);
    rvec_[2]      = list.sublist("Problem").get("r2",10.0);
    zvec_[2]      = list.sublist("Problem").get("z2",-3.0);
    rvec_[3]      = list.sublist("Problem").get("r3",10.0);
    zvec_[3]      = list.sublist("Problem").get("z3",3.0);
    rvec_[4]      = list.sublist("Problem").get("r4",3.0);
    zvec_[4]      = list.sublist("Problem").get("z4",9.0);
    rvec_[5]      = list.sublist("Problem").get("r5",0.6875);
    zvec_[5]      = list.sublist("Problem").get("z5",14.001);
    yvel_         = list.sublist("Problem").get("Target Axial Velocity",15.0);
    optimizeFrac_ = list.sublist("Problem").get("Optimize Domain Fraction",false);
    invertFrac_   = list.sublist("Problem").get("Invert Domain Fraction",false);
    frac_         = list.sublist("Problem").get("Integration Domain Fraction",1.00);
    polyWeight_   = list.sublist("Problem").get("Use Polynomial Weight",false);
    veloWeight_   = list.sublist("Problem").get("Use Target Velocity Weight",false);
    onlyAxial_    = list.sublist("Problem").get("Only Use Axial Velocity",false);
    targetType_   = list.sublist("Problem").get("Target Type", 1);
    twpow_        = list.sublist("Problem").get("Target Weighting Power",0.0);
    Real xWScal   = list.sublist("Problem").get("Radial Tracking Scale",1.0);
    Real yWScal   = list.sublist("Problem").get("Axial Tracking Scale",1.0);
    bool useNorm  = list.sublist("Problem").get("Use Normalized Misfit",false);
    useNorm = onlyAxial_ ? false : useNorm;
    xWScal  = onlyAxial_ ? static_cast<Real>(0) : xWScal;
    const int c = fePrs_->gradN().extent_int(0);
    const int p = fePrs_->gradN().extent_int(2);
    target_ = scalar_view("target",c,p,2);
    weight_ = scalar_view("weight",c,p,2);
    std::vector<Real> x(2);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        x[0] = (fePrs_->cubPts())(i,j,0);
        x[1] = (fePrs_->cubPts())(i,j,1);
        target_(i,j,0) = xTarget(x);
        target_(i,j,1) = yTarget(x);
        if (useNorm && yWeight(x)) {
          xWScal = static_cast<Real>(1)
                  /(std::pow(target_(i,j,0),2) + std::pow(target_(i,j,1),2));
          yWScal = xWScal;
        }
        weight_(i,j,0) = x[0] * xWScal * xWeight(x);
        weight_(i,j,1) = x[0] * yWScal * yWeight(x);
      }
    }
  }

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c = fePrs_->gradN().extent_int(0);
    const int p = fePrs_->gradN().extent_int(2);
    const int d = fePrs_->gradN().extent_int(3);
    Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
    // Initialize output val
    val = scalar_view("val",c);
    // Compute cost integral
    scalar_view gradU("gradU",c,p,d);
    scalar_view vel("vel",c,p,d);
    scalar_view wvel("wvel",c,p,d);
    scalar_view valZ("valZ",c,p);
    scalar_view alpha("alpha",c,p);
    fePrs_->evaluateGradient(gradU, u_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);

    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          vel(i,j,k)  = alpha(i,j)*gradU(i,j,k)+yvel*target_(i,j,k);
          wvel(i,j,k) = weight_(i,j,k)*vel(i,j,k);
        }
      }
    }

    fePrs_->computeIntegral(val,vel,wvel);
    // Scale by one half
    rst::scale(val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c = fePrs_->gradN().extent_int(0);
    const int f = fePrs_->gradN().extent_int(1);
    const int p = fePrs_->gradN().extent_int(2);
    const int d = fePrs_->gradN().extent_int(3);
    Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
    // Initialize output grad
    grad = scalar_view("grad", c, f);
    // Compute cost integral
    scalar_view gradU("gradU",c,p,d);
    scalar_view awvel("awvel",c,p,d);
    scalar_view valZ("valZ",c,p);
    scalar_view alpha("alpha",c,p);
    fePrs_->evaluateGradient(gradU, u_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);

    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          awvel(i,j,k) = weight_(i,j,k)*alpha(i,j)*(alpha(i,j)*gradU(i,j,k)+yvel*target_(i,j,k));
        }
      }
    }

    fst::integrate(grad,awvel,fePrs_->gradNdetJ(),false);
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c  = fePrs_->gradN().extent_int(0);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int p  = fePrs_->gradN().extent_int(2);
    const int d  = fePrs_->gradN().extent_int(3);
    Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
    // Initialize output grad
    grad = scalar_view("grad", c, fc);
    // Compute cost integral
    scalar_view gradU("gradU",c,p,d);
    scalar_view valZ("valZ",c,p);
    scalar_view alpha("alpha",c,p);
    scalar_view dalpha("dalpha",c,p);
    scalar_view deriv("deriv",c,p);
    fePrs_->evaluateGradient(gradU, u_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);
    perm_->compute(dalpha, valZ, fePrs_->cubPts(), 1);

    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          deriv(i,j) += weight_(i,j,k)*(alpha(i,j)*gradU(i,j,k)+yvel*target_(i,j,k))
                        *dalpha(i,j)*gradU(i,j,k);
        }
      }
    }

    fst::integrate(grad,deriv,feCtrl_->NdetJ(),false);
  }

  std::vector<Real> gradient_3(std::vector<scalar_view> & grad,
                               const scalar_view u_coeff,
                               const scalar_view z_coeff = ROL::nullPtr,
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN().extent_int(0);
      const int p = fePrs_->gradN().extent_int(2);
      const int d = fePrs_->gradN().extent_int(3);
      Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
      // Initialize output val
      std::vector<Real> g_param(d,static_cast<Real>(0));
      grad.clear(); grad.resize(d);
      for (int i = 0; i < d; ++i)
        grad[i] = scalar_view("grad", c);
      // Compute cost integral
      scalar_view gradU("gradU",c,p,d);
      scalar_view wvel("wvel",c,p,d);
      scalar_view valZ("valZ",c,p);
      scalar_view alpha("alpha",c,p);
      fePrs_->evaluateGradient(gradU, u_coeff);
      feCtrl_->evaluateValue(valZ, z_coeff);
      perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);

      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          for (int k = 0; k < d; ++k) {
            wvel(i,j,k) = weight_(i,j,k)*(alpha(i,j)*gradU(i,j,k)+yvel*target_(i,j,k));
          }
        }
      }

      fePrs_->computeIntegral(grad[0],wvel,target_);
      
      return g_param;
    }
    else {
      throw Exception::Zero(">>> QoI_Velocity_Darcy::gradient_3 is zero.");
    }
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c = fePrs_->gradN().extent_int(0);
    const int f = fePrs_->gradN().extent_int(1);
    const int p = fePrs_->gradN().extent_int(2);
    const int d = fePrs_->gradN().extent_int(3);
    // Initialize output hess
    hess = scalar_view("hess", c, f);
    // Compute cost integral
    scalar_view gradV("gradV",c,p,d);
    scalar_view awvel("awvel",c,p,d);
    scalar_view valZ("valZ",c,p);
    scalar_view alpha("alpha",c,p);
    fePrs_->evaluateGradient(gradV, v_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);

    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          awvel(i,j,k) = alpha(i,j)*alpha(i,j)*weight_(i,j,k)*gradV(i,j,k);
        }
      }
    }

    fst::integrate(hess,awvel,fePrs_->gradNdetJ(),false);
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c = fePrs_->gradN().extent_int(0);
    const int f = fePrs_->gradN().extent_int(1);
    const int p = fePrs_->gradN().extent_int(2);
    const int d = fePrs_->gradN().extent_int(3);
    Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
    // Initialize output grad
    hess = scalar_view("hess", c, f);
    // Compute cost integral
    scalar_view gradU("gradU",c,p,d);
    scalar_view awvel("awvel",c,p,d);
    scalar_view valZ("valZ",c,p);
    scalar_view valV("valV",c,p);
    scalar_view alpha("alpha",c,p);
    scalar_view dalpha("dalpha",c,p);
    fePrs_->evaluateGradient(gradU, u_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    feCtrl_->evaluateValue(valV, v_coeff);
    perm_->compute( alpha, valZ, fePrs_->cubPts(), 0);
    perm_->compute(dalpha, valZ, fePrs_->cubPts(), 1);

    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          awvel(i,j,k)  = static_cast<Real>(2)*alpha(i,j)*gradU(i,j,k)+yvel*target_(i,j,k);
          awvel(i,j,k) *= weight_(i,j,k)*dalpha(i,j)*valV(i,j);
        }
      }
    }

    fst::integrate(hess,awvel,fePrs_->gradNdetJ(),false);
  }

  void HessVec_13(scalar_view & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN().extent_int(0);
      const int f = fePrs_->gradN().extent_int(1);
      const int p = fePrs_->gradN().extent_int(2);
      const int d = fePrs_->gradN().extent_int(3);
      // Initialize output val
      hess = scalar_view("hess", c, f);
      // Compute cost integral
      scalar_view wvel("wvel",c,p,d);
      scalar_view valZ("valZ",c,p);
      scalar_view alpha("alpha",c,p);
      feCtrl_->evaluateValue(valZ, z_coeff);
      perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);

      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          for (int k = 0; k < d; ++k) {
            wvel(i,j,k) = weight_(i,j,k)*alpha(i,j)*target_(i,j,k)*(*v_param)[0];
          }
        }
      }

      fst::integrate(hess,wvel,fePrs_->gradNdetJ(),false);
    }
    else {
      throw Exception::NotImplemented(">>> HessVec_13 not implemented.");
    }
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c  = fePrs_->gradN().extent_int(0);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int p  = fePrs_->gradN().extent_int(2);
    const int d  = fePrs_->gradN().extent_int(3);
    Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
    // Initialize output grad
    hess = scalar_view("hess", c, fc);
    // Compute cost integral
    scalar_view gradU("gradU",c,p,d);
    scalar_view gradV("gradV",c,p,d);
    scalar_view valZ("valZ",c,p);
    scalar_view alpha("alpha",c,p);
    scalar_view dalpha("dalpha",c,p);
    scalar_view deriv("deriv",c,p);
    fePrs_->evaluateGradient(gradU, u_coeff);
    fePrs_->evaluateGradient(gradV, v_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);
    perm_->compute(dalpha, valZ, fePrs_->cubPts(), 1);

    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          deriv(i,j) += weight_(i,j,k)*(alpha(i,j)*gradU(i,j,k)+yvel*target_(i,j,k))
                        *dalpha(i,j)*gradV(i,j,k);
          deriv(i,j) += weight_(i,j,k)*dalpha(i,j)*gradU(i,j,k)*alpha(i,j)*gradV(i,j,k);
        }
      }
    }

    fst::integrate(hess,deriv,feCtrl_->NdetJ(),false);
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c  = fePrs_->gradN().extent_int(0);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int p  = fePrs_->gradN().extent_int(2);
    const int d  = fePrs_->gradN().extent_int(3);
    Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
    // Initialize output hess
    hess = scalar_view("hess", c, fc);
    // Compute cost integral
    scalar_view gradU("gradU",c,p,d);
    scalar_view valZ("valZ",c,p);
    scalar_view valV("valV",c,p);
    scalar_view alpha("alpha",c,p);
    scalar_view dalpha("dalpha",c,p);
    scalar_view ddalpha("ddalpha",c,p);
    scalar_view deriv("deriv",c,p);
    fePrs_->evaluateGradient(gradU, u_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    feCtrl_->evaluateValue(valV, v_coeff);
    perm_->compute(  alpha, valZ, fePrs_->cubPts(), 0);
    perm_->compute( dalpha, valZ, fePrs_->cubPts(), 1);
    perm_->compute(ddalpha, valZ, fePrs_->cubPts(), 2);

    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          deriv(i,j) += weight_(i,j,k)*(alpha(i,j)*gradU(i,j,k)+yvel*target_(i,j,k))
                        *ddalpha(i,j)*valV(i,j)*gradU(i,j,k);
          deriv(i,j) += weight_(i,j,k)*dalpha(i,j)*gradU(i,j,k)
                        *dalpha(i,j)*valV(i,j)*gradU(i,j,k);
        }
      }
    }

    fst::integrate(hess,deriv,feCtrl_->NdetJ(),false);
  }

  void HessVec_23(scalar_view & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c  = fePrs_->gradN().extent_int(0);
      const int fc = feCtrl_->N().extent_int(1);
      const int p  = fePrs_->gradN().extent_int(2);
      const int d  = fePrs_->gradN().extent_int(3);
      // Initialize output val
      hess = scalar_view("hess", c, fc);
      // Compute cost integral
      scalar_view gradU("gradU",c,p,d);
      scalar_view wvel("wvel",c,p);
      scalar_view valZ("valZ",c,p);
      scalar_view dalpha("dalpha",c,p);
      fePrs_->evaluateGradient(gradU, u_coeff);
      feCtrl_->evaluateValue(valZ, z_coeff);
      perm_->compute(dalpha, valZ, fePrs_->cubPts(), 1);

      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          for (int k = 0; k < d; ++k) {
            wvel(i,j) += weight_(i,j,k)*dalpha(i,j)*gradU(i,j,k)*target_(i,j,k)*(*v_param)[0];
          }
        }
      }

      fst::integrate(hess,wvel,feCtrl_->NdetJ(),false);
    }
    else {
      throw Exception::Zero(">>> QoI_Velocity_Darcy::HessVec_23 is zero.");
    }
  }

  std::vector<Real> HessVec_31(std::vector<scalar_view> & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN().extent_int(0);
      const int p = fePrs_->gradN().extent_int(2);
      const int d = fePrs_->gradN().extent_int(3);
      // Initialize output val
      std::vector<Real> h_param(d,static_cast<Real>(0));
      hess.clear(); hess.resize(d);
      for (int k = 0; k < d; ++k)
        hess[k] = scalar_view("hess",c);
      // Compute cost integral
      scalar_view gradV("gradV",c,p,d);
      scalar_view wvel("wvel",c,p,d);
      scalar_view valZ("valZ",c,p);
      scalar_view alpha("alpha",c,p);
      fePrs_->evaluateGradient(gradV, v_coeff);
      feCtrl_->evaluateValue(valZ, z_coeff);
      perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);

      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          for (int k = 0; k < d; ++k) {
            wvel(i,j,k) = weight_(i,j,k)*alpha(i,j)*gradV(i,j,k);
          }
        }
      }

      fePrs_->computeIntegral(hess[0],wvel,target_);

      return h_param;
    }
    else {
      throw Exception::Zero(">>> QoI_Velocity_Darcy::HessVec_31 is zero.");
    }
  }

  std::vector<Real> HessVec_32(std::vector<scalar_view> & hess,
                               const scalar_view v_coeff,
                               const scalar_view u_coeff,
                               const scalar_view z_coeff = ROL::nullPtr,
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN().extent_int(0);
      const int p = fePrs_->gradN().extent_int(2);
      const int d = fePrs_->gradN().extent_int(3);
      // Initialize output val
      std::vector<Real> h_param(d,static_cast<Real>(0));
      hess.clear(); hess.resize(d);
      for (int k = 0; k < d; ++k)
        hess[k] = scalar_view("hess",c);
      // Compute cost integral
      scalar_view gradU("gradU",c,p,d);
      scalar_view wvel("wvel",c,p,d);
      scalar_view valZ("valZ",c,p);
      scalar_view valV("valV",c,p);
      scalar_view dalpha("alpha",c,p);
      fePrs_->evaluateGradient(gradU, u_coeff);
      feCtrl_->evaluateValue(valZ, z_coeff);
      feCtrl_->evaluateValue(valV, v_coeff);
      perm_->compute(dalpha, valZ, fePrs_->cubPts(), 1);

      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          for (int k = 0; k < d; ++k) {
            wvel(i,j,k) = weight_(i,j,k)*dalpha(i,j)*gradU(i,j,k)*valV(i,j);
          }
        }
      }

      fePrs_->computeIntegral(hess[0],wvel,target_);

      return h_param;
    }
    else {
      throw Exception::Zero(">>> HessVec_32 is zero.");
    }
  }

  std::vector<Real> HessVec_33(std::vector<scalar_view> & hess,
                               const ROL::Ptr<const std::vector<Real>> & v_param,
                               const scalar_view u_coeff,
                               const scalar_view z_coeff = ROL::nullPtr,
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN().extent_int(0);
      const int p = fePrs_->gradN().extent_int(2);
      const int d = fePrs_->gradN().extent_int(3);
      // Initialize output val
      std::vector<Real> h_param(d,static_cast<Real>(0));
      hess.clear(); hess.resize(d);
      for (int k = 0; k < d; ++k)
        hess[k] = scalar_view("hess",c);
      // Compute cost integral
      scalar_view wtarget("wtarget",c,p,d);
      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          for (int k = 0; k < d; ++k) {
            wtarget(i,j,k) = weight_(i,j,k)*target_(i,j,k);
          }
        }
      }
      fePrs_->computeIntegral(hess[0],wtarget,target_);
      rst::scale(hess[0],(*v_param)[0]);
      
      return h_param;
    }
    else {
      throw Exception::Zero(">>> QoI_Velocity_Darcy::HessVec_33 is zero.");
    }
  }

}; // QoI_Velocity_Darcy
#endif
