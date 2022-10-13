// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  obj.hpp
    \brief Provides the interface for local (cell-based) objective function computations.
*/

#ifndef OBJ_DARCY_HPP
#define OBJ_DARCY_HPP

#include "../../../../TOOLS/qoi.hpp"
#include "pde_darcy.hpp"
#include "permeability.hpp"

template <class Real>
class QoI_Velocity_Darcy : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>>                    fePrs_, feCtrl_;
  const std::vector<ROL::Ptr<FE<Real>>>       fePrsBdry_, feCtrlBdry_;
  const std::vector<std::vector<int>>         bdryCellLocIds_;
  const ROL::Ptr<Permeability<Real>>          perm_;
  std::vector<Real>                           target_;
  bool                                        onlyAxial_;

  ROL::Ptr<Intrepid::FieldContainer<Real>> getBoundaryCoeff(
      const Intrepid::FieldContainer<Real> & cell_coeff,
      int locSideId,
      const ROL::Ptr<FE<Real>> &fe) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[locSideId];
    const int numCellsSide = bdryCellLocId.size();
    const int f = fe->N()->dimension(1);
    
    ROL::Ptr<Intrepid::FieldContainer<Real>> bdry_coeff = 
      ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        (*bdry_coeff)(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

public:
  QoI_Velocity_Darcy(Teuchos::ParameterList   &list,
                      const ROL::Ptr<FE<Real>> &fePrs,
                      const ROL::Ptr<FE<Real>> &feCtrl,
                      const std::vector<ROL::Ptr<FE<Real>>> &fePrsBdry,
                      const std::vector<ROL::Ptr<FE<Real>>> &feCtrlBdry,
                      const std::vector<std::vector<int>>   &bdryCellLocIds,
                      const ROL::Ptr<Permeability<Real>>    &perm)
    : fePrs_(fePrs), feCtrl_(feCtrl),
      fePrsBdry_(fePrsBdry), feCtrlBdry_(feCtrlBdry),
      bdryCellLocIds_(bdryCellLocIds),
      perm_(perm) {
    target_.clear(); target_.resize(2);
    target_[0] = list.sublist("Problem").get("Target Radial Velocity",0.0);
    target_[1] = list.sublist("Problem").get("Target Axial Velocity",-15.0);
    onlyAxial_ = list.sublist("Problem").get("Only Use Axial Velocity",false);
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fePrs_->gradN()->dimension(0);
    const int d = fePrs_->gradN()->dimension(3);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Compute cost integral
    std::vector<Real> target = (z_param == ROL::nullPtr) ? target_ : *z_param;
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = fePrsBdry_[l]->cubPts()->dimension(1);
        ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff_bdry, z_coeff_bdry, valU_eval, valZ_eval, alpha, intVal;
        valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide, d);
        valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        alpha     = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        intVal    = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide);
        u_coeff_bdry = getBoundaryCoeff(*u_coeff, l, fePrs_);
        z_coeff_bdry = getBoundaryCoeff(*z_coeff, l, feCtrl_);
        fePrsBdry_[l]->evaluateGradient(valU_eval, u_coeff_bdry);
        feCtrlBdry_[l]->evaluateValue(valZ_eval, z_coeff_bdry);
        perm_->compute(alpha, valZ_eval, fePrsBdry_[l]->cubPts(), 0);
        for (int i = 0; i < numCellsSide; ++i) {
          for (int j = 0; j < numCubPerSide; ++j) {
            for (int k = 0; k < d; ++k) {
              (*valU_eval)(i,j,k) *= -(*alpha)(i,j);
              if (k==0 && onlyAxial_) (*valU_eval)(i,j,k) = static_cast<Real>(0);
              else                    (*valU_eval)(i,j,k) -= target[k];
            }
          }
        }
        fePrsBdry_[l]->computeIntegral(intVal,valU_eval,valU_eval,true);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          (*val)(cidx) += static_cast<Real>(0.5)*(*intVal)(i);
        }
      }
    }
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fePrs_->gradN()->dimension(0);
    const int f = fePrs_->gradN()->dimension(1);
    const int d = fePrs_->gradN()->dimension(3);
    // Initialize output grad
    grad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Compute cost integral
    std::vector<Real> target = (z_param == ROL::nullPtr) ? target_ : *z_param;
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = fePrsBdry_[l]->cubPts()->dimension(1);
        ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff_bdry, z_coeff_bdry, valU_eval, valZ_eval, alpha, intGrad;
        valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide, d);
        valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        alpha     = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        intGrad   = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f);
        u_coeff_bdry = getBoundaryCoeff(*u_coeff, l, fePrs_);
        z_coeff_bdry = getBoundaryCoeff(*z_coeff, l, feCtrl_);
        fePrsBdry_[l]->evaluateGradient(valU_eval, u_coeff_bdry);
        feCtrlBdry_[l]->evaluateValue(valZ_eval, z_coeff_bdry);
        perm_->compute(alpha, valZ_eval, fePrsBdry_[l]->cubPts(), 0);
        for (int i = 0; i < numCellsSide; ++i) {
          for (int j = 0; j < numCubPerSide; ++j) {
            for (int k = 0; k < d; ++k) {
              (*valU_eval)(i,j,k) *= -(*alpha)(i,j);
              if (k==0 && onlyAxial_) (*valU_eval)(i,j,k) = static_cast<Real>(0);
              else                    (*valU_eval)(i,j,k) -= target[k];
              (*valU_eval)(i,j,k) *= -(*alpha)(i,j);
            }
          }
        }
        // Compute gradient of squared L2-norm of diff
        Intrepid::FunctionSpaceTools::integrate<Real>(*intGrad,
                                                      *valU_eval,
                                                      *fePrsBdry_[l]->gradNdetJ(),
                                                      Intrepid::COMP_CPP, false);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          for (int j = 0; j < f; ++j) {
            (*grad)(cidx,j) += (*intGrad)(i,j);
          }
        }
      }
    }
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c  = fePrs_->gradN()->dimension(0);
    const int fc = feCtrl_->gradN()->dimension(1);
    const int d  = fePrs_->gradN()->dimension(3);
    // Initialize output grad
    grad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc);
    // Compute cost integral
    std::vector<Real> target = (z_param == ROL::nullPtr) ? target_ : *z_param;
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = fePrsBdry_[l]->cubPts()->dimension(1);
        ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff_bdry, z_coeff_bdry;
        ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval, valZ_eval, alpha, dalpha;
        ROL::Ptr<Intrepid::FieldContainer<Real>> integrand, intGrad;
        valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide, d);
        valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        alpha     = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        dalpha    = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        integrand = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        intGrad   = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, fc);
        u_coeff_bdry = getBoundaryCoeff(*u_coeff, l, fePrs_);
        z_coeff_bdry = getBoundaryCoeff(*z_coeff, l, feCtrl_);
        fePrsBdry_[l]->evaluateGradient(valU_eval, u_coeff_bdry);
        feCtrlBdry_[l]->evaluateValue(valZ_eval, z_coeff_bdry);
        perm_->compute(alpha, valZ_eval, fePrsBdry_[l]->cubPts(), 0);
        perm_->compute(dalpha, valZ_eval, fePrsBdry_[l]->cubPts(), 1);
        Real dalphaU(0), misfit(0);
        for (int i = 0; i < numCellsSide; ++i) {
          for (int j = 0; j < numCubPerSide; ++j) {
            for (int k = 0; k < d; ++k) {
              dalphaU = -(*dalpha)(i,j) * (*valU_eval)(i,j,k);
              if (k==0 && onlyAxial_) misfit = static_cast<Real>(0);
              else                    misfit = -(*alpha)(i,j) * (*valU_eval)(i,j,k) - target[k];
              (*integrand)(i,j) += dalphaU * misfit;
            }
          }
        }
        // Compute gradient of squared L2-norm of diff
        Intrepid::FunctionSpaceTools::integrate<Real>(*intGrad,
                                                      *integrand,
                                                      *feCtrlBdry_[l]->NdetJ(),
                                                      Intrepid::COMP_CPP, false);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          for (int j = 0; j < fc; ++j) {
            (*grad)(cidx,j) += (*intGrad)(i,j);
          }
        }
      }
    }
  }

  std::vector<Real> gradient_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & grad,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN()->dimension(0);
      const int d = fePrs_->gradN()->dimension(3);
      // Initialize output val
      std::vector<Real> g_param(d,static_cast<Real>(0));
      grad.clear(); grad.resize(d);
      for (int k = 0; k < d; ++k)
       grad[k] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
      // Compute cost integral
      const int numLocSides = bdryCellLocIds_.size();
      for (int l = 0; l < numLocSides; ++l) {
        const int numCellsSide  = bdryCellLocIds_[l].size();
        if ( numCellsSide ) {
          const int numCubPerSide = fePrsBdry_[l]->cubPts()->dimension(1);
          std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> intVal(d);
          ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff_bdry, z_coeff_bdry, valU_eval, valZ_eval, integrand, weight, alpha;
          valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide, d);
          valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          integrand = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          weight    = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          alpha     = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          u_coeff_bdry = getBoundaryCoeff(*u_coeff, l, fePrs_);
          z_coeff_bdry = getBoundaryCoeff(*z_coeff, l, feCtrl_);
          fePrsBdry_[l]->evaluateGradient(valU_eval, u_coeff_bdry);
          feCtrlBdry_[l]->evaluateValue(valZ_eval, z_coeff_bdry);
          perm_->compute(alpha, valZ_eval, fePrsBdry_[l]->cubPts(), 0);
          for (int k = 0; k < d; ++k) {
            intVal[k] = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide);
            if ((k==0 && !onlyAxial_) || k==1) {
              integrand->initialize();
              for (int i = 0; i < numCellsSide; ++i) {
                for (int j = 0; j < numCubPerSide; ++j) {
                  (*weight)(i,j) = static_cast<Real>(1);
                  (*integrand)(i,j) = (*z_param)[k] + (*alpha)(i,j) * (*valU_eval)(i,j,k);
                }
              }
              fePrsBdry_[l]->computeIntegral(intVal[k],weight,integrand,true);
            }
          }
          // Add to integral value
          for (int k = 0; k < d; ++k) {
            for (int i = 0; i < numCellsSide; ++i) {
              int cidx = bdryCellLocIds_[l][i];
              (*grad[k])(cidx) += (*intVal[k])(i);
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

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fePrs_->gradN()->dimension(0);
    const int f = fePrs_->gradN()->dimension(1);
    const int d = fePrs_->gradN()->dimension(3);
    // Initialize output grad
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Compute cost integral
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = fePrsBdry_[l]->cubPts()->dimension(1);
        ROL::Ptr<Intrepid::FieldContainer<Real>> v_coeff_bdry, z_coeff_bdry , valV_eval, valZ_eval, alpha, intHess;
        valV_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide, d);
        valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        alpha     = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        intHess   = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f);
        v_coeff_bdry = getBoundaryCoeff(*v_coeff, l, fePrs_);
        z_coeff_bdry = getBoundaryCoeff(*z_coeff, l, feCtrl_);
        fePrsBdry_[l]->evaluateGradient(valV_eval, v_coeff_bdry);
        feCtrlBdry_[l]->evaluateValue(valZ_eval, z_coeff_bdry);
        perm_->compute(alpha, valZ_eval, fePrsBdry_[l]->cubPts(), 0);
        for (int i = 0; i < numCellsSide; ++i) {
          for (int j = 0; j < numCubPerSide; ++j) {
            for (int k = 0; k < d; ++k) {
              (*valV_eval)(i,j,k) *= (-(*alpha)(i,j)) * (-(*alpha)(i,j));
                if (k==0 && onlyAxial_) (*valV_eval)(i,j,k) = static_cast<Real>(0);
            }
          }
        }
        // Compute hessian of squared L2-norm of diff
        Intrepid::FunctionSpaceTools::integrate<Real>(*intHess,
                                                      *valV_eval,
                                                      *(fePrsBdry_[l]->gradNdetJ()),
                                                      Intrepid::COMP_CPP, false);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          for (int j = 0; j < f; ++j) {
            (*hess)(cidx,j) += (*intHess)(i,j);
          }
        }
      }
    }
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c  = fePrs_->gradN()->dimension(0);
    const int f  = fePrs_->gradN()->dimension(1);
    const int fc = feCtrl_->gradN()->dimension(1);
    const int d  = fePrs_->gradN()->dimension(3);
    // Initialize output grad
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Compute cost integral
    std::vector<Real> target = (z_param == ROL::nullPtr) ? target_ : *z_param;
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = fePrsBdry_[l]->cubPts()->dimension(1);
        ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff_bdry, z_coeff_bdry, v_coeff_bdry;
        ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval, valZ_eval, valV_eval;
        ROL::Ptr<Intrepid::FieldContainer<Real>> alpha, dalpha, integrand, intHess;
        valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide, d);
        valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        valV_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        alpha     = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        dalpha    = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        integrand = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide, d);
        intHess   = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f);
        u_coeff_bdry = getBoundaryCoeff(*u_coeff, l, fePrs_);
        z_coeff_bdry = getBoundaryCoeff(*z_coeff, l, feCtrl_);
        v_coeff_bdry = getBoundaryCoeff(*v_coeff, l, feCtrl_);
        fePrsBdry_[l]->evaluateGradient(valU_eval, u_coeff_bdry);
        feCtrlBdry_[l]->evaluateValue(valZ_eval, z_coeff_bdry);
        feCtrlBdry_[l]->evaluateValue(valV_eval, v_coeff_bdry);
        perm_->compute(alpha, valZ_eval, fePrsBdry_[l]->cubPts(), 0);
        perm_->compute(dalpha, valZ_eval, fePrsBdry_[l]->cubPts(), 1);
        Real dalphaV(0), misfit(0), dmisfit(0);
        for (int i = 0; i < numCellsSide; ++i) {
          for (int j = 0; j < numCubPerSide; ++j) {
            for (int k = 0; k < d; ++k) {
              dalphaV = -(*dalpha)(i,j) * (*valV_eval)(i,j);
              if (k==0 && onlyAxial_) {
                misfit  = static_cast<Real>(0);
                dmisfit = static_cast<Real>(0);
              }
              else {
                misfit  = -(*alpha)(i,j)*(*valU_eval)(i,j,k) - target[k];
                dmisfit = -(*alpha)(i,j)*(*valU_eval)(i,j,k);
              }
              (*integrand)(i,j,k) = dalphaV * (misfit + dmisfit);
            }
          }
        }
        // Compute gradient of squared L2-norm of diff
        Intrepid::FunctionSpaceTools::integrate<Real>(*intHess,
                                                      *integrand,
                                                      *fePrsBdry_[l]->gradNdetJ(),
                                                      Intrepid::COMP_CPP, false);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          for (int j = 0; j < f; ++j) {
            (*hess)(cidx,j) += (*intHess)(i,j);
          }
        }
      }
    }
  }

  void HessVec_13(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN()->dimension(0);
      const int f = fePrs_->gradN()->dimension(1);
      const int d = fePrs_->gradN()->dimension(3);
      // Initialize output grad
      hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      // Compute cost integral
      const int numLocSides = bdryCellLocIds_.size();
      for (int l = 0; l < numLocSides; ++l) {
        const int numCellsSide  = bdryCellLocIds_[l].size();
        if ( numCellsSide ) {
          const int numCubPerSide = fePrsBdry_[l]->cubPts()->dimension(1);
          ROL::Ptr<Intrepid::FieldContainer<Real>> z_coeff_bdry, valZ_eval, alpha, weight, intHess;
          weight    = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide, d);
          valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          alpha     = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          intHess   = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f);
          z_coeff_bdry = getBoundaryCoeff(*z_coeff, l, feCtrl_);
          feCtrlBdry_[l]->evaluateValue(valZ_eval, z_coeff_bdry);
          perm_->compute(alpha, valZ_eval, fePrsBdry_[l]->cubPts(), 0);
          for (int k = 0; k < d; ++k) {
            if ((k==0 && !onlyAxial_) || k==1) {
              for (int i = 0; i < numCellsSide; ++i) {
                for (int j = 0; j < numCubPerSide; ++j) {
                  (*weight)(i,j,k) = (*v_param)[k] * (-(*alpha)(i,j));
                }
              }
              // Compute gradient of squared L2-norm of diff
              Intrepid::FunctionSpaceTools::integrate<Real>(*intHess,
                                                            *weight,
                                                            *fePrsBdry_[l]->gradNdetJ(),
                                                            Intrepid::COMP_CPP, false);
              // Add to integral value
              for (int i = 0; i < numCellsSide; ++i) {
                int cidx = bdryCellLocIds_[l][i];
                for (int j = 0; j < f; ++j) {
                  (*hess)(cidx,j) += (*intHess)(i,j);
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

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c  = fePrs_->gradN()->dimension(0);
    const int f  = fePrs_->gradN()->dimension(1);
    const int fc = feCtrl_->gradN()->dimension(1);
    const int d  = fePrs_->gradN()->dimension(3);
    // Initialize output grad
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc);
    // Compute cost integral
    std::vector<Real> target = (z_param == ROL::nullPtr) ? target_ : *z_param;
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = fePrsBdry_[l]->cubPts()->dimension(1);
        ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff_bdry, z_coeff_bdry, v_coeff_bdry;
        ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval, valZ_eval, valV_eval;
        ROL::Ptr<Intrepid::FieldContainer<Real>> alpha, dalpha, integrand, intHess;
        valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide, d);
        valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        valV_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide, d);
        alpha     = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        dalpha    = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        integrand = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        intHess   = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, fc);
        u_coeff_bdry = getBoundaryCoeff(*u_coeff, l, fePrs_);
        z_coeff_bdry = getBoundaryCoeff(*z_coeff, l, feCtrl_);
        v_coeff_bdry = getBoundaryCoeff(*v_coeff, l, fePrs_);
        fePrsBdry_[l]->evaluateGradient(valU_eval, u_coeff_bdry);
        feCtrlBdry_[l]->evaluateValue(valZ_eval, z_coeff_bdry);
        fePrsBdry_[l]->evaluateGradient(valV_eval, v_coeff_bdry);
        perm_->compute(alpha, valZ_eval, fePrsBdry_[l]->cubPts(), 0);
        perm_->compute(dalpha, valZ_eval, fePrsBdry_[l]->cubPts(), 1);
        Real dalphaV(0), misfit(0), dmisfit(0);
        for (int i = 0; i < numCellsSide; ++i) {
          for (int j = 0; j < numCubPerSide; ++j) {
            for (int k = 0; k < d; ++k) {
              dalphaV = -(*dalpha)(i,j) * (*valV_eval)(i,j,k);
              if (k==0 && onlyAxial_) {
                misfit  = static_cast<Real>(0);
                dmisfit = static_cast<Real>(0);
              }
              else {
                misfit  = -(*alpha)(i,j)*(*valU_eval)(i,j,k) - target[k];
                dmisfit = -(*alpha)(i,j)*(*valU_eval)(i,j,k);
              }
              (*integrand)(i,j) += dalphaV * (misfit + dmisfit);
            }
          }
        }
        // Compute gradient of squared L2-norm of diff
        Intrepid::FunctionSpaceTools::integrate<Real>(*intHess,
                                                      *integrand,
                                                      *feCtrlBdry_[l]->NdetJ(),
                                                      Intrepid::COMP_CPP, false);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          for (int j = 0; j < fc; ++j) {
            (*hess)(cidx,j) += (*intHess)(i,j);
          }
        }
      }
    }
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c  = fePrs_->gradN()->dimension(0);
    const int f  = fePrs_->gradN()->dimension(1);
    const int fc = feCtrl_->gradN()->dimension(1);
    const int d  = fePrs_->gradN()->dimension(3);
    // Initialize output grad
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc);
    // Compute cost integral
    std::vector<Real> target = (z_param == ROL::nullPtr) ? target_ : *z_param;
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = fePrsBdry_[l]->cubPts()->dimension(1);
        ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff_bdry, z_coeff_bdry, v_coeff_bdry;
        ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval, valZ_eval, valV_eval;
        ROL::Ptr<Intrepid::FieldContainer<Real>> alpha, dalpha, ddalpha, integrand, intHess;
        valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide, d);
        valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        valV_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        alpha     = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        dalpha    = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        ddalpha   = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        integrand = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        intHess   = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, fc);
        u_coeff_bdry = getBoundaryCoeff(*u_coeff, l, fePrs_);
        z_coeff_bdry = getBoundaryCoeff(*z_coeff, l, feCtrl_);
        v_coeff_bdry = getBoundaryCoeff(*v_coeff, l, feCtrl_);
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
              dalphaV  = -(*dalpha)(i,j) * (*valV_eval)(i,j) * (*valU_eval)(i,j,k);
              ddalphaV = -(*ddalpha)(i,j) * (*valV_eval)(i,j) * (*valU_eval)(i,j,k);
              if (k==0 && onlyAxial_) {
                misfit  = static_cast<Real>(0);
                dmisfit = static_cast<Real>(0);
              }
              else {
                misfit  = -(*alpha)(i,j)*(*valU_eval)(i,j,k) - target[k];
                dmisfit = -(*dalpha)(i,j)*(*valU_eval)(i,j,k);
              }
              (*integrand)(i,j) += dalphaV * dmisfit + ddalphaV * misfit;
            }
          }
        }
        // Compute gradient of squared L2-norm of diff
        Intrepid::FunctionSpaceTools::integrate<Real>(*intHess,
                                                      *integrand,
                                                      *feCtrlBdry_[l]->NdetJ(),
                                                      Intrepid::COMP_CPP, false);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          for (int j = 0; j < fc; ++j) {
            (*hess)(cidx,j) += (*intHess)(i,j);
          }
        }
      }
    }
  }

  void HessVec_23(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Velocity_Darcy::HessVec_23 is zero.");
  }

  std::vector<Real> HessVec_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN()->dimension(0);
      const int d = fePrs_->gradN()->dimension(3);
      // Initialize output val
      std::vector<Real> h_param(d,static_cast<Real>(0));
      hess.clear(); hess.resize(d);
      for (int k = 0; k < d; ++k)
        hess[k] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
      // Compute cost integral
      const int numLocSides = bdryCellLocIds_.size();
      for (int l = 0; l < numLocSides; ++l) {
        const int numCellsSide  = bdryCellLocIds_[l].size();
        if ( numCellsSide ) {
          const int numCubPerSide = fePrsBdry_[l]->cubPts()->dimension(1);
          std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> intVal(d);
          ROL::Ptr<Intrepid::FieldContainer<Real>> v_coeff_bdry, valV_eval, weight;
          valV_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          weight    = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          for (int k = 0; k < d; ++k) {
            intVal[k] = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide);
            if ((k==0 && !onlyAxial_) || k==1) {
              v_coeff_bdry = getBoundaryCoeff(*v_coeff, l, fePrs_);
              valV_eval->initialize();
              fePrsBdry_[l]->evaluateValue(valV_eval, v_coeff_bdry);
              for (int i = 0; i < numCellsSide; ++i) {
                for (int j = 0; j < numCubPerSide; ++j) {
                  (*weight)(i,j) = static_cast<Real>(1);
                  (*valV_eval)(i,j) *= static_cast<Real>(-1);
                }
              }
              fePrsBdry_[l]->computeIntegral(intVal[k],weight,valV_eval,true);
            }
          }
          // Add to integral value
          for (int k = 0; k < d; ++k) {
            for (int i = 0; i < numCellsSide; ++i) {
              int cidx = bdryCellLocIds_[l][i];
              (*hess[k])(cidx) += (*intVal[k])(i);
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

  std::vector<Real> HessVec_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> HessVec_32 is zero.");
  }

  std::vector<Real> HessVec_33(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                               const ROL::Ptr<const std::vector<Real>> & v_param,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN()->dimension(0);
      const int d = fePrs_->gradN()->dimension(3);
      // Initialize output val
      std::vector<Real> h_param(d,static_cast<Real>(0));
      hess.clear(); hess.resize(d);
      for (int k = 0; k < d; ++k)
        hess[k] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
      // Compute cost integral
      const int numLocSides = bdryCellLocIds_.size();
      for (int l = 0; l < numLocSides; ++l) {
        const int numCellsSide  = bdryCellLocIds_[l].size();
        if ( numCellsSide ) {
          const int numCubPerSide = fePrsBdry_[l]->cubPts()->dimension(1);
          std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> intVal(d);
          ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval, weight;
          valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          weight    = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          for (int k = 0; k < d; ++k) {
            intVal[k] = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide);
            if ((k==0 && !onlyAxial_) || k==1) {
              valU_eval->initialize();
              for (int i = 0; i < numCellsSide; ++i) {
                for (int j = 0; j < numCubPerSide; ++j) {
                  (*weight)(i,j) = static_cast<Real>(1);
                  (*valU_eval)(i,j) = (*v_param)[k];
                }
              }
              fePrsBdry_[l]->computeIntegral(intVal[k],weight,valU_eval,true);
            }
          }
          // Add to integral value
          for (int k = 0; k < d; ++k) {
            for (int i = 0; i < numCellsSide; ++i) {
              int cidx = bdryCellLocIds_[l][i];
              (*hess[k])(cidx) += (*intVal[k])(i);
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


template <class Real>
class QoI_VelocityTracking_Darcy : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>>                    fePrs_, feCtrl_;
  const ROL::Ptr<Permeability<Real>>          perm_;
  ROL::Ptr<Intrepid::FieldContainer<Real>>    target_, weight_;
  Real                                        rad_, yvel_, frac_, twpow_, top_;
  std::vector<Real>                           rvec_, zvec_;
  int                                         targetType_;
  bool                                        onlyAxial_, optimizeFrac_, invertFrac_, polyWeight_, veloWeight_;

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
  QoI_VelocityTracking_Darcy(Teuchos::ParameterList             &list,
                             const ROL::Ptr<FE<Real>>           &fePrs,
                             const ROL::Ptr<FE<Real>>           &feCtrl,
                             const ROL::Ptr<Permeability<Real>> &perm)
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
    const int c = fePrs_->gradN()->dimension(0);
    const int p = fePrs_->gradN()->dimension(2);
    target_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,2);
    weight_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,2);
    std::vector<Real> x(2);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        x[0] = (*fePrs_->cubPts())(i,j,0);
        x[1] = (*fePrs_->cubPts())(i,j,1);
        (*target_)(i,j,0) = xTarget(x);
        (*target_)(i,j,1) = yTarget(x);
        if (useNorm && yWeight(x)) {
          xWScal = static_cast<Real>(1)
                  /(std::pow((*target_)(i,j,0),2) + std::pow((*target_)(i,j,1),2));
          yWScal = xWScal;
        }
        (*weight_)(i,j,0) = x[0] * xWScal * xWeight(x);
        (*weight_)(i,j,1) = x[0] * yWScal * yWeight(x);
      }
    }
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fePrs_->gradN()->dimension(0);
    const int p = fePrs_->gradN()->dimension(2);
    const int d = fePrs_->gradN()->dimension(3);
    Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Compute cost integral
    ROL::Ptr<Intrepid::FieldContainer<Real>> gradU, valZ, alpha, vel, wvel;
    gradU      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
    vel        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
    wvel       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
    valZ       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    alpha      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    fePrs_->evaluateGradient(gradU, u_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);

    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          (*vel)(i,j,k)  = (*alpha)(i,j)*(*gradU)(i,j,k)+yvel*(*target_)(i,j,k);
          (*wvel)(i,j,k) = (*weight_)(i,j,k)*(*vel)(i,j,k);
        }
      }
    }

    fePrs_->computeIntegral(val,vel,wvel);
    // Scale by one half
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fePrs_->gradN()->dimension(0);
    const int f = fePrs_->gradN()->dimension(1);
    const int p = fePrs_->gradN()->dimension(2);
    const int d = fePrs_->gradN()->dimension(3);
    Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
    // Initialize output grad
    grad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Compute cost integral
    ROL::Ptr<Intrepid::FieldContainer<Real>> gradU, valZ, alpha, awvel;
    gradU      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
    awvel      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
    valZ       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    alpha      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    fePrs_->evaluateGradient(gradU, u_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);

    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          (*awvel)(i,j,k) = (*weight_)(i,j,k)*((*alpha)(i,j)*(*gradU)(i,j,k)+yvel*(*target_)(i,j,k))
                            *(*alpha)(i,j);
        }
      }
    }

    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  *awvel,
                                                  *(fePrs_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c  = fePrs_->gradN()->dimension(0);
    const int fc = feCtrl_->gradN()->dimension(1);
    const int p  = fePrs_->gradN()->dimension(2);
    const int d  = fePrs_->gradN()->dimension(3);
    Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
    // Initialize output grad
    grad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc);
    // Compute cost integral
    ROL::Ptr<Intrepid::FieldContainer<Real>> gradU, valZ, alpha, dalpha, deriv;
    gradU      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
    valZ       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    alpha      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    dalpha     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    deriv      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    fePrs_->evaluateGradient(gradU, u_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);
    perm_->compute(dalpha, valZ, fePrs_->cubPts(), 1);

    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          (*deriv)(i,j) += (*weight_)(i,j,k)*((*alpha)(i,j)*(*gradU)(i,j,k)+yvel*(*target_)(i,j,k))
                           *(*dalpha)(i,j)*(*gradU)(i,j,k);
        }
      }
    }

    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  *deriv,
                                                  *feCtrl_->NdetJ(),
                                                  Intrepid::COMP_CPP, false);
  }

  std::vector<Real> gradient_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & grad,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN()->dimension(0);
      const int p = fePrs_->gradN()->dimension(2);
      const int d = fePrs_->gradN()->dimension(3);
      Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
      // Initialize output val
      std::vector<Real> g_param(d,static_cast<Real>(0));
      grad.clear(); grad.resize(d);
      for (int i = 0; i < d; ++i)
        grad[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
      // Compute cost integral
      ROL::Ptr<Intrepid::FieldContainer<Real>> gradU, wvel, valZ, alpha;
      gradU      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
      wvel       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
      valZ       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
      alpha      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
      fePrs_->evaluateGradient(gradU, u_coeff);
      feCtrl_->evaluateValue(valZ, z_coeff);
      perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);

      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          for (int k = 0; k < d; ++k) {
            (*wvel)(i,j,k) = (*weight_)(i,j,k)*((*alpha)(i,j)*(*gradU)(i,j,k)+yvel*(*target_)(i,j,k));
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

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fePrs_->gradN()->dimension(0);
    const int f = fePrs_->gradN()->dimension(1);
    const int p = fePrs_->gradN()->dimension(2);
    const int d = fePrs_->gradN()->dimension(3);
    // Initialize output hess
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Compute cost integral
    ROL::Ptr<Intrepid::FieldContainer<Real>> gradV, awvel, valZ, alpha;
    gradV      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
    awvel      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
    valZ       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    alpha      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    fePrs_->evaluateGradient(gradV, v_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);

    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          (*awvel)(i,j,k) = (*alpha)(i,j)*(*alpha)(i,j)*(*weight_)(i,j,k)*(*gradV)(i,j,k);
        }
      }
    }

    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *awvel,
                                                  *(fePrs_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fePrs_->gradN()->dimension(0);
    const int f = fePrs_->gradN()->dimension(1);
    const int p = fePrs_->gradN()->dimension(2);
    const int d = fePrs_->gradN()->dimension(3);
    Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
    // Initialize output grad
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Compute cost integral
    ROL::Ptr<Intrepid::FieldContainer<Real>> gradU, valZ, alpha, awvel, dalpha, valV;
    gradU      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
    awvel      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
    valZ       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    valV       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    alpha      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    dalpha     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    fePrs_->evaluateGradient(gradU, u_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    feCtrl_->evaluateValue(valV, v_coeff);
    perm_->compute( alpha, valZ, fePrs_->cubPts(), 0);
    perm_->compute(dalpha, valZ, fePrs_->cubPts(), 1);

    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          (*awvel)(i,j,k)  = static_cast<Real>(2)*(*alpha)(i,j)*(*gradU)(i,j,k)+yvel*(*target_)(i,j,k);
          (*awvel)(i,j,k) *= (*weight_)(i,j,k)*(*dalpha)(i,j)*(*valV)(i,j);
        }
      }
    }

    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *awvel,
                                                  *fePrs_->gradNdetJ(),
                                                  Intrepid::COMP_CPP, false);
  }

  void HessVec_13(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN()->dimension(0);
      const int f = fePrs_->gradN()->dimension(1);
      const int p = fePrs_->gradN()->dimension(2);
      const int d = fePrs_->gradN()->dimension(3);
      // Initialize output val
      hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      // Compute cost integral
      ROL::Ptr<Intrepid::FieldContainer<Real>> wvel, valZ, alpha;
      wvel       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
      valZ       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
      alpha      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
      feCtrl_->evaluateValue(valZ, z_coeff);
      perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);

      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          for (int k = 0; k < d; ++k) {
            (*wvel)(i,j,k) = (*weight_)(i,j,k)*(*alpha)(i,j)*(*target_)(i,j,k)*(*v_param)[0];
          }
        }
      }

      Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                    *wvel,
                                                    *fePrs_->gradNdetJ(),
                                                    Intrepid::COMP_CPP, false);
    }
    else {
      throw Exception::NotImplemented(">>> HessVec_13 not implemented.");
    }
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c  = fePrs_->gradN()->dimension(0);
    const int fc = feCtrl_->gradN()->dimension(1);
    const int p  = fePrs_->gradN()->dimension(2);
    const int d  = fePrs_->gradN()->dimension(3);
    Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
    // Initialize output grad
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc);
    // Compute cost integral
    ROL::Ptr<Intrepid::FieldContainer<Real>> gradU, gradV, valZ, alpha, dalpha, deriv;
    gradU      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
    gradV      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
    valZ       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    alpha      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    dalpha     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    deriv      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    fePrs_->evaluateGradient(gradU, u_coeff);
    fePrs_->evaluateGradient(gradV, v_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);
    perm_->compute(dalpha, valZ, fePrs_->cubPts(), 1);

    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          (*deriv)(i,j) += (*weight_)(i,j,k)*((*alpha)(i,j)*(*gradU)(i,j,k)+yvel*(*target_)(i,j,k))
                           *(*dalpha)(i,j)*(*gradV)(i,j,k);
          (*deriv)(i,j) += (*weight_)(i,j,k)*(*dalpha)(i,j)*(*gradU)(i,j,k)*(*alpha)(i,j)*(*gradV)(i,j,k);
        }
      }
    }

    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *deriv,
                                                  *feCtrl_->NdetJ(),
                                                  Intrepid::COMP_CPP, false);
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c  = fePrs_->gradN()->dimension(0);
    const int fc = feCtrl_->gradN()->dimension(1);
    const int p  = fePrs_->gradN()->dimension(2);
    const int d  = fePrs_->gradN()->dimension(3);
    Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
    // Initialize output hess
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc);
    // Compute cost integral
    ROL::Ptr<Intrepid::FieldContainer<Real>> gradU, valZ, alpha, dalpha, ddalpha, valV, deriv;
    gradU      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
    valZ       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    valV       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    alpha      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    dalpha     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    ddalpha    = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    deriv      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    fePrs_->evaluateGradient(gradU, u_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    feCtrl_->evaluateValue(valV, v_coeff);
    perm_->compute(  alpha, valZ, fePrs_->cubPts(), 0);
    perm_->compute( dalpha, valZ, fePrs_->cubPts(), 1);
    perm_->compute(ddalpha, valZ, fePrs_->cubPts(), 2);

    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          (*deriv)(i,j) += (*weight_)(i,j,k)*((*alpha)(i,j)*(*gradU)(i,j,k)+yvel*(*target_)(i,j,k))
                           *(*ddalpha)(i,j)*(*valV)(i,j)*(*gradU)(i,j,k);
          (*deriv)(i,j) += (*weight_)(i,j,k)*(*dalpha)(i,j)*(*gradU)(i,j,k)
                           *(*dalpha)(i,j)*(*valV)(i,j)*(*gradU)(i,j,k);
        }
      }
    }

    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *deriv,
                                                  *feCtrl_->NdetJ(),
                                                  Intrepid::COMP_CPP, false);
  }

  void HessVec_23(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c  = fePrs_->gradN()->dimension(0);
      const int fc = feCtrl_->N()->dimension(1);
      const int p  = fePrs_->gradN()->dimension(2);
      const int d  = fePrs_->gradN()->dimension(3);
      // Initialize output val
      hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc);
      // Compute cost integral
      ROL::Ptr<Intrepid::FieldContainer<Real>> gradU, wvel, valZ, dalpha;
      gradU      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
      wvel       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
      valZ       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
      dalpha     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
      fePrs_->evaluateGradient(gradU, u_coeff);
      feCtrl_->evaluateValue(valZ, z_coeff);
      perm_->compute(dalpha, valZ, fePrs_->cubPts(), 1);

      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          for (int k = 0; k < d; ++k) {
            (*wvel)(i,j) += (*weight_)(i,j,k)*(*dalpha)(i,j)*(*gradU)(i,j,k)*(*target_)(i,j,k)*(*v_param)[0];
          }
        }
      }

      Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                    *wvel,
                                                    *feCtrl_->NdetJ(),
                                                    Intrepid::COMP_CPP, false);
    }
    else {
      throw Exception::Zero(">>> QoI_Velocity_Darcy::HessVec_23 is zero.");
    }
  }

  std::vector<Real> HessVec_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN()->dimension(0);
      const int p = fePrs_->gradN()->dimension(2);
      const int d = fePrs_->gradN()->dimension(3);
      // Initialize output val
      std::vector<Real> h_param(d,static_cast<Real>(0));
      hess.clear(); hess.resize(d);
      for (int k = 0; k < d; ++k)
        hess[k] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
      // Compute cost integral
      ROL::Ptr<Intrepid::FieldContainer<Real>> gradV, wvel, valZ, alpha;
      gradV      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
      wvel       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
      valZ       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
      alpha      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
      fePrs_->evaluateGradient(gradV, v_coeff);
      feCtrl_->evaluateValue(valZ, z_coeff);
      perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);

      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          for (int k = 0; k < d; ++k) {
            (*wvel)(i,j,k) = (*weight_)(i,j,k)*(*alpha)(i,j)*(*gradV)(i,j,k);
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

  std::vector<Real> HessVec_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN()->dimension(0);
      const int p = fePrs_->gradN()->dimension(2);
      const int d = fePrs_->gradN()->dimension(3);
      // Initialize output val
      std::vector<Real> h_param(d,static_cast<Real>(0));
      hess.clear(); hess.resize(d);
      for (int k = 0; k < d; ++k)
        hess[k] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
      // Compute cost integral
      ROL::Ptr<Intrepid::FieldContainer<Real>> gradU, wvel, valZ, valV, dalpha;
      gradU      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
      wvel       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
      valZ       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
      valV       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
      dalpha     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
      fePrs_->evaluateGradient(gradU, u_coeff);
      feCtrl_->evaluateValue(valZ, z_coeff);
      feCtrl_->evaluateValue(valV, v_coeff);
      perm_->compute(dalpha, valZ, fePrs_->cubPts(), 1);

      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          for (int k = 0; k < d; ++k) {
            (*wvel)(i,j,k) = (*weight_)(i,j,k)*(*dalpha)(i,j)*(*gradU)(i,j,k)*(*valV)(i,j);
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

  std::vector<Real> HessVec_33(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                               const ROL::Ptr<const std::vector<Real>> & v_param,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN()->dimension(0);
      const int p = fePrs_->gradN()->dimension(2);
      const int d = fePrs_->gradN()->dimension(3);
      // Initialize output val
      std::vector<Real> h_param(d,static_cast<Real>(0));
      hess.clear(); hess.resize(d);
      for (int k = 0; k < d; ++k)
        hess[k] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
      // Compute cost integral
      ROL::Ptr<Intrepid::FieldContainer<Real>> wtarget;
      wtarget = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          for (int k = 0; k < d; ++k) {
            (*wtarget)(i,j,k) = (*weight_)(i,j,k)*(*target_)(i,j,k);
          }
        }
      }
      fePrs_->computeIntegral(hess[0],wtarget,target_);
      Intrepid::RealSpaceTools<Real>::scale(*hess[0],(*v_param)[0]);
      
      return h_param;
    }
    else {
      throw Exception::Zero(">>> QoI_Velocity_Darcy::HessVec_33 is zero.");
    }
  }

}; // QoI_Velocity_Darcy
#endif
