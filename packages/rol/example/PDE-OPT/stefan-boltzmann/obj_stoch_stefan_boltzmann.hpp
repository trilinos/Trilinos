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

#ifndef PDEOPT_QOI_STOCHASTIC_STEFAN_BOLTZMANN_HPP
#define PDEOPT_QOI_STOCHASTIC_STEFAN_BOLTZMANN_HPP

#include "../TOOLS/qoi.hpp"
#include "pde_stoch_stefan_boltzmann.hpp"

template <class Real>
class QoI_StateCost : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real> > fe_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > weight_;
  Real T_;
  Real rate_;
  bool exp_;
  Real xmid_;
  Real vol_;

  Real weightFunc(const std::vector<Real> & x) const {
    return ((x[1] < xmid_) ? static_cast<Real>(0) : static_cast<Real>(1));
  }

  Real cost(const Real u, const int deriv = 0) const {
    if ( exp_ ) {
      if ( deriv == 1 ) {
        return std::exp(rate_ * (u - T_)) / vol_;
      }
      if ( deriv == 2 ) {
        return rate_ * std::exp(rate_ * (u - T_)) / vol_;
      }
      return (std::exp(rate_ * (u - T_)) - static_cast<Real>(1)) / (rate_ * vol_);
    }
    else {    
      if ( u > T_ ) {
        if ( deriv == 0 ) {
          return static_cast<Real>(0.5)*(u-T_)*(u-T_)/vol_;
        }
        if ( deriv == 1 ) {
          return (u-T_)/vol_;
        }
        if ( deriv == 2 ) {
          return static_cast<Real>(1)/vol_;
        }
      }
      return static_cast<Real>(0);
    }
  }

public:
  QoI_StateCost(const ROL::Ptr<FE<Real> > &fe,
                Teuchos::ParameterList &parlist) : fe_(fe) {
    T_    = parlist.sublist("Problem").get("Desired engine temperature",373.0);
    exp_  = parlist.sublist("Problem").get("Use exponential engine cost model",false);
    rate_ = parlist.sublist("Problem").get("Exponential engine cost rate",1.0);
    
    xmid_ = parlist.sublist("Geometry").get("Step height",0.5);
    Real height = parlist.sublist("Geometry").get("Channel height",1.0);
    Real width  = parlist.sublist("Geometry").get("Channel width",5.0);
    vol_ = (height - xmid_) * width;

    int c = fe_->cubPts()->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int d = fe_->cubPts()->dimension(2);
    std::vector<Real> pt(d);
    weight_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (*fe_->cubPts())(i,j,k);
        }
        (*weight_)(i,j) = weightFunc(pt);
      }
    }
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int p = fe_->gradN()->dimension(2);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Evaluate state on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real> > valU_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute cost
    ROL::Ptr<Intrepid::FieldContainer<Real> > costU
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*costU)(i,j) = cost((*valU_eval)(i,j),0);
      }
    }
    // Compute squared L2-norm of diff
    fe_->computeIntegral(val,costU,weight_);
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    // Initialize output grad
    grad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Evaluate state on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real> > valU_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute cost
    ROL::Ptr<Intrepid::FieldContainer<Real> > costU
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*costU)(i,j) = cost((*valU_eval)(i,j),1);
      }
    }
    // Multiply cost and weight
    Intrepid::FieldContainer<Real> WcostU(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(WcostU,
                                                               *weight_,
                                                               *costU);
    // Compute gradient of squared L2-norm of diff
    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  WcostU,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_StateCost::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    // Initialize output grad
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Evaluate state on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real> > valU_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    // Evaluate direction on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real> > valV_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valV_eval, v_coeff);
    // Compute cost
    ROL::Ptr<Intrepid::FieldContainer<Real> > costU
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*costU)(i,j) = cost((*valU_eval)(i,j),2);
      }
    }
    // Multiply cost and weight
    Intrepid::FieldContainer<Real> WcostU(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(WcostU,
                                                               *weight_,
                                                               *costU);
    // Multiply weighted cost and direction
    Intrepid::FieldContainer<Real> WcostUV(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(WcostUV,
                                                               WcostU,
                                                               *valV_eval);
    // Compute gradient of squared L2-norm of diff
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  WcostUV,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_StateCost::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_StateCost::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_StateCost::HessVec_22 is zero.");
  }

}; // QoI_L2Tracking

template <class Real>
class QoI_ControlCost : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real> > fe_vol_;
  const std::vector<ROL::Ptr<FE<Real> > > fe_bdry_;
  const std::vector<std::vector<int> > bdryCellLocIds_;
  std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > weight_;
  Real T_;
  bool exp_;
  Real rate1_;
  Real rate2_;
  Real area_;

  Real weightFunc(const std::vector<Real> & x) const {
    return static_cast<Real>(1);
  }

  Real cost(const Real z, const int deriv = 0) const {
    if ( exp_ ) {
      if ( deriv == 1 ) {
        return (std::exp(rate1_ * (z - T_)) - std::exp(rate2_ * (T_ - z))) / area_;
      }
      if ( deriv == 2 ) {
        return (rate1_ * std::exp(rate1_ * (z - T_))
              + rate2_ * std::exp(rate2_ * (T_ - z))) / area_;
      }
      return (std::exp(rate1_ * (z - T_)) - static_cast<Real>(1)) / (rate1_ * area_)
            +(std::exp(rate2_ * (T_ - z)) - static_cast<Real>(1)) / (rate2_ * area_);
    }
    else {
      if ( deriv == 0 ) {
        return static_cast<Real>(0.5)*(z-T_)*(z-T_)/area_;
      }
      if ( deriv == 1 ) {
        return (z-T_)/area_;
      }
      if ( deriv == 2 ) {
        return static_cast<Real>(1)/area_;
      }
      return static_cast<Real>(0);
    }
  }

  ROL::Ptr<Intrepid::FieldContainer<Real> > getBoundaryCoeff(
      const Intrepid::FieldContainer<Real> & cell_coeff,
      int locSideId) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[locSideId];
    const int numCellsSide = bdryCellLocId.size();
    const int f = fe_vol_->N()->dimension(1);
    
    ROL::Ptr<Intrepid::FieldContainer<Real > > bdry_coeff = 
      ROL::makePtr<Intrepid::FieldContainer<Real >>(numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        (*bdry_coeff)(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

public:
  QoI_ControlCost(const ROL::Ptr<FE<Real> > &fe_vol,
                  const std::vector<ROL::Ptr<FE<Real> > > &fe_bdry,
                  const std::vector<std::vector<int> > &bdryCellLocIds,
                  Teuchos::ParameterList &parlist)
    : fe_vol_(fe_vol), fe_bdry_(fe_bdry), bdryCellLocIds_(bdryCellLocIds) {
    T_     = parlist.sublist("Problem").get("Desired control temperature",293.0);
    exp_   = parlist.sublist("Problem").get("Use exponential control cost model",false);
    rate1_ = parlist.sublist("Problem").get("Exponential upper control cost rate",1.0);
    rate2_ = parlist.sublist("Problem").get("Exponential lower control cost rate",1.0);
    Real w1 = parlist.sublist("Geometry").get("Channel width",5.0);
    Real w2 = parlist.sublist("Geometry").get("Step width",1.0);
    area_ = w1 - w2;
    const int numLocSides = bdryCellLocIds_.size();
    const int d = fe_vol_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    weight_.resize(numLocSides);
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = fe_bdry_[l]->cubPts()->dimension(1);
        weight_[l] = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        for (int i = 0; i < numCellsSide; ++i) {
          for (int j = 0; j < numCubPerSide; ++j) {
            for (int k = 0; k < d; ++k) {
              pt[k] = (*fe_bdry_[l]->cubPts())(i,j,k);
            }
            (*weight_[l])(i,j) = weightFunc(pt);
          }
        }
      }
    }
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    if ( z_coeff != ROL::nullPtr ) {
      const int c = fe_vol_->gradN()->dimension(0);
      // Initialize output val
      val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
      // Compute cost integral
      const int numLocSides = bdryCellLocIds_.size();
      for (int l = 0; l < numLocSides; ++l) {
        const int numCellsSide  = bdryCellLocIds_[l].size();
        if ( numCellsSide ) {
          const int numCubPerSide = fe_bdry_[l]->cubPts()->dimension(1);
          // Evaluate control on FE basis
          ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff_bdry
            = getBoundaryCoeff(*z_coeff, l);
          ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          fe_bdry_[l]->evaluateValue(valZ_eval, z_coeff_bdry);
          // Compute cost
          ROL::Ptr<Intrepid::FieldContainer<Real> > costZ
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          for (int i = 0; i < numCellsSide; ++i) {
            for (int j = 0; j < numCubPerSide; ++j) {
              (*costZ)(i,j) = cost((*valZ_eval)(i,j),0);
            }
          }
          // Integrate cell cost
          ROL::Ptr<Intrepid::FieldContainer<Real> > intVal
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide);
          fe_bdry_[l]->computeIntegral(intVal,costZ,weight_[l]);
          // Add to integral value
          for (int i = 0; i < numCellsSide; ++i) {
            int cidx = bdryCellLocIds_[l][i];
            (*val)(cidx) += (*intVal)(i);
          }
        }
      }
    }
    else {
      throw Exception::Zero(">>> QoI_ControlCost: Value is zero.");
    }
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_ControlCost::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    if ( z_coeff != ROL::nullPtr ) {
      const int c = fe_vol_->gradN()->dimension(0);
      const int f = fe_vol_->gradN()->dimension(1);
      // Initialize output val
      grad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      // Compute cost integral
      const int numLocSides = bdryCellLocIds_.size();
      for (int l = 0; l < numLocSides; ++l) {
        const int numCellsSide  = bdryCellLocIds_[l].size();
        if ( numCellsSide ) {
          const int numCubPerSide = fe_bdry_[l]->cubPts()->dimension(1);
          // Evaluate control on FE basis
          ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff_bdry
            = getBoundaryCoeff(*z_coeff, l);
          ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          fe_bdry_[l]->evaluateValue(valZ_eval, z_coeff_bdry);
          // Compute cost
          ROL::Ptr<Intrepid::FieldContainer<Real> > costZ
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          for (int i = 0; i < numCellsSide; ++i) {
            for (int j = 0; j < numCubPerSide; ++j) {
              (*costZ)(i,j) = cost((*valZ_eval)(i,j),1);
            }
          }
          // Multiply cost and weight
          Intrepid::FieldContainer<Real> WcostZ(numCellsSide, numCubPerSide);
          Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(WcostZ,
                                                                     *(weight_[l]),
                                                                     *costZ);
          // Compute gradient of squared L2-norm of diff
          ROL::Ptr<Intrepid::FieldContainer<Real> > intGrad
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f);
          Intrepid::FunctionSpaceTools::integrate<Real>(*intGrad,
                                                        WcostZ,
                                                        *(fe_bdry_[l]->NdetJ()),
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
    else {
      throw Exception::Zero(">>> QoI_ControlCost: Gradient_2 is zero.");
    }
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_ControlCost::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_ControlCost::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_ControlCost::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    if ( z_coeff != ROL::nullPtr ) {
      const int c = fe_vol_->gradN()->dimension(0);
      const int f = fe_vol_->gradN()->dimension(1);
      // Initialize output val
      hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      // Compute cost integral
      const int numLocSides = bdryCellLocIds_.size();
      for (int l = 0; l < numLocSides; ++l) {
        const int numCellsSide  = bdryCellLocIds_[l].size();
        if ( numCellsSide ) {
          const int numCubPerSide = fe_bdry_[l]->cubPts()->dimension(1);
          // Evaluate control on FE basis
          ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff_bdry
            = getBoundaryCoeff(*z_coeff, l);
          ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          fe_bdry_[l]->evaluateValue(valZ_eval, z_coeff_bdry);
          // Evaluate direction on FE basis
          ROL::Ptr<Intrepid::FieldContainer<Real> > v_coeff_bdry
            = getBoundaryCoeff(*v_coeff, l);
          ROL::Ptr<Intrepid::FieldContainer<Real> > valV_eval
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          fe_bdry_[l]->evaluateValue(valV_eval, v_coeff_bdry);
          // Compute cost
          ROL::Ptr<Intrepid::FieldContainer<Real> > costZ
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          for (int i = 0; i < numCellsSide; ++i) {
            for (int j = 0; j < numCubPerSide; ++j) {
              (*costZ)(i,j) = cost((*valZ_eval)(i,j),2);
            }
          }
          // Multiply cost and weight
          Intrepid::FieldContainer<Real> WcostZ(numCellsSide, numCubPerSide);
          Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(WcostZ,
                                                                     *(weight_[l]),
                                                                     *costZ);
          // Multiply weighted cost and direction
          Intrepid::FieldContainer<Real> WcostZV(numCellsSide, numCubPerSide);
          Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(WcostZV,
                                                                     WcostZ,
                                                                     *valV_eval);
          // Compute hessian times a vector of cost
          ROL::Ptr<Intrepid::FieldContainer<Real> > intHess
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f);
          Intrepid::FunctionSpaceTools::integrate<Real>(*intHess,
                                                        WcostZV,
                                                        *(fe_bdry_[l]->NdetJ()),
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
    else {
      throw Exception::Zero(">>> QoI_ControlCost: HessVec_22 is zero.");
    }
  }

}; // QoI_L2Penalty

template <class Real>
class QoI_AdvectionCost : public QoI<Real> {
public:
  QoI_AdvectionCost(void) {}

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    if ( z_param != ROL::nullPtr ) {
      const int size = z_param->size();
      Real sum(0);
      for (int i = 0; i < size; ++i) {
        sum += std::pow((*z_param)[i], 2);
      }
      val = ROL::nullPtr;
      return static_cast<Real>(0.5)*sum;
    }
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_AdvectionCost::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_AdvectionCost::gradient_2 is zero.");
  }

  std::vector<Real> gradient_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    if ( z_param != ROL::nullPtr ) {
      return *z_param;
    }
    else {
      throw Exception::Zero(">>> QoI_AdvectionCost::gradient_3 is zero.");
    }
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_AdvectionCost::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_AdvectionCost::HessVec_12 is zero.");
  }

  void HessVec_13(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const std::vector<Real> > & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_AdvectionCost::HessVec_13 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_AdvectionCost::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_AdvectionCost::HessVec_22 is zero.");
  }

  void HessVec_23(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const std::vector<Real> > & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_AdvectionCost::HessVec_23 is zero.");
  }

  std::vector<Real> HessVec_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_AdvectionCost::HessVec_31 is zero.");
  }

  std::vector<Real> HessVec_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_AdvectionCost::HessVec_32 is zero.");
  }

  std::vector<Real> HessVec_33(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const std::vector<Real> > & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    if ( z_param != ROL::nullPtr ) {
      return *v_param;
    }
    else {
      throw Exception::Zero(">>> QoI_AdvectionCost::HessVec_33 is zero.");
    }
  }

}; // QoI_Advection_Cost

template <class Real>
class StochasticStefanBoltzmannStdObjective : public ROL::StdObjective<Real> {
private:
  Real alpha0_;
  Real alpha1_;

public:
  StochasticStefanBoltzmannStdObjective(Teuchos::ParameterList &parlist) {
    alpha0_ = parlist.sublist("Problem").get("State Cost",1.0);
    alpha1_ = parlist.sublist("Problem").get("Control Cost",1.0);
  }

  Real value(const std::vector<Real> &x, Real &tol) {
    return alpha0_*x[0] + alpha1_*x[1];
  }

  void gradient(std::vector<Real> &g, const std::vector<Real> &x, Real &tol) {
    g[0] = alpha0_;
    g[1] = alpha1_;
  }

  void hessVec(std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    const Real zero(0);
    hv[0] = zero;
    hv[1] = zero;
  }

}; // OBJ_SCALAR

template <class Real>
class StochasticStefanBoltzmannStdObjective3 : public ROL::StdObjective<Real> {
private:
  Real alpha0_;
  Real alpha1_;
  Real alpha2_;

public:
  StochasticStefanBoltzmannStdObjective3(Teuchos::ParameterList &parlist) {
    alpha0_ = parlist.sublist("Problem").get("State Cost",1.0);
    alpha1_ = parlist.sublist("Problem").get("Control Cost",1.0);
    alpha2_ = parlist.sublist("Problem").get("Advection Cost",1.0);
  }

  Real value(const std::vector<Real> &x, Real &tol) {
    return alpha0_*x[0] + alpha1_*x[1] + alpha2_*x[2];
  }

  void gradient(std::vector<Real> &g, const std::vector<Real> &x, Real &tol) {
    g[0] = alpha0_;
    g[1] = alpha1_;
    g[2] = alpha2_;
  }

  void hessVec(std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    const Real zero(0);
    hv[0] = zero;
    hv[1] = zero;
    hv[2] = zero;
  }

}; // OBJ_SCALAR

#endif
