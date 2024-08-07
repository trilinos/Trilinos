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

#ifndef PDEOPT_QOI_L2TRACKING_ALLEN_CAHN_HPP
#define PDEOPT_QOI_L2TRACKING_ALLEN_CAHN_HPP

#include "../TOOLS/qoi.hpp"
#include "pde_allen_cahn.hpp"

template <class Real>
class QoI_State_Cost_Allen_Cahn : public QoI<Real> {
private:
  ROL::Ptr<FE<Real> > fe_;

  ROL::Ptr<Intrepid::FieldContainer<Real> > target_;

  Real targetFunc(const std::vector<Real> & x) const {
    int size = x.size();
    Real val(0);
    for (int i = 0; i < size; ++i) {
      val += x[i]*x[i];
    }
    return val;
  }

public:
  QoI_State_Cost_Allen_Cahn(const ROL::Ptr<FE<Real> > &fe) : fe_(fe) {
    int c = fe_->cubPts()->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int d = fe_->cubPts()->dimension(2);
    std::vector<Real> pt(d);
    target_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (*fe_->cubPts())(i,j,k);
        }
        (*target_)(i,j) = targetFunc(pt);
      }
    }
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = u_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Evaluate state on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real> > valU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute difference between state and target
    Intrepid::RealSpaceTools<Real>::subtract(*valU_eval,*target_);
    // Compute squared L2-norm of diff
    fe_->computeIntegral(val,valU_eval,valU_eval);
    // Scale by one half
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = u_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int f = fe_->N()->dimension(1);
    // Initialize output grad
    grad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Evaluate state on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real> > valU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute difference between state and target
    Intrepid::RealSpaceTools<Real>::subtract(*valU_eval,*target_);
    // Compute gradient of squared L2-norm of diff
    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  *valU_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_State_Cost_Allen_Cahn::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    int c = v_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int f = fe_->N()->dimension(1);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valV_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    fe_->evaluateValue(valV_eval, v_coeff);
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *valV_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_State_Cost_Allen_Cahn::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_State_Cost_Allen_Cahn::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_State_Cost_Allen_Cahn::HessVec_22 is zero.");
  }

}; // QoI_L2Tracking

template <class Real>
class QoI_Control_Cost_Allen_Cahn : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real> > fe_vol_;
  const std::vector<std::vector<ROL::Ptr<FE<Real> > > > fe_bdry_;
  const std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;

  ROL::Ptr<Intrepid::FieldContainer<Real> > getBoundaryCoeff(
      const Intrepid::FieldContainer<Real> & cell_coeff,
      const int sideset, const int locSideId) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideset][locSideId];
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
  QoI_Control_Cost_Allen_Cahn(const ROL::Ptr<FE<Real> > &fe_vol,
                  const std::vector<std::vector<ROL::Ptr<FE<Real> > > > &fe_bdry,
                  const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds)
    : fe_vol_(fe_vol), fe_bdry_(fe_bdry), bdryCellLocIds_(bdryCellLocIds) {}

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    const int c = fe_vol_->gradN()->dimension(0);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Compute cost integral
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      const int numLocSides = bdryCellLocIds_[i].size();
      for (int j = 0; j < numLocSides; ++j) {
        const int numCellsSide  = bdryCellLocIds_[i][j].size();
        if ( numCellsSide ) {
          const int numCubPerSide = fe_bdry_[i][j]->cubPts()->dimension(1);
          // Evaluate control on FE basis
          ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff_bdry
            = getBoundaryCoeff(*z_coeff, i, j);
          ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          fe_bdry_[i][j]->evaluateValue(valZ_eval, z_coeff_bdry);
          // Integrate cell cost
          ROL::Ptr<Intrepid::FieldContainer<Real> > intVal
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide);
          fe_bdry_[i][j]->computeIntegral(intVal,valZ_eval,valZ_eval);
          // Add to integral value
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
            (*val)(cidx) += static_cast<Real>(0.5)*(*intVal)(k);
          }
        }
      }
    }
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_Allen_Cahn::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    const int c = fe_vol_->gradN()->dimension(0);
    const int f = fe_vol_->gradN()->dimension(1);
    // Initialize output val
    grad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Compute cost integral
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      const int numLocSides = bdryCellLocIds_[i].size();
      for (int j = 0; j < numLocSides; ++j) {
        const int numCellsSide  = bdryCellLocIds_[i][j].size();
        if ( numCellsSide ) {
          const int numCubPerSide = fe_bdry_[i][j]->cubPts()->dimension(1);
          // Evaluate control on FE basis
          ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff_bdry
            = getBoundaryCoeff(*z_coeff, i, j);
          ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          fe_bdry_[i][j]->evaluateValue(valZ_eval, z_coeff_bdry);
          // Compute gradient of squared L2-norm of diff
          ROL::Ptr<Intrepid::FieldContainer<Real> > intGrad
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f);
          Intrepid::FunctionSpaceTools::integrate<Real>(*intGrad,
                                                        *valZ_eval,
                                                        *(fe_bdry_[i][j]->NdetJ()),
                                                        Intrepid::COMP_CPP, false);
          // Add to integral value
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
            for (int l = 0; l < f; ++l) {
              (*grad)(cidx,l) += (*intGrad)(k,l);
            }
          }
        }
      }
    }
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_Allen_Cahn::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_Allen_Cahn::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_Allen_Cahn::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    const int c = fe_vol_->gradN()->dimension(0);
    const int f = fe_vol_->gradN()->dimension(1);
    // Initialize output val
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Compute cost integral
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      const int numLocSides = bdryCellLocIds_[i].size();
      for (int j = 0; j < numLocSides; ++j) {
        const int numCellsSide  = bdryCellLocIds_[i][j].size();
        if ( numCellsSide ) {
          const int numCubPerSide = fe_bdry_[i][j]->cubPts()->dimension(1);
          // Evaluate direction on FE basis
          ROL::Ptr<Intrepid::FieldContainer<Real> > v_coeff_bdry
            = getBoundaryCoeff(*v_coeff, i, j);
          ROL::Ptr<Intrepid::FieldContainer<Real> > valV_eval
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          fe_bdry_[i][j]->evaluateValue(valV_eval, v_coeff_bdry);
          // Compute hessian times a vector of cost
          ROL::Ptr<Intrepid::FieldContainer<Real> > intHess
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f);
          Intrepid::FunctionSpaceTools::integrate<Real>(*intHess,
                                                        *valV_eval,
                                                        *(fe_bdry_[i][j]->NdetJ()),
                                                        Intrepid::COMP_CPP, false);
          // Add to integral value
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
            for (int l = 0; l < f; ++l) {
              (*hess)(cidx,l) += (*intHess)(k,l);
            }
          }
        }
      }
    }
  }

}; // QoI_L2Penalty

#endif
