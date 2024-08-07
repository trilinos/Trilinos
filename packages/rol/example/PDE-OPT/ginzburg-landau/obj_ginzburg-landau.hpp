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

#ifndef PDEOPT_QOI_GINZBURGLANDAU_HPP
#define PDEOPT_QOI_GINZBURGLANDAU_HPP

#include "../TOOLS/qoi.hpp"
#include "pde_ginzburg-landau.hpp"

template <class Real>
class QoI_GinzburgLandau_StateTracking : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> fe_;
  const ROL::Ptr<FieldHelper<Real>> fieldHelper_;

  std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> target_;

  Real epsilon0_;
  Real lambda_;

protected:
  void computeTarget(void) {
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);

    target_.clear(); target_.resize(2);
    target_[0] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    target_[1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
   
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (*fe_->cubPts())(i,j,k);
        }
        (*target_[0])(i,j) = evaluateRealTarget(x);
        (*target_[1])(i,j) = evaluateImagTarget(x);
      }
    } 
  }
  
public:
  QoI_GinzburgLandau_StateTracking(const ROL::Ptr<FE<Real>> &fe,
                                   const ROL::Ptr<FieldHelper<Real>> &fieldHelper,
                                   Teuchos::ParameterList &parlist)
    : fe_(fe), fieldHelper_(fieldHelper) {
    lambda_   = parlist.sublist("Problem").get("Current Loading",1.0);
    epsilon0_ = parlist.sublist("Problem").get("State Scaling",1.0);
  }

  virtual Real evaluateRealTarget(const std::vector<Real> &x) const = 0;

  virtual Real evaluateImagTarget(const std::vector<Real> &x) const = 0;

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>>             &val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int p = fe_->gradN()->dimension(2);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval;
    valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i=0; i<2; ++i) {
      valU_eval->initialize();
      fe_->evaluateValue(valU_eval, U[i]);
      Intrepid::RealSpaceTools<Real>::subtract(*valU_eval,*target_[i]);
      fe_->computeIntegral(val,valU_eval,valU_eval,true);
    }
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5)*lambda_/epsilon0_);
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>>             &grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G(2);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval;
    valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i=0; i<2; ++i) {
      valU_eval->initialize();
      fe_->evaluateValue(valU_eval, U[i]);
      Intrepid::RealSpaceTools<Real>::subtract(*valU_eval,*target_[i]);
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      Intrepid::FunctionSpaceTools::integrate<Real>(*G[i],
                                                    *valU_eval,
                                                    *(fe_->NdetJ()),
                                                    Intrepid::COMP_CPP,
                                                    false);
      Intrepid::RealSpaceTools<Real>::scale(*G[i],lambda_/epsilon0_);
    }

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>>             &grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_GinzburgLandau_StateTracking::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>>             &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    // Initialize output hessvec
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H(2);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> valV_eval;
    valV_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i=0; i<2; ++i) {
      valV_eval->initialize();
      fe_->evaluateValue(valV_eval, V[i]);
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      Intrepid::FunctionSpaceTools::integrate<Real>(*H[i],
                                                    *valV_eval,
                                                    *(fe_->NdetJ()),
                                                    Intrepid::COMP_CPP,
                                                    false);
      Intrepid::RealSpaceTools<Real>::scale(*H[i],lambda_/epsilon0_);
    }

    fieldHelper_->combineFieldCoeff(hess, H);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>>             &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_TopoOpt::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>>             &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_TopoOpt::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>>             &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_TopoOpt::HessVec_22 is zero.");
  }

}; // QoI_TopoOpt


template <class Real>
class QoI_GinzburgLandau_ControlPenalty : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> fe_;
  const std::vector<ROL::Ptr<FE<Real>>> feBdry_;
  const std::vector<std::vector<int>> bdryCellLocIds_;
  const ROL::Ptr<FieldHelper<Real>> fieldHelper_;
  Real delta0_;
  Real lambda_;

  ROL::Ptr<Intrepid::FieldContainer<Real>> getBoundaryCoeff(
      const Intrepid::FieldContainer<Real> &cell_coeff,
      int locSideId) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[locSideId];
    const int numCellsSide = bdryCellLocId.size();
    const int f = fe_->N()->dimension(1);
    
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
  QoI_GinzburgLandau_ControlPenalty(const ROL::Ptr<FE<Real>> &fe,
                                    const std::vector<ROL::Ptr<FE<Real>>> &feBdry,
                                    const std::vector<std::vector<int>> &bdryCellLocIds,
                                    const ROL::Ptr<FieldHelper<Real>> &fieldHelper,
                                    Teuchos::ParameterList &parlist)
  : fe_(fe), feBdry_(feBdry), bdryCellLocIds_(bdryCellLocIds), fieldHelper_(fieldHelper) {
    delta0_ = parlist.sublist("Problem").get("Control Scaling",1.0);
    lambda_ = parlist.sublist("Problem").get("Current Loading",1.0);
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>>             &val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    const int c = fe_->gradN()->dimension(0);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Compute cost integral
    const int numLocSides = bdryCellLocIds_.size();
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < numLocSides; ++j) {
        const int numCellsSide = bdryCellLocIds_[j].size();
        if ( numCellsSide ) {
          const int numCubPerSide = feBdry_[j]->cubPts()->dimension(1);
          // Evaluate control on FE basis
          ROL::Ptr<Intrepid::FieldContainer<Real>> z_coeff_bdry
            = getBoundaryCoeff(*Z[i], j);
          ROL::Ptr<Intrepid::FieldContainer<Real>> valZ_eval
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          feBdry_[j]->evaluateValue(valZ_eval, z_coeff_bdry);
          // Integrate cell cost
          ROL::Ptr<Intrepid::FieldContainer<Real> > intVal
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide);
          feBdry_[j]->computeIntegral(intVal,valZ_eval,valZ_eval,false);
          // Add to integral value
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[j][k];
            (*val)(cidx) += static_cast<Real>(0.5)*delta0_*lambda_*(*intVal)(k);
          }
        }
      }
    }
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>>             &grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_GinzburgLandau_ControlCost::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>>             &grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G(2);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Compute cost integral
    const int numLocSides = bdryCellLocIds_.size();
    for (int i = 0; i < 2; ++i) {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      for (int j = 0; j < numLocSides; ++j) {
        const int numCellsSide  = bdryCellLocIds_[j].size();
        if ( numCellsSide ) {
          const int numCubPerSide = feBdry_[j]->cubPts()->dimension(1);
          // Evaluate control on FE basis
          ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff_bdry
            = getBoundaryCoeff(*Z[i], j);
          ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          feBdry_[j]->evaluateValue(valZ_eval, z_coeff_bdry);
          // Compute gradient of squared L2-norm
          Intrepid::FieldContainer<Real> intGrad(numCellsSide, f);
          Intrepid::FunctionSpaceTools::integrate<Real>(intGrad,
                                                        *valZ_eval,
                                                        *(feBdry_[j]->NdetJ()),
                                                        Intrepid::COMP_CPP, false);
          // Add to integral value
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[j][k];
            for (int l = 0; l < f; ++l) {
              (*G[i])(cidx,l) += lambda_*delta0_*intGrad(k,l);
            }
          }
        }
      }
    }

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>>             &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_GinzburgLandau_ControlCost::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>>             &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_GinzburgLandau_ControlCost::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>>             &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_GinzburgLandau_ControlCost::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>>             &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H(2);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Compute cost integral
    const int numLocSides = bdryCellLocIds_.size();
    for (int i = 0; i < 2; ++i) {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      for (int j = 0; j < numLocSides; ++j) {
        const int numCellsSide  = bdryCellLocIds_[j].size();
        if ( numCellsSide ) {
          const int numCubPerSide = feBdry_[j]->cubPts()->dimension(1);
          // Evaluate control on FE basis
          ROL::Ptr<Intrepid::FieldContainer<Real> > v_coeff_bdry
            = getBoundaryCoeff(*V[i], j);
          ROL::Ptr<Intrepid::FieldContainer<Real> > valV_eval
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          feBdry_[j]->evaluateValue(valV_eval, v_coeff_bdry);
          // Compute gradient of squared L2-norm of diff
          Intrepid::FieldContainer<Real> intHess(numCellsSide, f);
          Intrepid::FunctionSpaceTools::integrate<Real>(intHess,
                                                        *valV_eval,
                                                        *(feBdry_[j]->NdetJ()),
                                                        Intrepid::COMP_CPP, false);
          // Add to integral value
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[j][k];
            for (int l = 0; l < f; ++l) {
              (*H[i])(cidx,l) += lambda_*delta0_*intHess(k,l);
            }
          }
        }
      }
    }

    fieldHelper_->combineFieldCoeff(hess, H);
  }

}; // QoI_GinzburgLandau_ControlPenalty

#endif
