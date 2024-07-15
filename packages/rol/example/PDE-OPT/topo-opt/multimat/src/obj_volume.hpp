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

#ifndef PDEOPT_QOI_VOLUME_MULTIMAT_HPP
#define PDEOPT_QOI_VOLUME_MULTIMAT_HPP

#include "../../../TOOLS/qoi.hpp"
#include "pde_elasticity.hpp"

template <class Real>
class QoI_MultiMat_Weight : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> fe_;
  const ROL::Ptr<FieldUtils::FieldInfo> fieldInfo_;
  ROL::Ptr<Intrepid::FieldContainer<Real>> ones_;
  ROL::Ptr<Intrepid::FieldContainer<Real>> volFrac_;
  std::vector<Real> w_;

public:
  QoI_MultiMat_Weight(const ROL::Ptr<FE<Real>> &fe,
                      const ROL::Ptr<FieldUtils::FieldInfo> &fieldInfo,
                      ROL::ParameterList &list)
  : fe_(fe), fieldInfo_(fieldInfo) {
    Real w0 = list.sublist("Problem").get("Maximum Weight Fraction", 0.5);
    w_  = ROL::getArrayFromStringParameter<Real>(list.sublist("Problem"), "Density");
    Real mw(0);
    int T = w_.size();
    for (int i = 0; i < T; ++i) {
      mw = std::max(mw,w_[i]);
    }
    // Get relevant dimensions
    int c = fe_->cubPts()->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    ones_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    volFrac_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    ones_->initialize(static_cast<Real>(1));
    volFrac_->initialize(w0 * mw);
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int T = w_.size();

    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Split z_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> Z;
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfo_);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> valZ_eval, sumZ;
    valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    sumZ      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int t=0; t<T; ++t) {
      valZ_eval->initialize();
      fe_->evaluateValue(valZ_eval, Z[t]);
      // Compute weight
      Intrepid::RealSpaceTools<Real>::scale(*valZ_eval,w_[t]);
      Intrepid::RealSpaceTools<Real>::add(*sumZ,*valZ_eval);
    }
    Intrepid::RealSpaceTools<Real>::subtract(*sumZ,*volFrac_);
    fe_->computeIntegral(val,sumZ,ones_,true);
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Volume_TopoOpt::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int T = w_.size();

    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G(T);
    for (int t=0; t<T; ++t) {
      G[t] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      // Compute gradient of energy
      Intrepid::FunctionSpaceTools::integrate<Real>(*G[t],
                                                    *ones_,
                                                    *fe_->NdetJ(),
                                                    Intrepid::COMP_CPP, false);
      Intrepid::RealSpaceTools<Real>::scale(*G[t],w_[t]);
    }
    // Combine the gradients.
    FieldUtils::combineFieldCoeff<Real>(grad, G, fieldInfo_);
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Volume_TopoOpt::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Volume_TopoOpt::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Volume_TopoOpt::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Volume_TopoOpt::HessVec_22 is zero.");
  }

}; // QoI_MultiMat_Weight

#endif
