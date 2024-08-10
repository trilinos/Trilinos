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

#ifndef PDEOPT_QOI_HELMHOLTZ_HPP
#define PDEOPT_QOI_HELMHOLTZ_HPP

#include "../TOOLS/qoi.hpp"
#include "pde_helmholtz.hpp"

template <class Real>
class QoI_Helmholtz_StateTracking : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real> > fe_;
  const ROL::Ptr<FieldHelper<Real> > fieldHelper_;

  ROL::Ptr<Intrepid::FieldContainer<Real> > weight_;
  std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > target_;

  Real RoiRadius_;
  Real waveNumber_;
  Real angle_;

protected:
  void computeDomainWeight(void) {
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);

    weight_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
   
    const Real zero(0), one(1);
    bool inside(false);
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      inside = false;
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (*fe_->cubPts())(i,j,k);
        }
        if ( insideDomain(x) ) {
          inside = true;
          break;
        }
      }
      for (int j = 0; j < p; ++j) {
        (*weight_)(i,j) = (inside ? one : zero);
      }
    }
  }

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
        (*target_[0])(i,j) = evaluateTarget(x,0);
        (*target_[1])(i,j) = evaluateTarget(x,1);
      }
    }
  }

  Real DegreesToRadians(const Real deg) const {
    return deg * static_cast<Real>(M_PI) / static_cast<Real>(180);
  }
  
public:
  QoI_Helmholtz_StateTracking(const ROL::Ptr<FE<Real> > &fe,
                                   const ROL::Ptr<FieldHelper<Real> > &fieldHelper,
                                   Teuchos::ParameterList &parlist)
    : fe_(fe), fieldHelper_(fieldHelper) {
    RoiRadius_  = parlist.sublist("Problem").get("ROI Radius",2.0);
    waveNumber_ = parlist.sublist("Problem").get("Wave Number",10.0);
    angle_      = parlist.sublist("Problem").get("Target Angle",45.0);
    angle_      = DegreesToRadians(angle_);
    computeDomainWeight();
    computeTarget();
  }

  virtual bool insideDomain(const std::vector<Real> &x) const {
    Real xnorm(0);
    const int d = x.size();
    for (int i = 0; i < d; ++i) {
      xnorm += x[i]*x[i];
    }
    xnorm = std::sqrt(xnorm);
    return (xnorm <= RoiRadius_);
  }

  virtual Real evaluateTarget(const std::vector<Real> &x, const int component) const {
    const Real arg = waveNumber_ * (std::cos(angle_)*x[0] + std::sin(angle_)*x[1]);
    return (component==0) ? std::cos(arg) : std::sin(arg);
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate tracking term
    ROL::Ptr<Intrepid::FieldContainer<Real> > diffU
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > WdiffU
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i=0; i<2; ++i) {
      diffU->initialize(); WdiffU->initialize(0);
      fe_->evaluateValue(diffU, U[i]);
      Intrepid::RealSpaceTools<Real>::subtract(*diffU,*target_[i]);
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*WdiffU,*weight_,*diffU);
      fe_->computeIntegral(val,WdiffU,diffU,true);
    }
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > G(2);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate tracking term
    ROL::Ptr<Intrepid::FieldContainer<Real> > diffU
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > WdiffU
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i=0; i<2; ++i) {
      diffU->initialize(); WdiffU->initialize(0);
      fe_->evaluateValue(diffU, U[i]);
      Intrepid::RealSpaceTools<Real>::subtract(*diffU,*target_[i]);
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*WdiffU,*weight_,*diffU);
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      Intrepid::FunctionSpaceTools::integrate<Real>(*G[i],
                                                    *WdiffU,
                                                    *(fe_->NdetJ()),
                                                    Intrepid::COMP_CPP,
                                                    false);
    }
    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz_StateTracking::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    // Initialize output hessvec
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > H(2);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Evaluate tracking term
    ROL::Ptr<Intrepid::FieldContainer<Real> > valV
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > WvalV
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i=0; i<2; ++i) {
      valV->initialize(); WvalV->initialize(0);
      fe_->evaluateValue(valV, V[i]);
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*WvalV,*weight_,*valV);
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      Intrepid::FunctionSpaceTools::integrate<Real>(*H[i],
                                                    *WvalV,
                                                    *(fe_->NdetJ()),
                                                    Intrepid::COMP_CPP,
                                                    false);
    }
    fieldHelper_->combineFieldCoeff(hess, H);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_22 is zero.");
  }

}; // QoI_Helmholtz


template <class Real>
class QoI_Helmholtz_ControlPenalty : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real> > fe_;
  const ROL::Ptr<FieldHelper<Real> > fieldHelper_;

  Real innerAnnulusRadius_;
  Real outerAnnulusRadius_;
  Real RoiRadius_;
  Real alpha_;

  ROL::Ptr<Intrepid::FieldContainer<Real> > weight_;

  void computeDomainWeight(void) {
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);

    weight_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
   
    const Real zero(0);
    bool inside(false);
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      inside = false;
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (*fe_->cubPts())(i,j,k);
        }
        if ( insideDomain(x) ) {
          inside = true;
          break;
        }
      }
      for (int j = 0; j < p; ++j) {
        (*weight_)(i,j) = (inside ? alpha_ : zero);
      }
    }
  }

public:
  QoI_Helmholtz_ControlPenalty(const ROL::Ptr<FE<Real> > &fe,
                               const ROL::Ptr<FieldHelper<Real> > &fieldHelper,
                               Teuchos::ParameterList &parlist)
  : fe_(fe), fieldHelper_(fieldHelper) {
    Real dist2annulus   = parlist.sublist("Problem").get("Distance to Control Annulus",0.5);
    Real annulusWidth   = parlist.sublist("Problem").get("Control Annulus Width",0.1);
    RoiRadius_          = parlist.sublist("Problem").get("ROI Radius",2.0);
    innerAnnulusRadius_ = RoiRadius_ + dist2annulus;
    outerAnnulusRadius_ = innerAnnulusRadius_ + annulusWidth;
    alpha_              = parlist.sublist("Problem").get("Control Penalty",1e-4);
    computeDomainWeight();
  }

  virtual bool insideDomain(const std::vector<Real> &x) const {
    Real xnorm(0);
    const int d = x.size();
    for (int i = 0; i < d; ++i) {
      xnorm += x[i]*x[i];
    }
    xnorm = std::sqrt(xnorm);
    return (xnorm <= outerAnnulusRadius_ && xnorm >= innerAnnulusRadius_);
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Evaluate control penalty
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > WvalZ
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i=0; i<2; ++i) {
      valZ->initialize(); WvalZ->initialize();
      fe_->evaluateValue(valZ, Z[i]);
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*WvalZ,*weight_,*valZ);
      fe_->computeIntegral(val,WvalZ,valZ,true);
    }
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz_ControlCost::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > G(2);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Evaluate control penalty
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > WvalZ
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i=0; i<2; ++i) {
      valZ->initialize(); WvalZ->initialize();
      fe_->evaluateValue(valZ, Z[i]);
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*WvalZ,*weight_,*valZ);
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      Intrepid::FunctionSpaceTools::integrate<Real>(*G[i],
                                                    *WvalZ,
                                                    *(fe_->NdetJ()),
                                                    Intrepid::COMP_CPP,
                                                    false);
    }
    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz_ControlCost::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz_ControlCost::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz_ControlCost::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    // Initialize output hessvec
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > H(2);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Evaluate control penalty
    ROL::Ptr<Intrepid::FieldContainer<Real> > valV
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > WvalV
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i=0; i<2; ++i) {
      valV->initialize(); WvalV->initialize();
      fe_->evaluateValue(valV, V[i]);
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*WvalV,*weight_,*valV);
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      Intrepid::FunctionSpaceTools::integrate<Real>(*H[i],
                                                    *WvalV,
                                                    *(fe_->NdetJ()),
                                                    Intrepid::COMP_CPP,
                                                    false);
    }
    fieldHelper_->combineFieldCoeff(hess, H);
  }

}; // QoI_Helmholtz_ControlPenalty

#endif
