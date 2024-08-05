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

#include "../../TOOLS/qoi.hpp"
#include "pde_helmholtz_real.hpp"

template <class Real>
class QoI_Helmholtz_StateTracking : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> fe_;

  ROL::Ptr<Intrepid::FieldContainer<Real>> weight_;
  std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> target_;

  Real RoiRadius_;
  Real waveNumber_;
  Real angle_;
  int example_;

  unsigned comp_;

protected:
  bool insideDomain(const std::vector<Real> &x) const {
    bool val = true;
    if (example_==1) {
      const Real eps = std::sqrt(ROL::ROL_EPSILON<Real>());
      Real xnorm(0);
      const int d = x.size();
      for (int i = 0; i < d; ++i) {
        xnorm += x[i]*x[i];
      }
      xnorm = std::sqrt(xnorm);
      val = (xnorm <= RoiRadius_+eps);
    }
    return val;
  }

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

  Real evaluateTarget(const std::vector<Real> &x, const int component) const {
    const Real arg = waveNumber_ * (std::cos(angle_)*x[0] + std::sin(angle_)*x[1]);
    return (component==0) ? std::cos(arg) : std::sin(arg);
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
  QoI_Helmholtz_StateTracking(const ROL::Ptr<FE<Real>> &fe,
                              Teuchos::ParameterList &parlist,
                              const int comp = 0)
    : fe_(fe), RoiRadius_(2), comp_(comp) {
    waveNumber_ = parlist.sublist("Problem").get("Wave Number",10.0);
    angle_      = parlist.sublist("Problem").get("Target Angle",45.0);
    angle_      = DegreesToRadians(angle_);
    example_    = parlist.sublist("Problem").get("Example",1);
    computeDomainWeight();
    computeTarget();
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Evaluate tracking term
    ROL::Ptr<Intrepid::FieldContainer<Real>> diffU, WdiffU;
    diffU  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    WdiffU = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    diffU->initialize(); WdiffU->initialize(0);
    fe_->evaluateValue(diffU, u_coeff);
    Intrepid::RealSpaceTools<Real>::subtract(*diffU,*target_[comp_]);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*WdiffU,*weight_,*diffU);
    fe_->computeIntegral(val,WdiffU,diffU,true);
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    // Evaluate tracking term
    ROL::Ptr<Intrepid::FieldContainer<Real>> diffU, WdiffU;
    diffU  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    WdiffU = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    diffU->initialize(); WdiffU->initialize(0);
    fe_->evaluateValue(diffU, u_coeff);
    Intrepid::RealSpaceTools<Real>::subtract(*diffU,*target_[comp_]);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*WdiffU,*weight_,*diffU);
    grad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  *WdiffU,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP,
                                                  false);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz_StateTracking::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    // Evaluate tracking term
    ROL::Ptr<Intrepid::FieldContainer<Real>> valV, WvalV;
    valV  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    WvalV = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    valV->initialize(); WvalV->initialize(0);
    fe_->evaluateValue(valV, v_coeff);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*WvalV,*weight_,*valV);
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *WvalV,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP,
                                                  false);
  }
  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
            const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
            const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    // Build force/control term
    ROL::Ptr<Intrepid::FieldContainer<Real>> F
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*F, *weight_, *fe_->N());
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *F,
                                                  *fe_->NdetJ(),
                                                  Intrepid::COMP_CPP,
                                                  false);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_22 is zero.");
  }

}; // QoI_Helmholtz


template <class Real>
class QoI_Helmholtz_ControlPenalty : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> fe_;

  Real innerAnnulusRadius_;
  Real outerAnnulusRadius_;
  int example_;

  ROL::Ptr<Intrepid::FieldContainer<Real>> weight_;

  bool insideDomain(const std::vector<Real> &x) const {
    bool val = true;
    if (example_==1) {
      const Real eps = std::sqrt(ROL::ROL_EPSILON<Real>());
      Real xnorm(0);
      const int d = x.size();
      for (int i = 0; i < d; ++i) {
        xnorm += x[i]*x[i];
      }
      xnorm = std::sqrt(xnorm);
      val = (xnorm <= outerAnnulusRadius_+eps && xnorm >= innerAnnulusRadius_-eps);
    }
    return val;
  }

  void computeDomainWeight(void) {
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);

    weight_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
   
    const Real one(1), zero(0);
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

public:
  QoI_Helmholtz_ControlPenalty(const ROL::Ptr<FE<Real>> &fe,
                               Teuchos::ParameterList &parlist)
  : fe_(fe), innerAnnulusRadius_(2.5), outerAnnulusRadius_(2.6) {
    example_ = parlist.sublist("Problem").get("Example",1);
    computeDomainWeight();
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Evaluate control penalty
    ROL::Ptr<Intrepid::FieldContainer<Real>> valZ, WvalZ;
    valZ  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    WvalZ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    valZ->initialize(); WvalZ->initialize();
    fe_->evaluateValue(valZ, z_coeff);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*WvalZ,*weight_,*valZ);
    fe_->computeIntegral(val,WvalZ,valZ,true);
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz_ControlCost::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    // Evaluate control penalty
    ROL::Ptr<Intrepid::FieldContainer<Real>> valZ, WvalZ;
    valZ  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    WvalZ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    valZ->initialize(); WvalZ->initialize();
    fe_->evaluateValue(valZ, z_coeff);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*WvalZ,*weight_,*valZ);
    grad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  *WvalZ,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP,
                                                  false);
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz_ControlCost::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz_ControlCost::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz_ControlCost::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    // Evaluate control penalty
    ROL::Ptr<Intrepid::FieldContainer<Real>> valV, WvalV;
    valV  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    WvalV = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    valV->initialize(); WvalV->initialize();
    fe_->evaluateValue(valV, v_coeff);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*WvalV,*weight_,*valV);
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *WvalV,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP,
                                                  false);
  }
  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
            const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
            const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    // Build force/control term
    ROL::Ptr<Intrepid::FieldContainer<Real>> F
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*F, *weight_, *fe_->N());
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *F,
                                                  *fe_->NdetJ(),
                                                  Intrepid::COMP_CPP,
                                                  false);
  }

}; // QoI_Helmholtz_ControlPenalty

#endif
