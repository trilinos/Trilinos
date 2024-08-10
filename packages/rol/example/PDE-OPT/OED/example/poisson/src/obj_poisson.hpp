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

#ifndef PDEOPT_QOI_OED_POISSON_HPP
#define PDEOPT_QOI_OED_POISSON_HPP

#include "../../../../TOOLS/qoi.hpp"
#include "ROL_StdObjective.hpp"
#include "pde_poisson.hpp"

template <class Real>
class QoI_Poisson_Observation : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> fe_;
  Real width_, coeff_;

  Real observationFunc(const std::vector<Real> &x, const std::vector<Real> &loc, const int deriv = 0, const int dim1 = 0, const int dim2 = 0) const {
    const Real zero(0), half(0.5), one(1);
    const int d = x.size();
    Real dot(0), width2 = std::pow(width_,2);
    for (int i = 0; i < d; ++i) {
      dot += std::pow((x[i]-loc[i]),2)/width2;
    }
    Real val = std::exp(-half*dot);
    if (deriv == 1) {
      val *= (x[dim1]-loc[dim1])/width2;
    }
    else if (deriv == 2) {
      Real v0 = (dim1==dim2 ? -one/width2 : zero);
      Real v1 = (x[dim1]-loc[dim1])/width2;
      Real v2 = (x[dim2]-loc[dim2])/width2;
      val *= (v0 + v1*v2);
    }
    return coeff_ * val;
  }

  void evaluateObservation(ROL::Ptr<Intrepid::FieldContainer<Real>> &out, const std::vector<Real> &loc, const int deriv = 0, const int dim1 = 0, const int dim2 = 0) const {
    const int c = fe_->gradN()->dimension(0);
    const int p = fe_->gradN()->dimension(2);
    const int d = fe_->gradN()->dimension(3);
    out = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (*fe_->cubPts())(i,j,k);
        }
        (*out)(i,j) = observationFunc(x,loc,deriv,dim1,dim2);
      }
    }
  }

public:
  QoI_Poisson_Observation(const ROL::Ptr<FE<Real>> &fe,
                            Teuchos::ParameterList &parlist)
    : fe_(fe) {
    width_ = parlist.sublist("Problem").get("Microphone Width",5e-2);
    coeff_ = static_cast<Real>(1)/(static_cast<Real>(2.0*M_PI)*std::pow(width_,2));
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int p = fe_->gradN()->dimension(2);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Get components of the state
    ROL::Ptr<Intrepid::FieldContainer<Real>> u_real, phi;
    u_real = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(u_real, u_coeff);
    // Compute phi
    evaluateObservation(phi,QoI<Real>::getParameter());
    // Integrate observation
    fe_->computeIntegral(val,phi,u_real,false);
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    // Initialize output grad
    grad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Compute phi
    ROL::Ptr<Intrepid::FieldContainer<Real>> phi;
    evaluateObservation(phi,QoI<Real>::getParameter());
    // Integrate observation derivative
    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  *phi,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP,
                                                  false);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Poisson_Observation::gradient_2 is zero.");
  }

  std::vector<Real> gradient_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & grad,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Poisson_Observation::gradient_3 is zero.");
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int p = fe_->gradN()->dimension(2);
    const int d = fe_->gradN()->dimension(3);
    // Get components of the state
    ROL::Ptr<Intrepid::FieldContainer<Real>> u_real
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(u_real, u_coeff);
    for (int i = 0; i < d; ++i) {
      // Initialize output val
      grad[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
      // Compute phi'
      ROL::Ptr<Intrepid::FieldContainer<Real>> phi;
      evaluateObservation(phi,*z_param,1,i);
      // Integrate observation gradient
      fe_->computeIntegral(grad[i],phi,u_real,false);
    }
    std::vector<Real> empty(d);
    return empty;
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Poisson::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Poisson::HessVec_12 is zero.");
  }

  void HessVec_13(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                          const ROL::Ptr<const std::vector<Real>> & v_param,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Poisson::HessVec_13 is zero.");
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    const int d = fe_->gradN()->dimension(3);
    // Initialize hessian
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Get components of the state
    ROL::Ptr<Intrepid::FieldContainer<Real>> vphi
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < d; ++i) {
      // Compute phi'
      ROL::Ptr<Intrepid::FieldContainer<Real>> phi;
      evaluateObservation(phi,*z_param,1,i);
      // Weight phi'
      Intrepid::RealSpaceTools<Real>::scale(*vphi, *phi, (*v_param)[i]);
      // Integrate observation gradient
      Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                    *vphi,
                                                    *fe_->NdetJ(),
                                                    Intrepid::COMP_CPP,
                                                    true);
    }
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Poisson::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Poisson::HessVec_22 is zero.");
  }

  void HessVec_23(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                          const ROL::Ptr<const std::vector<Real>> & v_param,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QOI_Poisson::HessVec_23 is zero.");
  }

  std::vector<Real> HessVec_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Poisson::HessVec_31 is zero.");
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int p = fe_->gradN()->dimension(2);
    const int d = fe_->gradN()->dimension(3);
    // Get components of the state
    ROL::Ptr<Intrepid::FieldContainer<Real>> v_real
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(v_real, v_coeff);
    for (int i = 0; i < d; ++i) {
      // Initialize output hess
      hess[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
      // Compute phi
      ROL::Ptr<Intrepid::FieldContainer<Real>> phi;
      evaluateObservation(phi,*z_param,1,i);
      // Integrate observation
      fe_->computeIntegral(hess[i],phi,v_real,false);
    }
    std::vector<Real> empty(d);
    return empty;
  }

  std::vector<Real> HessVec_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Poisson::HessVec_32 is zero.");
  }

  std::vector<Real> HessVec_33(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                          const ROL::Ptr<const std::vector<Real>> & v_param,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Poisson::HessVec_33 is zero.");
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int p = fe_->gradN()->dimension(2);
    const int d = fe_->gradN()->dimension(3);
    // Get components of the state
    ROL::Ptr<Intrepid::FieldContainer<Real>> u_real
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(u_real, u_coeff);

    ROL::Ptr<Intrepid::FieldContainer<Real>> vphi
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < d; ++i) {
      // Initialize output val
      hess[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
      // Compute phi'
      vphi->initialize();
      for (int j = 0; j < d; ++j) {
        ROL::Ptr<Intrepid::FieldContainer<Real>> phi;
        evaluateObservation(phi,*z_param,2,i,j);
        Intrepid::RealSpaceTools<Real>::scale(*phi, (*v_param)[j]);
        Intrepid::RealSpaceTools<Real>::add(*vphi, *phi);
      }
      // Integrate observation gradient
      fe_->computeIntegral(hess[i],vphi,u_real,false);
    }
    std::vector<Real> empty(d);
    return empty;
  }

}; // QoI_Poisson_Observation

#endif
