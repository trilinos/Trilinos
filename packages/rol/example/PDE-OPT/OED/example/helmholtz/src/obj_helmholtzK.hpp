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

#ifndef PDEOPT_QOI_HELMHOLTZ_OCT_REALK_HPP
#define PDEOPT_QOI_HELMHOLTZ_OCT_REALK_HPP

#include "../../../../TOOLS/qoiK.hpp"
#include "ROL_StdObjective.hpp"
#include "pde_helmholtzK.hpp"

template <class Real,class DeviceType>
class QoI_Helmholtz_Observation : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> fe_;
  Real width_, coeff_;

  Real observationFunc(const std::vector<Real> &x, const std::vector<Real> &loc, int deriv = 0, int dim1 = 0, int dim2 = 0) const {
    const Real zero(0), half(0.5), one(1);
    const int d = x.size();
    Real dot(0), width2 = std::pow(width_,2);
    for (int i = 0; i < d; ++i)
      dot += std::pow((x[i]-loc[i]),2)/width2;
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

  void evaluateObservation(scalar_view &out, const std::vector<Real> &loc, int deriv = 0, int dim1 = 0, int dim2 = 0) const {
    const int c = fe_->gradN().extent_int(0);
    const int p = fe_->gradN().extent_int(2);
    const int d = fe_->gradN().extent_int(3);
    out = scalar_view("out",c,p);
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k)
          x[k] = (fe_->cubPts())(i,j,k);
        out(i,j) = observationFunc(x,loc,deriv,dim1,dim2);
      }
    }
  }

public:
  QoI_Helmholtz_Observation(const ROL::Ptr<fe_type> &fe,
                            ROL::ParameterList &parlist)
    : fe_(fe) {
    width_ = parlist.sublist("Problem").get("Microphone Width",5e-2);
    coeff_ = static_cast<Real>(1)/(static_cast<Real>(2.0*M_PI)*std::pow(width_,2));
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
    // Get components of the state
    scalar_view u_real("u_real", c, p);
    fe_->evaluateValue(u_real, u_coeff);
    // Compute phi
    scalar_view phi;
    evaluateObservation(phi,QoI<Real,DeviceType>::getParameter());
    // Integrate observation
    fe_->computeIntegral(val,phi,u_real,false);
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);
    // Initialize output grad
    grad = scalar_view("grad", c, f);
    // Compute phi
    scalar_view phi;
    evaluateObservation(phi,QoI<Real,DeviceType>::getParameter());
    // Integrate observation derivative
    fst::integrate(grad,phi,fe_->NdetJ(),false);
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz_Observation::gradient_2 is zero.");
  }

  std::vector<Real> gradient_3(std::vector<scalar_view> & grad,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz_Observation::gradient_3 is zero.");
    // Get relevant dimensions
    const int c = fe_->gradN().extent_int(0);
    const int p = fe_->gradN().extent_int(2);
    const int d = fe_->gradN().extent_int(3);
    // Get components of the state
    scalar_view u_real("u_real", c, p);
    fe_->evaluateValue(u_real, u_coeff);
    for (int i = 0; i < d; ++i) {
      // Initialize output val
      grad[i] = scalar_view("grad", c);
      // Compute phi'
      scalar_view phi;
      evaluateObservation(phi,*z_param,1,i);
      // Integrate observation gradient
      fe_->computeIntegral(grad[i],phi,u_real,false);
    }
    std::vector<Real> empty(d);
    return empty;
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_11 is zero.");
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_12 is zero.");
  }

  void HessVec_13(scalar_view & hess,
                          const ROL::Ptr<const std::vector<Real>> & v_param,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_13 is zero.");
    // Get relevant dimensions
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);
    const int p = fe_->gradN().extent_int(2);
    const int d = fe_->gradN().extent_int(3);
    // Initialize hessian
    hess = scalar_view("hess", c, f);
    // Get components of the state
    scalar_view vphi("vphi", c, p);
    for (int i = 0; i < d; ++i) {
      // Compute phi'
      scalar_view phi;
      evaluateObservation(phi,*z_param,1,i);
      // Weight phi'
      rst::scale(vphi, phi, (*v_param)[i]);
      // Integrate observation gradient
      fst::integrate(hess,vphi,fe_->NdetJ(),true);
    }
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_22 is zero.");
  }

  void HessVec_23(scalar_view & hess,
                          const ROL::Ptr<const std::vector<Real>> & v_param,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QOI_Helmholtz::HessVec_23 is zero.");
  }

  std::vector<Real> HessVec_31(std::vector<scalar_view> & hess,
                          const scalar_view v_coeff,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_31 is zero.");
    // Get relevant dimensions
    const int c = fe_->gradN().extent_int(0);
    const int p = fe_->gradN().extent_int(2);
    const int d = fe_->gradN().extent_int(3);
    // Get components of the state
    scalar_view v_real("v_real", c, p);
    fe_->evaluateValue(v_real, v_coeff);
    for (int i = 0; i < d; ++i) {
      // Initialize output hess
      hess[i] = scalar_view("hess", c);
      // Compute phi
      scalar_view phi;
      evaluateObservation(phi,*z_param,1,i);
      // Integrate observation
      fe_->computeIntegral(hess[i],phi,v_real,false);
    }
    std::vector<Real> empty(d);
    return empty;
  }

  std::vector<Real> HessVec_32(std::vector<scalar_view> & hess,
                          const scalar_view v_coeff,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_32 is zero.");
  }

  std::vector<Real> HessVec_33(std::vector<scalar_view> & hess,
                          const ROL::Ptr<const std::vector<Real>> & v_param,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_33 is zero.");
    // Get relevant dimensions
    const int c = fe_->gradN().extent_int(0);
    const int p = fe_->gradN().extent_int(2);
    const int d = fe_->gradN().extent_int(3);
    // Get components of the state
    scalar_view u_real("u_real", c, p);
    fe_->evaluateValue(u_real, u_coeff);

    scalar_view vphi("vphi", c, p);
    for (int i = 0; i < d; ++i) {
      // Initialize output val
      hess[i] = scalar_view("hess", c);
      // Compute phi'
      Kokkos::deep_copy(vphi,static_cast<Real>(0));
      for (int j = 0; j < d; ++j) {
        scalar_view phi;
        evaluateObservation(phi,*z_param,2,i,j);
        rst::scale(phi, (*v_param)[j]);
        rst::add(vphi, phi);
      }
      // Integrate observation gradient
      fe_->computeIntegral(hess[i],vphi,u_real,false);
    }
    std::vector<Real> empty(d);
    return empty;
  }

}; // QoI_Helmholtz_Observation

#endif
