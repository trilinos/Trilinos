// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  qoi.hpp
    \brief Provides the interface for local (cell-based) quantities of interest.
*/

#ifndef ROL_PDEOPT_QOIK_HPP
#define ROL_PDEOPT_QOIK_HPP

#include "pdeK.hpp"

template <class Real, class DeviceType>
class QoI {
public:
  
  using scalar_view = Kokkos::DynRankView<Real,DeviceType>;

  virtual ~QoI() {}

  virtual Real value(scalar_view & val,
                     const scalar_view u_coeff,
                     const scalar_view z_coeff = scalar_view(),
                     const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) = 0;

  virtual void gradient_1(scalar_view & grad,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> gradient_1 not implemented.");
  }

  virtual void gradient_2(scalar_view & grad,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> gradient_2 not implemented.");
  }

  virtual std::vector<Real> gradient_3(std::vector<scalar_view > & grad,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> gradient_3 not implemented.");
  }

  virtual void HessVec_11(scalar_view & hess,
                          const scalar_view v_coeff,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> HessVec_11 not implemented.");
  }

  virtual void HessVec_12(scalar_view & hess,
                          const scalar_view v_coeff,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> HessVec_12 not implemented.");
  }

  virtual void HessVec_13(scalar_view & hess,
                          const ROL::Ptr<const std::vector<Real> > & v_param,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> HessVec_13 not implemented.");
  }

  virtual void HessVec_21(scalar_view & hess,
                          const scalar_view v_coeff,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> HessVec_21 not implemented.");
  }

  virtual void HessVec_22(scalar_view & hess,
                          const scalar_view v_coeff,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> HessVec_22 not implemented.");
  }

  virtual void HessVec_23(scalar_view & hess,
                          const ROL::Ptr<const std::vector<Real> > & v_param,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> HessVec_23 not implemented.");
  }

  virtual std::vector<Real> HessVec_31(std::vector<scalar_view > & hess,
                          const scalar_view v_coeff,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> HessVec_31 not implemented.");
  }

  virtual std::vector<Real> HessVec_32(std::vector<scalar_view > & hess,
                          const scalar_view v_coeff,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> HessVec_32 not implemented.");
  }

  virtual std::vector<Real> HessVec_33(std::vector<scalar_view > & hess,
                          const ROL::Ptr<const std::vector<Real> > & v_param,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> HessVec_33 not implemented.");
  }

  virtual void Hessian_11(scalar_view & hess,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_11 not implemented.");
  }

  virtual void Hessian_12(scalar_view & hess,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_12 not implemented.");
  }

  virtual void Hessian_13(std::vector<scalar_view > & hess,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_13 not implemented.");
  }

  virtual void Hessian_21(scalar_view & hess,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_21 not implemented.");
  }

  virtual void Hessian_22(scalar_view & hess,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_22 not implemented.");
  }

  virtual void Hessian_23(std::vector<scalar_view > & hess,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_23 not implemented.");
  }

  virtual void Hessian_31(std::vector<scalar_view > & hess,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_31 not implemented.");
  }

  virtual void Hessian_32(std::vector<scalar_view > & hess,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_32 not implemented.");
  }

  virtual void Hessian_33(std::vector<std::vector<Real> > & hess,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_33 not implemented.");
  }

private:
  std::vector<Real> param_;

protected:
  std::vector<Real> getParameter(void) const {
    return param_;
  }

public:
  void setParameter(const std::vector<Real> &param) {
    param_.assign(param.begin(),param.end());
  }

}; // QOI

#endif
