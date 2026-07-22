// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  dynpdeK.hpp
    \brief Provides the interface for local (cell-based) PDE residual computations.
*/

#ifndef PDEOPT_DYNPDEK_HPP
#define PDEOPT_DYNPDEK_HPP

#include "ROL_Ptr.hpp"
#include "ROL_TimeStamp.hpp"
#include "Intrepid2_Basis.hpp"

template <class Real, class DeviceType>
class DynamicPDE {
public:
  using scalar_view = Kokkos::DynRankView<Real,DeviceType>;
  using basis_ptr = Intrepid2::BasisPtr<DeviceType, Real, Real>;

  virtual ~DynamicPDE() {}

  virtual void residual(scalar_view & res,
                        const ROL::TimeStamp<Real> &ts,
                        const scalar_view uo_coeff,
                        const scalar_view un_coeff,
                        const scalar_view z_coeff = scalar_view(),
                        const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) = 0;

  virtual void Jacobian_uo(scalar_view & jac,
                           const ROL::TimeStamp<Real> &ts,
                           const scalar_view uo_coeff,
                           const scalar_view un_coeff,
                           const scalar_view z_coeff = scalar_view(),
                           const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Jacobian_uo not implemented.");
  }

  virtual void Jacobian_un(scalar_view & jac,
                           const ROL::TimeStamp<Real> &ts,
                           const scalar_view uo_coeff,
                           const scalar_view un_coeff,
                           const scalar_view z_coeff = scalar_view(),
                           const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Jacobian_un not implemented.");
  }

  virtual void Jacobian_zf(scalar_view & jac,
                           const ROL::TimeStamp<Real> &ts,
                           const scalar_view uo_coeff,
                           const scalar_view un_coeff,
                           const scalar_view z_coeff = scalar_view(),
                           const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Jacobian_zf not implemented.");
  }

  virtual void Jacobian_zp(std::vector<scalar_view> & jac,
                           const ROL::TimeStamp<Real> &ts,
                           const scalar_view uo_coeff,
                           const scalar_view un_coeff,
                           const scalar_view z_coeff = scalar_view(),
                           const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Jacobian_3 not implemented.");
  }

  virtual void Hessian_uo_uo(scalar_view & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const scalar_view l_coeff,
                             const scalar_view uo_coeff,
                             const scalar_view un_coeff,
                             const scalar_view z_coeff = scalar_view(),
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_uo_uo not implemented.");
  }

  virtual void Hessian_uo_un(scalar_view & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const scalar_view l_coeff,
                             const scalar_view uo_coeff,
                             const scalar_view un_coeff,
                             const scalar_view z_coeff = scalar_view(),
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_uo_un not implemented.");
  }

  virtual void Hessian_uo_zf(scalar_view & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const scalar_view l_coeff,
                             const scalar_view uo_coeff,
                             const scalar_view un_coeff,
                             const scalar_view z_coeff = scalar_view(),
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_uo_zf not implemented.");
  }

  virtual void Hessian_uo_zp(std::vector<scalar_view> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const scalar_view l_coeff,
                             const scalar_view uo_coeff,
                             const scalar_view un_coeff,
                             const scalar_view z_coeff = scalar_view(),
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_uo_zp not implemented.");
  }

  virtual void Hessian_un_uo(scalar_view & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const scalar_view l_coeff,
                             const scalar_view uo_coeff,
                             const scalar_view un_coeff,
                             const scalar_view z_coeff = scalar_view(),
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_un_uo not implemented.");
  }

  virtual void Hessian_un_un(scalar_view & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const scalar_view l_coeff,
                             const scalar_view uo_coeff,
                             const scalar_view un_coeff,
                             const scalar_view z_coeff = scalar_view(),
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_un_un not implemented.");
  }

  virtual void Hessian_un_zf(scalar_view & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const scalar_view l_coeff,
                             const scalar_view uo_coeff,
                             const scalar_view un_coeff,
                             const scalar_view z_coeff = scalar_view(),
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_un_zf not implemented.");
  }

  virtual void Hessian_un_zp(std::vector<scalar_view> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const scalar_view l_coeff,
                             const scalar_view uo_coeff,
                             const scalar_view un_coeff,
                             const scalar_view z_coeff = scalar_view(),
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_un_zp not implemented.");
  }

  virtual void Hessian_zf_uo(scalar_view & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const scalar_view l_coeff,
                             const scalar_view uo_coeff,
                             const scalar_view un_coeff,
                             const scalar_view z_coeff = scalar_view(),
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_zf_uo not implemented.");
  }

  virtual void Hessian_zf_un(scalar_view & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const scalar_view l_coeff,
                             const scalar_view uo_coeff,
                             const scalar_view un_coeff,
                             const scalar_view z_coeff = ROL::nullPtr,
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_zf_un not implemented.");
  }

  virtual void Hessian_zf_zf(scalar_view & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const scalar_view l_coeff,
                             const scalar_view uo_coeff,
                             const scalar_view un_coeff,
                             const scalar_view z_coeff = scalar_view(),
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_zf_zf not implemented.");
  }

  virtual void Hessian_zf_zp(std::vector<scalar_view> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const scalar_view l_coeff,
                             const scalar_view uo_coeff,
                             const scalar_view un_coeff,
                             const scalar_view z_coeff = scalar_view(),
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_zf_zp not implemented.");
  }

  virtual void Hessian_zp_uo(std::vector<scalar_view> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const scalar_view l_coeff,
                             const scalar_view uo_coeff,
                             const scalar_view un_coeff,
                             const scalar_view z_coeff = scalar_view(),
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_zp_uo not implemented.");
  }

  virtual void Hessian_zp_un(std::vector<scalar_view> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const scalar_view l_coeff,
                             const scalar_view uo_coeff,
                             const scalar_view un_coeff,
                             const scalar_view z_coeff = scalar_view(),
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_zp_un not implemented.");
  }

  virtual void Hessian_zp_zf(std::vector<scalar_view> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const scalar_view l_coeff,
                             const scalar_view uo_coeff,
                             const scalar_view un_coeff,
                             const scalar_view z_coeff = scalar_view(),
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_zp_zf not implemented.");
  }

  virtual void Hessian_zp_zp(std::vector<std::vector<scalar_view>> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const scalar_view l_coeff,
                             const scalar_view uo_coeff,
                             const scalar_view un_coeff,
                             const scalar_view z_coeff = scalar_view(),
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_zp_zp not implemented.");
  }

  virtual void RieszMap_1(scalar_view &riesz) {
    throw Exception::NotImplemented(">>> RieszMap_1 not implemented.");
  }

  virtual void RieszMap_2(scalar_view &riesz) {
    throw Exception::NotImplemented(">>> RieszMap_2 not implemented.");
  }

  virtual std::vector<basis_ptr> getFields() = 0;
  //virtual std::vector<basis_ptr> getFields2() {
  //  return getFields();
  //}

  virtual void setCellNodes(const scalar_view &cellNodes,
                            const std::vector<std::vector<scalar_view>> &bdryCellNodes,
                            const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds) = 0;

  virtual void setFieldPattern(const std::vector<std::vector<int>> &fieldPattern) {}
  virtual void setFieldPattern(const std::vector<std::vector<int>> &fieldPattern1,
                               const std::vector<std::vector<int>> &fieldPattern2) {
    setFieldPattern(fieldPattern1);
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

}; // DynamicPDE

#endif
