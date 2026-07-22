// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde.hpp
    \brief Provides the interface for local (cell-based) PDE residual computations.
*/

#ifndef PDEOPT_DYNPDE_HPP
#define PDEOPT_DYNPDE_HPP

#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_Basis.hpp"
#include "ROL_Ptr.hpp"
#include "ROL_TimeStamp.hpp"
#include "pde.hpp"

template <class Real>
class DynamicPDE {
public:
  virtual ~DynamicPDE() {}

  virtual void residual(ROL::Ptr<Intrepid::FieldContainer<Real>> & res,
                        const ROL::TimeStamp<Real> &ts,
                        const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                        const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                        const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                        const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) = 0;

  virtual void Jacobian_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                           const ROL::TimeStamp<Real> &ts,
                           const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                           const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                           const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                           const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Jacobian_uo not implemented.");
  }

  virtual void Jacobian_un(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                           const ROL::TimeStamp<Real> &ts,
                           const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                           const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                           const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                           const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Jacobian_un not implemented.");
  }

  virtual void Jacobian_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                           const ROL::TimeStamp<Real> &ts,
                           const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                           const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                           const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                           const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Jacobian_zf not implemented.");
  }

  virtual void Jacobian_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & jac,
                           const ROL::TimeStamp<Real> &ts,
                           const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                           const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                           const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                           const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Jacobian_3 not implemented.");
  }

  virtual void Hessian_uo_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_uo_uo not implemented.");
  }

  virtual void Hessian_uo_un(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_uo_un not implemented.");
  }

  virtual void Hessian_uo_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_uo_zf not implemented.");
  }

  virtual void Hessian_uo_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_uo_zp not implemented.");
  }

  virtual void Hessian_un_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_un_uo not implemented.");
  }

  virtual void Hessian_un_un(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_un_un not implemented.");
  }

  virtual void Hessian_un_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_un_zf not implemented.");
  }

  virtual void Hessian_un_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_un_zp not implemented.");
  }

  virtual void Hessian_zf_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_zf_uo not implemented.");
  }

  virtual void Hessian_zf_un(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_zf_un not implemented.");
  }

  virtual void Hessian_zf_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_zf_zf not implemented.");
  }

  virtual void Hessian_zf_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_zf_zp not implemented.");
  }

  virtual void Hessian_zp_uo(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_zp_uo not implemented.");
  }

  virtual void Hessian_zp_un(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_zp_un not implemented.");
  }

  virtual void Hessian_zp_zf(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_zp_zf not implemented.");
  }

  virtual void Hessian_zp_zp(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> & hess,
                             const ROL::TimeStamp<Real> &ts,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_zp_zp not implemented.");
  }

  virtual void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real>> &riesz) {
    throw Exception::NotImplemented(">>> RieszMap_1 not implemented.");
  }

  virtual void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real>> &riesz) {
    throw Exception::NotImplemented(">>> RieszMap_2 not implemented.");
  }

  virtual std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields() = 0;
  virtual std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields2() {
    return getFields();
  }

  virtual void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real>> &cellNodes,
                            const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> &bdryCellNodes,
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
