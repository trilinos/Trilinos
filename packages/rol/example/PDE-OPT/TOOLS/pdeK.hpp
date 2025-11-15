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

#ifndef PDEOPT_PDEK_HPP
#define PDEOPT_PDEK_HPP

#include "Intrepid2_Basis.hpp"
#include "ROL_Ptr.hpp"

namespace Exception {

  class NotImplemented : public Teuchos::ExceptionBase {
    public:
      NotImplemented(const std::string & what_arg) : Teuchos::ExceptionBase(what_arg) {}
  }; // NotImplemented

  class Zero : public Teuchos::ExceptionBase {
    public:
      Zero(const std::string & what_arg) : Teuchos::ExceptionBase(what_arg) {}
  }; // Zero

} // Exception

template <class Real, class DeviceType>
class PDE {
public:
  using scalar_view = Kokkos::DynRankView<Real,DeviceType>;
  using basis_ptr = Intrepid2::BasisPtr<DeviceType, Real, Real>;

  virtual ~PDE() {}

  virtual void residual(scalar_view & res,
                        const scalar_view u_coeff,
                        const scalar_view z_coeff = scalar_view(),
                        const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) = 0;

  virtual void Jacobian_1(scalar_view & jac,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Jacobian_1 not implemented.");
  }

  virtual void Jacobian_2(scalar_view & jac,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Jacobian_2 not implemented.");
  }

  virtual void Jacobian_3(std::vector<scalar_view > & jac,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Jacobian_3 not implemented.");
  }

  virtual void applyJacobian_1(scalar_view &jv,
                               const scalar_view v_coeff,
                               const scalar_view u_coeff,
                               const scalar_view z_coeff,
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> applyJacobian_1 not implemented.");
  }

  virtual void applyAdjointJacobian_1(scalar_view &jv,
                                      const scalar_view v_coeff,
                                      const scalar_view u_coeff,
                                      const scalar_view z_coeff,
                                      const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> applyAdjointJacobian_1 not implemented.");
  }

  virtual void applyJacobian_2(scalar_view &jv,
                               const scalar_view v_coeff,
                               const scalar_view u_coeff,
                               const scalar_view z_coeff,
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> applyJacobian_2 not implemented.");
  }

  virtual void applyAdjointJacobian_2(scalar_view &jv,
                                      const scalar_view v_coeff,
                                      const scalar_view u_coeff,
                                      const scalar_view z_coeff,
                                      const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> applyAdjointJacobian_2 not implemented.");
  }

  virtual void Hessian_11(scalar_view & hess,
                          const scalar_view l_coeff,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_11 not implemented.");
  }

  virtual void Hessian_12(scalar_view & hess,
                          const scalar_view l_coeff,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_12 not implemented.");
  }

  virtual void Hessian_13(std::vector<scalar_view > & hess,
                          const scalar_view l_coeff,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_13 not implemented.");
  }

  virtual void Hessian_21(scalar_view & hess,
                          const scalar_view l_coeff,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_21 not implemented.");
  }

  virtual void Hessian_22(scalar_view & hess,
                          const scalar_view l_coeff,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_22 not implemented.");
  }

  virtual void Hessian_23(std::vector<scalar_view > & hess,
                          const scalar_view l_coeff,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_23 not implemented.");
  }

  virtual void Hessian_31(std::vector<scalar_view > & hess,
                          const scalar_view l_coeff,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_31 not implemented.");
  }

  virtual void Hessian_32(std::vector<scalar_view > & hess,
                          const scalar_view l_coeff,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_32 not implemented.");
  }

  virtual void Hessian_33(std::vector<std::vector<scalar_view > > & hess,
                          const scalar_view l_coeff,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_33 not implemented.");
  }

  virtual void applyHessian_11(scalar_view & hess,
                               const scalar_view v_coeff,
                               const scalar_view l_coeff,
                               const scalar_view u_coeff,
                               const scalar_view z_coeff = scalar_view(),
                               const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> applyHessian_11 not implemented.");
  }

  virtual void applyHessian_12(scalar_view & hess,
                               const scalar_view v_coeff,
                               const scalar_view l_coeff,
                               const scalar_view u_coeff,
                               const scalar_view z_coeff = scalar_view(),
                               const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> applyHessian_12 not implemented.");
  }

  virtual void applyHessian_21(scalar_view & hess,
                               const scalar_view v_coeff,
                               const scalar_view l_coeff,
                               const scalar_view u_coeff,
                               const scalar_view z_coeff = scalar_view(),
                               const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> applyHessian_21 not implemented.");
  }

  virtual void applyHessian_22(scalar_view & hess,
                               const scalar_view v_coeff,
                               const scalar_view l_coeff,
                               const scalar_view u_coeff,
                               const scalar_view z_coeff = scalar_view(),
                               const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> applyHessian_22 not implemented.");
  }

  virtual void RieszMap_1(scalar_view &riesz) {
    throw Exception::NotImplemented(">>> RieszMap_1 not implemented.");
  }

  virtual void RieszMap_2(scalar_view &riesz) {
    throw Exception::NotImplemented(">>> RieszMap_2 not implemented.");
  }

  virtual std::vector<basis_ptr> getFields() = 0;
  virtual std::vector<basis_ptr> getFields2() {
    return getFields();
  }

  virtual void setCellNodes(const scalar_view &cellNodes,
                            const std::vector<std::vector<scalar_view>> &bdryCellNodes,
                            const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds) = 0;

  virtual void setFieldPattern(const std::vector<std::vector<int>> &fieldPattern) {}
  virtual void setFieldPattern(const std::vector<std::vector<int>> &fieldPattern1,
                               const std::vector<std::vector<int>> &fieldPattern2) {
    setFieldPattern(fieldPattern1);
  }

  virtual void printData(std::string tag,
                         const scalar_view u_coeff,
                         const scalar_view z_coeff = scalar_view(),
                         const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {}

  virtual void printCellAverages(std::string tag,
                                 const scalar_view u_coeff,
                                 const scalar_view z_coeff = scalar_view(),
                                 const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {}

private:
  Real time_;
  std::vector<Real> param_;

protected:
  Real getTime(void) const {
    return time_;
  }

  std::vector<Real> getParameter(void) const {
    return param_;
  }

public:
  void setTime(const Real time) {
    time_ = time;
  }

  void setParameter(const std::vector<Real> &param) {
    param_.assign(param.begin(),param.end());
  }

}; // PDE

#endif
