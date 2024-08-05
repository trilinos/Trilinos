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

#ifndef PDEOPT_PDE_HPP
#define PDEOPT_PDE_HPP

#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_Basis.hpp"
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

template <class Real>
class PDE {
public:
  virtual ~PDE() {}

  virtual void residual(ROL::Ptr<Intrepid::FieldContainer<Real> > & res,
                        const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                        const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                        const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) = 0;

  virtual void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Jacobian_1 not implemented.");
  }

  virtual void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Jacobian_2 not implemented.");
  }

  virtual void Jacobian_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & jac,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Jacobian_3 not implemented.");
  }

  virtual void applyJacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real>> &jv,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff,
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> applyJacobian_1 not implemented.");
  }

  virtual void applyAdjointJacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real>> &jv,
                                      const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                                      const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                                      const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff,
                                      const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> applyAdjointJacobian_1 not implemented.");
  }

  virtual void applyJacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real>> &jv,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff,
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> applyJacobian_2 not implemented.");
  }

  virtual void applyAdjointJacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real>> &jv,
                                      const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                                      const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                                      const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff,
                                      const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> applyAdjointJacobian_2 not implemented.");
  }

  virtual void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_11 not implemented.");
  }

  virtual void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_12 not implemented.");
  }

  virtual void Hessian_13(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_13 not implemented.");
  }

  virtual void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_21 not implemented.");
  }

  virtual void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_22 not implemented.");
  }

  virtual void Hessian_23(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_23 not implemented.");
  }

  virtual void Hessian_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_31 not implemented.");
  }

  virtual void Hessian_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_32 not implemented.");
  }

  virtual void Hessian_33(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > & hess,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> Hessian_33 not implemented.");
  }

  virtual void applyHessian_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                               const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> applyHessian_11 not implemented.");
  }

  virtual void applyHessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                               const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> applyHessian_12 not implemented.");
  }

  virtual void applyHessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                               const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> applyHessian_21 not implemented.");
  }

  virtual void applyHessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                               const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> applyHessian_22 not implemented.");
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

  virtual void printData(std::string tag,
                         const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                         const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                         const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {}

  virtual void printCellAverages(std::string tag,
                                 const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                                 const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
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
