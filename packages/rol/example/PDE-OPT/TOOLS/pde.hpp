// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
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
