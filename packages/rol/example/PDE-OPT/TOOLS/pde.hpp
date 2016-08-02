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
#include "Teuchos_RCP.hpp"

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

  virtual void residual(Teuchos::RCP<Intrepid::FieldContainer<Real> > & res,
                        const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                        const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff) = 0;

  virtual void Jacobian_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & jac,
                          const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff) = 0;

  virtual void Jacobian_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & jac,
                          const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff) = 0;

  virtual void Hessian_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                          const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff,
                          const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff) {
    throw Exception::NotImplemented(">>> Hessian_11 not implemented.");
  }

  virtual void Hessian_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                          const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff,
                          const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff) {
    throw Exception::NotImplemented(">>> Hessian_12 not implemented.");
  }

  virtual void Hessian_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                          const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff,
                          const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff) {
    throw Exception::NotImplemented(">>> Hessian_21 not implemented.");
  }

  virtual void Hessian_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                          const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff,
                          const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff) {
    throw Exception::NotImplemented(">>> Hessian_22 not implemented.");
  }

  virtual std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > getFields() = 0;

  virtual void setCellNodes(const Teuchos::RCP<Intrepid::FieldContainer<Real> > &cellNodes,
                            const std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > &bdryCellNodes,
                            const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds) = 0;

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

}; // PDE

#endif
