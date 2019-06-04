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

/*! \file  obj.hpp
    \brief Provides the interface for local (cell-based) objective function computations.
*/

#ifndef PDEOPT_QOI_ENERGY_POISSON_TOPOPT_HPP
#define PDEOPT_QOI_ENERGY_POISSON_TOPOPT_HPP

#include "../../TOOLS/qoi.hpp"
#include "pde_poisson_topOpt.hpp"

template <class Real>
class QoI_Energy_Poisson_TopOpt : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real> > fe_;
  const ROL::Ptr<Intrepid::FieldContainer<Real> > force_eval_;
  const Real scale_;

public:
  QoI_Energy_Poisson_TopOpt(const ROL::Ptr<FE<Real> > &fe,
                            const ROL::Ptr<Intrepid::FieldContainer<Real> > &force_eval,
                            const Real scale = static_cast<Real>(1))
    : fe_(fe), force_eval_(force_eval), scale_(scale) {}

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->cubPts()->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Evaluate state on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real> > valU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute energy
    fe_->computeIntegral(val,valU_eval,force_eval_);
    Intrepid::RealSpaceTools<Real>::scale(*val, scale_);
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->N()->dimension(0);
    int f = fe_->N()->dimension(1);
    // Initialize output grad
    grad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Compute gradient of energy
    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  *force_eval_,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*grad, scale_);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Tracking_Poisson::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Energy_Poisson_TopOpt::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Energy_Poisson_TopOpt::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Energy_Poisson_TopOpt::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Energy_Poisson_TopOpt::HessVec_22 is zero.");
  }

}; // QoI_Energy

template <class Real>
class QoI_Volume_Poisson_TopOpt : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real> > fe_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > ones_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > volFrac_;

public:
  QoI_Volume_Poisson_TopOpt(const ROL::Ptr<FE<Real> > &fe,
                            Teuchos::ParameterList &parlist)
    : fe_(fe) {
    Real v0 = parlist.sublist("Problem").get("Volume Fraction",0.4);
    // Get relevant dimensions
    int c = fe_->cubPts()->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    ones_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    volFrac_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*ones_)(i,j) = static_cast<Real>(1);
        (*volFrac_)(i,j) = v0;
      }
    }
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->cubPts()->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Evaluate state on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valZ_eval, z_coeff);
    Intrepid::RealSpaceTools<Real>::subtract(*valZ_eval,*volFrac_);
    // Compute energy
    fe_->computeIntegral(val,valZ_eval,ones_);
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Volume_Poisson::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->N()->dimension(0);
    int f = fe_->N()->dimension(1);
    // Initialize output grad
    grad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Compute gradient of energy
    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  *ones_,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Volume_Poisson_TopOpt::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Volume_Poisson_TopOpt::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Volume_Poisson_TopOpt::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Volume_Poisson_TopOpt::HessVec_22 is zero.");
  }

}; // QoI_Volume

template <class Real>
class StdObjective_Poisson_TopOpt : public ROL::StdObjective<Real> {
public:
  StdObjective_Poisson_TopOpt() {}

  Real value(const std::vector<Real> &x, Real &tol) {
    return x[0];
  }

  void gradient(std::vector<Real> &g, const std::vector<Real> &x, Real &tol) {
    const Real one(1);
    g[0] = one;
  }

  void hessVec(std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    const Real zero(0);
    hv[0] = zero;
  }

}; // OBJ_SCALAR

#endif
