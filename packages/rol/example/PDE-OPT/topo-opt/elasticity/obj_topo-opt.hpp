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

#ifndef PDEOPT_QOI_TOPOOPT_HPP
#define PDEOPT_QOI_TOPOOPT_HPP

#include "../../TOOLS/qoi.hpp"
#include "pde_topo-opt.hpp"

template <class Real>
class QoI_TopoOpt : public QoI<Real> {
private:
  const Teuchos::RCP<FE<Real> > fe_;
  const Teuchos::RCP<Load<Real> > load_;
  const Teuchos::RCP<FieldHelper<Real> > fieldHelper_;
  const Real scale_;

public:
  QoI_TopoOpt(const Teuchos::RCP<FE<Real> > &fe,
              const Teuchos::RCP<Load<Real> > &load,
              const Teuchos::RCP<FieldHelper<Real> > &fieldHelper,
              const Real scale = 1.0)
    : fe_(fe), load_(load), fieldHelper_(fieldHelper), scale_(scale) {}

  Real value(Teuchos::RCP<Intrepid::FieldContainer<Real> > & val,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
             const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    // Initialize output val
    val = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c));
    // Get components of the control
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > valU_eval(d);
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > load_eval(d);
    for (int i=0; i<d; ++i) {
      valU_eval[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      load_eval[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    }
    for (int i=0; i<d; ++i) {
      fe_->evaluateValue(valU_eval[i], U[i]);
    }
    load_->compute(load_eval, fe_, QoI<Real>::getParameter());
    for (int i=0; i<d; ++i) {
      fe_->computeIntegral(val,load_eval[i],valU_eval[i],true);
    }
    Intrepid::RealSpaceTools<Real>::scale(*val, scale_);
    return static_cast<Real>(0);
  }

  void gradient_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > G(d);
    for (int i=0; i<d; ++i) {
      G[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
    }
    // Evaluate on FE basis
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > load_eval(d);
    for (int i=0; i<d; ++i) {
      load_eval[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    }
    load_->compute(load_eval, fe_, QoI<Real>::getParameter());
    // Build local gradient of state tracking term
    for (int i=0; i<d; ++i) {
      Intrepid::FunctionSpaceTools::integrate<Real>(*G[i],
                                                    *load_eval[i],
                                                    *(fe_->NdetJ()),
                                                    Intrepid::COMP_CPP,
                                                    false);
      Intrepid::RealSpaceTools<Real>::scale(*G[i], scale_);
    }

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_TopoOpt::gradient_2 is zero.");
  }

  void HessVec_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_TopoOpt::HessVec_12 is zero.");
  }

  void HessVec_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_TopoOpt::HessVec_12 is zero.");
  }

  void HessVec_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_TopoOpt::HessVec_21 is zero.");
  }

  void HessVec_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_TopoOpt::HessVec_22 is zero.");
  }

}; // QoI_TopoOpt


template <class Real>
class QoI_Volume_TopoOpt : public QoI<Real> {
private:
  const Teuchos::RCP<FE<Real> > fe_;
  const Teuchos::RCP<FieldHelper<Real> > fieldHelper_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > ones_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > volFrac_;

public:
  QoI_Volume_TopoOpt(const Teuchos::RCP<FE<Real> > &fe,
                     const Teuchos::RCP<FieldHelper<Real> > &fieldHelper,
                     Teuchos::ParameterList &parlist)
  : fe_(fe), fieldHelper_(fieldHelper) {
    Real v0 = parlist.sublist("Problem").get("Volume Fraction",0.4);
    // Get relevant dimensions
    int c = fe_->cubPts()->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    ones_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,p));
    volFrac_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,p));
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*ones_)(i,j) = static_cast<Real>(1);
        (*volFrac_)(i,j) = v0;
      }
    }
  }

  Real value(Teuchos::RCP<Intrepid::FieldContainer<Real> > & val,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
             const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);

    // Initialize output val
    val = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c));

    // Get components of the control
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Evaluate on FE basis
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_->evaluateValue(valZ_eval, Z[0]);
    Intrepid::RealSpaceTools<Real>::subtract(*valZ_eval,*volFrac_);

    // Compute energy
    fe_->computeIntegral(val,valZ_eval,ones_);
    return static_cast<Real>(0);
  }

  void gradient_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_Volume_TopoOpt::gradient_1 is zero.");
  }

  void gradient_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int d = fe_->gradN()->dimension(3);

    // Initialize output grad
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > G(d);
    for (int i=0; i<d; ++i) {
      G[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
    }

    // Compute gradient of energy
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[0],
                                                  *ones_,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void HessVec_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_Volume_TopoOpt::HessVec_11 is zero.");
  }

  void HessVec_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_Volume_TopoOpt::HessVec_12 is zero.");
  }

  void HessVec_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_Volume_TopoOpt::HessVec_21 is zero.");
  }

  void HessVec_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_Volume_TopoOpt::HessVec_22 is zero.");
  }

}; // QoI_Volume


template <class Real>
class StdObjective_TopoOpt : public ROL::StdObjective<Real> {
public:
  StdObjective_TopoOpt() {}

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

}; // StdObjective_TopOpt

#endif
