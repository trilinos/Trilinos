// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde.hpp
    \brief Implements the local PDE interface for the Poisson control problem.
*/

#ifndef PDE_ADV_DIFF_HPP
#define PDE_ADV_DIFF_HPP

#include "../../TOOLS/pde.hpp"
#include "../../TOOLS/fe.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real>
class PDE_adv_diff : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>> basisPtr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real>> cellCub_;
  // Finite element definition
  ROL::Ptr<FE<Real>> fe_vol_;

  bool isLTI_;
  Real T_;

public:
  PDE_adv_diff(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("Problem").get("Basis Order",1);
    if (basisOrder == 1) {
      basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
    }
    else if (basisOrder == 2) {
      basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real>>>();
    }
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();    // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                   // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree",2); // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                   // create default cubature
    // Time-dependent coefficients
    isLTI_ = !parlist.sublist("Problem").get("Time Varying Coefficients",false);
    T_     = parlist.sublist("Time Discretization").get("End Time",1.0);
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real>> & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    // INITIALIZE RESIDUAL
    res = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // COMPUTE PDE COEFFICIENTS
    ROL::Ptr<Intrepid::FieldContainer<Real>> kappa
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real>> V
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    ROL::Ptr<Intrepid::FieldContainer<Real>> rhs
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    computeCoefficients(kappa,V,rhs);
    // COMPUTE DIFFUSION TERM
    // Compute grad(U)
    ROL::Ptr<Intrepid::FieldContainer<Real>> gradU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    // Multiply kappa * grad(U)
    Intrepid::FieldContainer<Real> kappa_gradU(c, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(kappa_gradU,
                                                               *kappa,
                                                               *gradU_eval);
    // Integrate (kappa * grad(U)) . grad(N)
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  kappa_gradU,
                                                  *(fe_vol_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // ADD ADVECTION TERM TO RESIDUAL
    // Multiply V . grad(U)
    Intrepid::FieldContainer<Real> V_gradU(c, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(V_gradU,
                                                            *V,
                                                            *gradU_eval);
    // Integrate (V . grad(U)) * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  V_gradU,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // ADD RHS TO RESIDUAL
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *rhs,
                                                  (*fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);

    // ADD CONTROL TERM TO RESIDUAL
    int size = z_param->size();
    ROL::Ptr<Intrepid::FieldContainer<Real>> ctrl
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < size; ++i) {
      computeControlOperator(ctrl,(*z_param)[i],i);
      Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                    *ctrl,
                                                    *(fe_vol_->NdetJ()),
                                                    Intrepid::COMP_CPP, true);
    }
  }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    // INITILAIZE JACOBIAN
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // COMPUTE PDE COEFFICIENTS
    ROL::Ptr<Intrepid::FieldContainer<Real>> kappa
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real>> V
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    ROL::Ptr<Intrepid::FieldContainer<Real>> rhs
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    computeCoefficients(kappa,V,rhs);
    // COMPUTE DIFFUSION TERM
    // Multiply kappa * grad(N)
    Intrepid::FieldContainer<Real> kappa_gradN(c, f, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(kappa_gradN,
                                                                *kappa,
                                                                *(fe_vol_->gradN()));
    // Integrate (kappa * grad(N)) . grad(N)
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  kappa_gradN,
                                                  *(fe_vol_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // ADD ADVECTION TERM TO JACOBIAN
    // Multiply V . grad(N)
    Intrepid::FieldContainer<Real> V_gradN(c, f, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(V_gradN,
                                                             *V,
                                                             *(fe_vol_->gradN()));
    // Integrate (V . grad(U)) * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *(fe_vol_->NdetJ()),
                                                  V_gradN,
                                                  Intrepid::COMP_CPP, true);
  }

  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Jacobian_2): Jacobian is zero.");
  }

  void Jacobian_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    // ADD CONTROL TERM TO RESIDUAL
    int size = z_param->size();
    ROL::Ptr<Intrepid::FieldContainer<Real>> ctrl
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < size; ++i) {
      jac[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      computeControlOperator(ctrl,static_cast<Real>(1),i);
      Intrepid::FunctionSpaceTools::integrate<Real>(*(jac[i]),
                                                    *ctrl,
                                                    *(fe_vol_->NdetJ()),
                                                    Intrepid::COMP_CPP, false);
    }
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_12): Hessian is zero.");
  }

  void Hessian_13(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_13): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_22): Hessian is zero.");
  }

  void Hessian_23(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_23): Hessian is zero.");
  }

  void Hessian_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_31): Hessian is zero.");
  }

  void Hessian_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_32): Hessian is zero.");
  }

  void Hessian_33(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_33): Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    // GET DIMENSIONS
    int c = fe_vol_->N()->dimension(0);
    int f = fe_vol_->N()->dimension(1);
    // INITIALIZE RIESZ
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    *riesz = *fe_vol_->stiffMat();
    Intrepid::RealSpaceTools<Real>::add(*riesz,*(fe_vol_->massMat()));
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    // GET DIMENSIONS
    int c = fe_vol_->N()->dimension(0);
    int f = fe_vol_->N()->dimension(1);
    // INITIALIZE RIESZ
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    *riesz = *fe_vol_->massMat();
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields() {
    return basisPtrs_;
  }

  void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real>> &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds) {
    // Finite element definition.
    fe_vol_ = ROL::makePtr<FE<Real>>(volCellNodes,basisPtr_,cellCub_);
  }

  const ROL::Ptr<FE<Real>> getFE(void) const {
    return fe_vol_;
  }

private:
  Real evaluateDiffusivity(const std::vector<Real> &x) const {
    return static_cast<Real>(1e-1);
  }

  void evaluateVelocity(std::vector<Real> &adv, const std::vector<Real> &x) const {
    const Real half(0.5), five(5), ten(10);
    const Real a = five*half;
    const Real b = five + (ten-five)*half;
    adv[0] = b - a*x[0];
    adv[1] =     a*x[1];
    // Get time
    if (!isLTI_) {
      const Real one(1), two(2), pi2(2.0*M_PI), ten(10), c(5e-2);
      const Real t = PDE<Real>::getTime()/T_;
      // U0(t,x,y) =  cos((x-t)*2*pi).*sin((y-t)*2*pi).*exp(-2*t)/(2*pi) + 5e-2*(7.5 - 2.5*x).*exp(20*(t-1));
      // V0(t,x,y) = -sin((x-t)*2*pi).*cos((y-t)*2*pi).*exp(-2*t)/(2*pi) + 5e-2*(      2.5*y).*exp(20*(t-1));
      adv[0] *= c*std::exp(two*ten*(t-one));
      adv[1] *= c*std::exp(two*ten*(t-one));
      adv[0] += std::cos(pi2*(x[0]-t))*std::sin(pi2*(x[1]-t))*std::exp(-two*t)/pi2;
      adv[1] -= std::sin(pi2*(x[0]-t))*std::cos(pi2*(x[1]-t))*std::exp(-two*t)/pi2;
    }
  }

  Real evaluateRHS(const std::vector<Real> &x) const {
    const int ns = 5;             
    const Real zero(0), half(0.5), one(1), two(2);
    Real source(0), arg1(0), arg2(0), mag(0), x0(0), y0(0), sx(0), sy(0);
    // Upper and lower bounds on source magintudes
    const std::vector<Real> ml = {1.5, 1.2, 1.5, 1.2, 1.1};
    const std::vector<Real> mu = {2.5, 1.8, 1.9, 2.6, 1.5};
    // Period and phase of sources
    const std::vector<Real> st = {0.0, 0.2, 0.4, 0.6, 0.8};
    const std::vector<Real> dr = {1.0, 0.2, 0.5, 0.1, 0.2};
    // Upper and lower bounds on source locations
    const std::vector<Real> xl = {0.45, 0.75, 0.40, 0.05, 0.85};
    const std::vector<Real> xu = {0.55, 0.85, 0.60, 0.35, 0.95};
    const std::vector<Real> yl = {0.25, 0.55, 0.50, 0.45, 0.45};
    const std::vector<Real> yu = {0.35, 0.65, 0.70, 0.65, 0.55};
    // Upper and lower bounds on source widths
    const std::vector<Real> sxl = {0.03, 0.02, 0.01, 0.02, 0.015};
    const std::vector<Real> sxu = {0.07, 0.04, 0.05, 0.04, 0.025};
    const std::vector<Real> syl = {0.04, 0.01, 0.02, 0.02, 0.01};
    const std::vector<Real> syu = {0.12, 0.05, 0.04, 0.04, 0.03};
    // Get time
    Real time(0);
    if (!isLTI_) {
      time = PDE<Real>::getTime();
    }
    for (int i=0; i<ns; ++i) {
      mag = half*((mu[i]+ml[i]) + (mu[i]-ml[i]));
      if (!isLTI_) {
        mag *= (time/T_ > st[i] && time/T_ < st[i]+dr[i]) ? one : zero;
      }
      x0   = half*(( xl[i]+ xu[i]) + ( xu[i]- xl[i]));
      y0   = half*(( yl[i]+ yu[i]) + ( yu[i]- yl[i]));
      sx   = half*((sxl[i]+sxu[i]) + (sxu[i]-sxl[i]));
      sy   = half*((syl[i]+syu[i]) + (syu[i]-syl[i]));
      arg1 = std::pow((x[0]-x0)/sx, two);
      arg2 = std::pow((x[1]-y0)/sy, two);
      source += mag*std::exp(-half*(arg1+arg2));
    }
    return source;
  }

  void computeCoefficients(ROL::Ptr<Intrepid::FieldContainer<Real>> &kappa,
                           ROL::Ptr<Intrepid::FieldContainer<Real>> &V,
                           ROL::Ptr<Intrepid::FieldContainer<Real>> &rhs) const {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    std::vector<Real> pt(d), adv(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*fe_vol_->cubPts())(i,j,k);
        }
        // Compute diffusivity kappa
        (*kappa)(i,j) = evaluateDiffusivity(pt);
        // Compute advection velocity field V
        evaluateVelocity(adv,pt);
        for (int k = 0; k < d; ++k) {
          (*V)(i,j,k) = adv[k];
        }
        // Compute forcing term f
        (*rhs)(i,j) = -evaluateRHS(pt);
      }
    }
  }

  Real evaluateControlOperator(const std::vector<Real> &x, const int i) const {
    const Real sx(0.05), sy(0.05), half(0.5);
    const std::vector<Real> xl = {0.25, 0.50, 0.75, 0.25, 0.50, 0.75, 0.25, 0.50, 0.75};
    const std::vector<Real> yl = {0.25, 0.25, 0.25, 0.50, 0.50, 0.50, 0.75, 0.75, 0.75};
    return -std::exp(- half*(x[0]-xl[i])*(x[0]-xl[i]) / (sx*sx)
                     - half*(x[1]-yl[i])*(x[1]-yl[i]) / (sy*sy));
  }
  
  void computeControlOperator(ROL::Ptr<Intrepid::FieldContainer<Real>> &ctrl,
                              const Real z, const int I) const {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    std::vector<Real> pt(d), adv(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*fe_vol_->cubPts())(i,j,k);
        }
        // Compute control operator
        (*ctrl)(i,j) = -z*evaluateControlOperator(pt,I);
      }
    }
  }

}; // PDE_adv_diff

#endif
