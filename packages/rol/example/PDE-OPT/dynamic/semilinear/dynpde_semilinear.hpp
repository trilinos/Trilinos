// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde.hpp
    \brief Implements the local PDE interface for the Poisson-Boltzmann control problem.
*/

#ifndef PDEOPT_DYNAMIC_SEMILINEAR_HPP
#define PDEOPT_DYNAMIC_SEMILINEAR_HPP

#include "../../TOOLS/dynpde.hpp"
#include "../../TOOLS/fe.hpp"
#include "pde_semilinear.hpp"

template <class Real>
class DynamicPDE_Semilinear : public DynamicPDE<Real> {
private:
  // Finite element definition
  ROL::Ptr<FE<Real>> fe_vol_;
  // Steady PDE without Dirichlet BC
  ROL::Ptr<PDE_Semilinear<Real>> pde_;
  Real theta_;

public:
  DynamicPDE_Semilinear(Teuchos::ParameterList &parlist) {
    pde_ = ROL::makePtr<PDE_Semilinear<Real>>(parlist);
    // Time-dependent coefficients
    theta_  = parlist.sublist("Time Discretization").get("Theta",1.0);
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real>> & res,
                const ROL::TimeStamp<Real> & ts,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real one(1);
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew-told;
    // INITIALIZE STORAGE
    ROL::Ptr<Intrepid::FieldContainer<Real>> pde
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    // COMPUTE OLD RESIDUAL
    pde_->setTime(told);
    pde_->residual(res,uo_coeff,z_coeff,z_param);
    // Integrate Uold * N
    fe_vol_->evaluateValue(valU_eval, uo_coeff);
    Intrepid::FunctionSpaceTools::integrate<Real>(*pde,
                                                  *valU_eval,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*res, (one-theta_)*dt);
    Intrepid::RealSpaceTools<Real>::subtract(*res, *pde);
    // COMPUTE NEW RESIDUAL
    pde_->setTime(tnew);
    pde_->residual(pde,un_coeff,z_coeff,z_param);
    // Integrate Uold * N
    valU_eval->initialize();
    fe_vol_->evaluateValue(valU_eval, un_coeff);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *valU_eval,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    Intrepid::RealSpaceTools<Real>::scale(*pde, theta_*dt);
    Intrepid::RealSpaceTools<Real>::add(*res, *pde);
  }

  void Jacobian_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                   const ROL::TimeStamp<Real> & ts,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                   const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real one(1);
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew-told;
    // INITILAIZE JACOBIAN
    ROL::Ptr<Intrepid::FieldContainer<Real>> pde
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // COMPUTE OLD RESIDUAL
    pde_->setTime(told);
    pde_->Jacobian_1(jac,uo_coeff,z_coeff,z_param);
    // Integrate N * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*pde,
                                                  *(fe_vol_->N()),
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*jac, (one-theta_)*dt);
    Intrepid::RealSpaceTools<Real>::subtract(*jac, *pde);
  }

  void Jacobian_un(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                   const ROL::TimeStamp<Real> & ts,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                   const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew-told;
    // INITILAIZE JACOBIAN
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    ROL::Ptr<Intrepid::FieldContainer<Real>> pde;
    // COMPUTE NEW RESIDUAL
    pde_->setTime(tnew);
    pde_->Jacobian_1(pde,un_coeff,z_coeff,z_param);
    // Integrate N * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *(fe_vol_->N()),
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*pde, theta_*dt);
    Intrepid::RealSpaceTools<Real>::add(*jac, *pde);
  }

  void Jacobian_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                   const ROL::TimeStamp<Real> & ts,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                   const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real one(1);
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew-told;
    // INITILAIZE JACOBIAN
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    ROL::Ptr<Intrepid::FieldContainer<Real>> pde;
    // ASSEMBLE OLD TIME JACOBIAN
    pde_->setTime(told);
    pde_->Jacobian_2(pde,uo_coeff,z_coeff,z_param);
    // ASSEMBLE NEW TIME JACOBIAN
    pde_->setTime(tnew);
    pde_->Jacobian_2(jac,un_coeff,z_coeff,z_param);
    // COMBINE JACOBIANS
    Intrepid::RealSpaceTools<Real>::scale(*pde, (one-theta_)*dt);
    Intrepid::RealSpaceTools<Real>::scale(*jac, theta_*dt);
    Intrepid::RealSpaceTools<Real>::add(*jac, *pde);
  }

  void Jacobian_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & jac,
                   const ROL::TimeStamp<Real> & ts,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                   const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Semilinear::Jacobian_zp): Jacobian is zero.");
  }

  void Hessian_uo_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real one(1);
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew-told;
    // COMPUTE OLD RESIDUAL
    pde_->setTime(told);
    pde_->Hessian_11(hess,l_coeff,uo_coeff,z_coeff,z_param);
    Intrepid::RealSpaceTools<Real>::scale(*hess, (one-theta_)*dt);
  }

  void Hessian_uo_un(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Semilinear::Hessian_uo_un): Hessian is zero.");
  }

  void Hessian_uo_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Semilinear::Hessian_uo_zf): Hessian is zero.");
  }

  void Hessian_uo_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Semilinear::Hessian_uo_zp): Hessian is zero.");
  }

  void Hessian_un_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Semilinear::Hessian_un_uo): Hessian is zero.");
  }

  void Hessian_un_un(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew-told;
    // COMPUTE OLD RESIDUAL
    pde_->setTime(tnew);
    pde_->Hessian_11(hess,l_coeff,un_coeff,z_coeff,z_param);
    Intrepid::RealSpaceTools<Real>::scale(*hess, theta_*dt);
  }

  void Hessian_un_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Semilinear::Hessian_un_zf): Hessian is zero.");
  }

  void Hessian_un_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Semilinear::Hessian_un_zp): Hessian is zero.");
  }

  void Hessian_zf_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Semilinear::Hessian_zf_uo): Hessian is zero.");
  }

  void Hessian_zf_un(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Semilinear::Hessian_zf_un): Hessian is zero.");
  }

  void Hessian_zf_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Semilinear::Hessian_zf_zf): Hessian is zero.");
  }

  void Hessian_zf_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Semilinear::Hessian_zf_zp): Hessian is zero.");
  }

  void Hessian_zp_uo(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Semilinear::Hessian_zp_uo): Hessian is zero.");
  }

  void Hessian_zp_un(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Semilinear::Hessian_zp_un): Hessian is zero.");
  }

  void Hessian_zp_zf(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Semilinear::Hessian_zp_zf): Hessian is zero.");
  }

  void Hessian_zp_zp(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Semilinear::Hessian_zp_zp): Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    pde_->RieszMap_1(riesz);
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    pde_->RieszMap_2(riesz);
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields() {
    return pde_->getFields();
  }

  void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real>> &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds) {
    pde_->setCellNodes(volCellNodes,bdryCellNodes,bdryCellLocIds);
    fe_vol_ = pde_->getFE();
  }

  const ROL::Ptr<FE<Real>> getFE(void) const {
    return fe_vol_;
  }
}; // DynamicPDE_Semilinear

#endif
