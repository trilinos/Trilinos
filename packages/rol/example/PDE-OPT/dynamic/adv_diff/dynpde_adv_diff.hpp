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

#ifndef PDEOPT_DYNAMIC_ADV_DIFF_HPP
#define PDEOPT_DYNAMIC_ADV_DIFF_HPP

#include "../../TOOLS/dynpde.hpp"
#include "../../TOOLS/fe.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"

#include "pde_adv_diff.hpp"

template <class Real>
class DynamicPDE_adv_diff : public DynamicPDE<Real> {
private:
  // Cell node information
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real>> fe_vol_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellDofValues_;

  // Steady PDE without Dirichlet BC
  ROL::Ptr<PDE_adv_diff<Real>> pde_;
  Real theta_;
  Real useDBC_;

public:
  DynamicPDE_adv_diff(Teuchos::ParameterList &parlist) {
    pde_ = ROL::makePtr<PDE_adv_diff<Real>>(parlist);
    // Time-dependent coefficients
    theta_  = parlist.sublist("Time Discretization").get("Theta",1.0);
    useDBC_ = parlist.sublist("Problem").get("Use Dirichlet",false);
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
    // APPLY DIRICHLET CONDITIONS
    if (useDBC_) {
      int numLocalSideIds = bdryCellLocIds_[0].size();
      for (int j = 0; j < numLocalSideIds; ++j) {
        int numCellsSide = bdryCellLocIds_[0][j].size();
        int numBdryDofs = fidx_[j].size();
        for (int k = 0; k < numCellsSide; ++k) {
          int cidx = bdryCellLocIds_[0][j][k];
          for (int l = 0; l < numBdryDofs; ++l) {
            (*res)(cidx,fidx_[j][l])
              = (*un_coeff)(cidx,fidx_[j][l]) - (*bdryCellDofValues_[0][j])(k,fidx_[j][l]);
          }
        }
      }
    }
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
    // APPLY DIRICHLET CONDITIONS
    if (useDBC_) {
      int numLocalSideIds = bdryCellLocIds_[0].size();
      for (int j = 0; j < numLocalSideIds; ++j) {
        int numCellsSide = bdryCellLocIds_[0][j].size();
        int numBdryDofs = fidx_[j].size();
        for (int k = 0; k < numCellsSide; ++k) {
          int cidx = bdryCellLocIds_[0][j][k];
          for (int l = 0; l < numBdryDofs; ++l) {
            //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
            for (int m = 0; m < f; ++m) {
              (*jac)(cidx,fidx_[j][l],m) = static_cast<Real>(0);
            }
          }
        }
      }
    }
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
    // APPLY DIRICHLET CONDITIONS
    if (useDBC_) {
      int numLocalSideIds = bdryCellLocIds_[0].size();
      for (int j = 0; j < numLocalSideIds; ++j) {
        int numCellsSide = bdryCellLocIds_[0][j].size();
        int numBdryDofs = fidx_[j].size();
        for (int k = 0; k < numCellsSide; ++k) {
          int cidx = bdryCellLocIds_[0][j][k];
          for (int l = 0; l < numBdryDofs; ++l) {
            //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
            for (int m = 0; m < f; ++m) {
              (*jac)(cidx,fidx_[j][l],m) = static_cast<Real>(0);
            }
            (*jac)(cidx,fidx_[j][l],fidx_[j][l]) = static_cast<Real>(1);
          }
        }
      }
    }
  }

  void Jacobian_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                   const ROL::TimeStamp<Real> & ts,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                   const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Jacobian_zf): Jacobian is zero.");
  }

  void Jacobian_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & jac,
                   const ROL::TimeStamp<Real> & ts,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                   const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real one(1);
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew-told;
    // ADD CONTROL TERM TO RESIDUAL
    int size = z_param->size();
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> pde(size,ROL::nullPtr);
    // EVALUATE OLD COMPONENT
    pde_->setTime(told); 
    pde_->Jacobian_3(jac,uo_coeff,z_coeff,z_param);
    // EVALUATE NEW COMPONENT
    pde_->setTime(tnew);
    pde_->Jacobian_3(pde,un_coeff,z_coeff,z_param);
    for (int i = 0; i < size; ++i) {
      Intrepid::RealSpaceTools<Real>::scale(*(jac[i]), (one-theta_)*dt);
      Intrepid::RealSpaceTools<Real>::scale(*(pde[i]), theta_*dt);
      Intrepid::RealSpaceTools<Real>::add(*(jac[i]), *(pde[i]));
      // APPLY DIRICHLET CONDITIONS
      if (useDBC_) {
        int numLocalSideIds = bdryCellLocIds_[0].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          int numCellsSide = bdryCellLocIds_[0][j].size();
          int numBdryDofs = fidx_[j].size();
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[0][j][k];
            for (int l = 0; l < numBdryDofs; ++l) {
              (*(jac[i]))(cidx,fidx_[j][l]) = static_cast<Real>(0);
            }
          }
        }
      }
    }
  }

  void Hessian_uo_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_uo_uo): Hessian is zero.");
  }

  void Hessian_uo_un(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_uo_un): Hessian is zero.");
  }

  void Hessian_uo_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_uo_zf): Hessian is zero.");
  }

  void Hessian_uo_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_uo_zp): Hessian is zero.");
  }

  void Hessian_un_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_un_uo): Hessian is zero.");
  }

  void Hessian_un_un(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_un_un): Hessian is zero.");
  }

  void Hessian_un_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_un_zf): Hessian is zero.");
  }

  void Hessian_un_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_un_zp): Hessian is zero.");
  }

  void Hessian_zf_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_zf_uo): Hessian is zero.");
  }

  void Hessian_zf_un(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_zf_un): Hessian is zero.");
  }

  void Hessian_zf_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_zf_zf): Hessian is zero.");
  }

  void Hessian_zf_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_zf_zp): Hessian is zero.");
  }

  void Hessian_zp_uo(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_zp_uo): Hessian is zero.");
  }

  void Hessian_zp_un(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_zp_un): Hessian is zero.");
  }

  void Hessian_zp_zf(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_zp_zf): Hessian is zero.");
  }

  void Hessian_zp_zp(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_adv_diff::Hessian_zp_zp): Hessian is zero.");
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
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fe_vol_ = pde_->getFE();
    if (useDBC_) {
      // Set local boundary DOFs.
      fidx_ = fe_vol_->getBoundaryDofs();
      // Compute Dirichlet values at DOFs.
      int d = pde_->getFields()[0]->getBaseCellTopology().getDimension();
      int numSidesets = bdryCellLocIds_.size();
      bdryCellDofValues_.resize(numSidesets);
      for (int i=0; i<numSidesets; ++i) {
        int numLocSides = bdryCellLocIds_[i].size();
        bdryCellDofValues_[i].resize(numLocSides);
        for (int j=0; j<numLocSides; ++j) {
          int c = bdryCellLocIds_[i][j].size();
          int f = pde_->getFields()[0]->getCardinality();
          bdryCellDofValues_[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
          ROL::Ptr<Intrepid::FieldContainer<Real>> coords =
            ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, d);
          if (c > 0) {
            fe_vol_->computeDofCoords(coords, bdryCellNodes[i][j]);
          }
          for (int k=0; k<c; ++k) {
            for (int l=0; l<f; ++l) {
              std::vector<Real> dofpoint(d);
              for (int m=0; m<d; ++m) {
                dofpoint[m] = (*coords)(k, l, m);
              }
              (*bdryCellDofValues_[i][j])(k, l) = evaluateDirichlet(dofpoint, i, j);
            }
          }
        }
      }
    }
  }

  const ROL::Ptr<FE<Real>> getFE(void) const {
    return fe_vol_;
  }

private:
  Real evaluateDirichlet(const std::vector<Real> & coords, int sideset, int locSideId) const {
    return static_cast<Real>(0);
  }
}; // DynamicPDE_adv_diff

#endif
