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

#ifndef PDEOPT_DYNAMIC_NAVIERSTOKES_HPP
#define PDEOPT_DYNAMIC_NAVIERSTOKES_HPP

#include "../../TOOLS/dynpde.hpp"
#include "../../TOOLS/fe.hpp"
#include "pde_navier-stokes.hpp"

template <class Real>
class DynamicPDE_NavierStokes : public DynamicPDE<Real> {
private:
  // Finite element basis information
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrs_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real> > volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;
  std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real>> feVel_;
  ROL::Ptr<FE<Real>> fePrs_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fvidx_;
  std::vector<std::vector<int>> fpidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellDofValues_;
  // Steady PDE without Dirichlet BC
  ROL::Ptr<PDE_NavierStokes<Real>> pde_;
  Real theta_;   // Time integration factor
  Real cx_, cy_; // Cylinder center
  Real r_;       // Cylinder radius
  bool useParametricControl_;
  bool useParabolicInflow_;
  bool useNonPenetratingWalls_;

  ROL::Ptr<FieldHelper<Real> > fieldHelper_;

public:
  DynamicPDE_NavierStokes(Teuchos::ParameterList &parlist) {
    pde_ = ROL::makePtr<PDE_NavierStokes<Real>>(parlist);
    // Problem data
    theta_ = parlist.sublist("Time Discretization").get("Theta", 1.0);
    useParametricControl_   = parlist.sublist("Problem").get("Use Parametric Control", false);
    useParabolicInflow_     = parlist.sublist("Problem").get("Use Parabolic Inflow", true);
    useNonPenetratingWalls_ = parlist.sublist("Problem").get("Use Non-Penetrating Walls", false);
    cx_                     = parlist.sublist("Problem").get("Cylinder Center X", 0.0);
    cy_                     = parlist.sublist("Problem").get("Cylinder Center Y", 0.0);
    r_                      = parlist.sublist("Problem").get("Cylinder Radius",   0.5);
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real>> & res,
                const ROL::TimeStamp<Real> & ts,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real one(1);
    // GET DIMENSIONS
    int  c = feVel_->gradN()->dimension(0);
    int fv = feVel_->gradN()->dimension(1);
    int  p = feVel_->gradN()->dimension(2);
    int  d = feVel_->gradN()->dimension(3);
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew-told;
    // INITIALIZE STORAGE
    ROL::Ptr<Intrepid::FieldContainer<Real>> r1, r2;
    Intrepid::FieldContainer<Real> velX_res(c, fv), velY_res(c, fv);
    ROL::Ptr<Intrepid::FieldContainer<Real>> valVelX_eval, valVelY_eval;
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U, Z, R;
    // COMPUTE STEADY RESIDUAL
    pde_->setTime(told);
    pde_->residual(r1,uo_coeff,z_coeff,z_param);
    Intrepid::RealSpaceTools<Real>::scale(*r1, (one-theta_)*dt);
    pde_->setTime(tnew);
    pde_->residual(r2,un_coeff,z_coeff,z_param);
    Intrepid::RealSpaceTools<Real>::scale(*r2, theta_*dt);
    Intrepid::RealSpaceTools<Real>::add(*r1, *r2);
    fieldHelper_->splitFieldCoeff(R, r1);
    // Integrate Uold * N
    valVelX_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    valVelY_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fieldHelper_->splitFieldCoeff(U, uo_coeff);
    feVel_->evaluateValue(valVelX_eval, U[0]);
    feVel_->evaluateValue(valVelY_eval, U[1]);
    Intrepid::FunctionSpaceTools::integrate<Real>(velX_res,
                                                  *valVelX_eval,        // UX
                                                  *(feVel_->NdetJ()),   // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(velY_res,
                                                  *valVelY_eval,        // UY
                                                  *(feVel_->NdetJ()),   // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::RealSpaceTools<Real>::subtract(*R[0], velX_res);
    Intrepid::RealSpaceTools<Real>::subtract(*R[1], velY_res);
    // Integrate Uold * N
    valVelX_eval->initialize();
    valVelY_eval->initialize();
    fieldHelper_->splitFieldCoeff(U, un_coeff);
    feVel_->evaluateValue(valVelX_eval, U[0]);
    feVel_->evaluateValue(valVelY_eval, U[1]);
    Intrepid::FunctionSpaceTools::integrate<Real>(velX_res,
                                                  *valVelX_eval,        // UX
                                                  *(feVel_->NdetJ()),   // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(velY_res,
                                                  *valVelY_eval,        // UY
                                                  *(feVel_->NdetJ()),   // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::RealSpaceTools<Real>::add(*R[0], velX_res);
    Intrepid::RealSpaceTools<Real>::add(*R[1], velY_res);
    // Apply boundary conditions
    int numSideSets = bdryCellLocIds_.size();
    Real bv(0);
    if (numSideSets > 0) {
      // APPLY DIRICHLET CONDITIONS
      for (int i = 0; i < numSideSets; ++i) {
        // Apply dirichlet conditions on inflow, top and bottom wall
        if (i==0 || i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                if (useNonPenetratingWalls_) {
                  (*R[1])(cidx,fvidx_[j][l]) = (*U[1])(cidx,fvidx_[j][l]);
                }
                else {
                  for (int m=0; m < d; ++m) {
                    bv = (*bdryCellDofValues_[i][j])(k,fvidx_[j][l],m);
                    (*R[m])(cidx,fvidx_[j][l]) = (*U[m])(cidx,fvidx_[j][l]) - bv;
                  }
                }
              }
            }
          }
        }
        // Apply dirichlet conditions on inflow, top and bottom wall
        if (i==3) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < d; ++m) {
                  bv = (*bdryCellDofValues_[i][j])(k,fvidx_[j][l],m);
                  (*R[m])(cidx,fvidx_[j][l]) = (*U[m])(cidx,fvidx_[j][l]) - bv;
                }
              }
            }
          }
        }
        // Apply Dirichlet control on cylinder
        if (i==4) {
          if (!useParametricControl_) {
            fieldHelper_->splitFieldCoeff(Z, z_coeff);
          }
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < d; ++m) {
                  if (!useParametricControl_) {
                    (*R[m])(cidx,fvidx_[j][l]) = (*U[m])(cidx,fvidx_[j][l]) - (*Z[m])(cidx,fvidx_[j][l]);
                  }
                  else {
                    bv  = (*bdryCellDofValues_[i][j])(k,fvidx_[j][l],m);
                    (*R[m])(cidx,fvidx_[j][l]) = (*U[m])(cidx,fvidx_[j][l]) - (*z_param)[0] * bv;
                  }
                }
              }
            }
          }
        }
      }
    }
    // Combine the residuals.
    fieldHelper_->combineFieldCoeff(res, R);
  }

  void Jacobian_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                   const ROL::TimeStamp<Real> & ts,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                   const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real zero(0), one(1);
    // GET DIMENSIONS
    int  c = feVel_->gradN()->dimension(0);
    int fv = feVel_->gradN()->dimension(1);
    int fp = fePrs_->gradN()->dimension(1);
    int  d = feVel_->gradN()->dimension(3);
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew-told;
    // INITIALIZE STORAGE
    ROL::Ptr<Intrepid::FieldContainer<Real>> j;
    Intrepid::FieldContainer<Real> velX_jac(c, fv, fv), velY_jac(c, fv, fv);
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J;
    // COMPUTE STEADY RESIDUAL
    pde_->setTime(told);
    pde_->Jacobian_1(j,uo_coeff,z_coeff,z_param);
    Intrepid::RealSpaceTools<Real>::scale(*j, (one-theta_)*dt);
    fieldHelper_->splitFieldCoeff(J, j);
    // Integrate Uold * N
    Intrepid::FunctionSpaceTools::integrate<Real>(velX_jac,
                                                  *(feVel_->N()),       // Phi
                                                  *(feVel_->NdetJ()),   // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(velY_jac,
                                                  *(feVel_->N()),       // Phi
                                                  *(feVel_->NdetJ()),   // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::RealSpaceTools<Real>::subtract(*J[0][0], velX_jac);
    Intrepid::RealSpaceTools<Real>::subtract(*J[1][1], velY_jac);
    // APPLY BOUNDARY CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==0 || i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < fv; ++m) {
                  if (useNonPenetratingWalls_) {
                    for (int p=0; p < d; ++p) {
                      (*J[1][p])(cidx,fvidx_[j][l],m) = zero;
                    }
                  }
                  else {
                    for (int n=0; n < d; ++n) {
                      for (int p=0; p < d; ++p) {
                        (*J[n][p])(cidx,fvidx_[j][l],m) = zero;
                      }
                    }
                  }
                }
                for (int m=0; m < fp; ++m) {
                  if (useNonPenetratingWalls_) {
                    (*J[1][2])(cidx,fvidx_[j][l],m) = zero;
                  }
                  else {
                    for (int n=0; n < d; ++n) {
                      (*J[n][2])(cidx,fvidx_[j][l],m) = zero;
                    }
                  }
                }
              }
            }
          }
        }
        if (i==3 || i==4) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < fv; ++m) {
                  for (int n=0; n < d; ++n) {
                    for (int p=0; p < d; ++p) {
                      (*J[n][p])(cidx,fvidx_[j][l],m) = zero;
                    }
                  }
                }
                for (int m=0; m < fp; ++m) {
                  for (int n=0; n < d; ++n) {
                    (*J[n][2])(cidx,fvidx_[j][l],m) = zero;
                  }
                }
              }
            }
          }
        }
      }
    }
    // Combine the Jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }

  void Jacobian_un(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                   const ROL::TimeStamp<Real> & ts,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                   const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real zero(0), one(1);
    // GET DIMENSIONS
    int  c = feVel_->gradN()->dimension(0);
    int fv = feVel_->gradN()->dimension(1);
    int fp = fePrs_->gradN()->dimension(1);
    int  d = feVel_->gradN()->dimension(3);
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew-told;
    // INITIALIZE STORAGE
    ROL::Ptr<Intrepid::FieldContainer<Real>> j;
    Intrepid::FieldContainer<Real> velX_jac(c, fv, fv), velY_jac(c, fv, fv);
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J;
    // COMPUTE STEADY RESIDUAL
    pde_->setTime(tnew);
    pde_->Jacobian_1(j,un_coeff,z_coeff,z_param);
    Intrepid::RealSpaceTools<Real>::scale(*j, theta_*dt);
    fieldHelper_->splitFieldCoeff(J, j);
    // Integrate Uold * N
    Intrepid::FunctionSpaceTools::integrate<Real>(velX_jac,
                                                  *(feVel_->N()),       // Phi
                                                  *(feVel_->NdetJ()),   // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(velY_jac,
                                                  *(feVel_->N()),       // Phi
                                                  *(feVel_->NdetJ()),   // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::RealSpaceTools<Real>::add(*J[0][0], velX_jac);
    Intrepid::RealSpaceTools<Real>::add(*J[1][1], velY_jac);
    // APPLY BOUNDARY CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==0 || i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < fv; ++m) {
                  if (useNonPenetratingWalls_) {
                    for (int p=0; p < d; ++p) {
                      (*J[1][p])(cidx,fvidx_[j][l],m) = zero;
                    }
                    (*J[1][1])(cidx,fvidx_[j][l],fvidx_[j][l]) = one;
                  }
                  else {
                    for (int n=0; n < d; ++n) {
                      for (int p=0; p < d; ++p) {
                        (*J[n][p])(cidx,fvidx_[j][l],m) = zero;
                      }
                      (*J[n][n])(cidx,fvidx_[j][l],fvidx_[j][l]) = one;
                    }
                  }
                }
                for (int m=0; m < fp; ++m) {
                  if (useNonPenetratingWalls_) {
                    (*J[1][2])(cidx,fvidx_[j][l],m) = zero;
                  }
                  else {
                    for (int n=0; n < d; ++n) {
                      (*J[n][2])(cidx,fvidx_[j][l],m) = zero;
                    }
                  }
                }
              }
            }
          }
        }
        if (i==3 || i==4) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < fv; ++m) {
                  for (int n=0; n < d; ++n) {
                    for (int p=0; p < d; ++p) {
                      (*J[n][p])(cidx,fvidx_[j][l],m) = zero;
                    }
                    (*J[n][n])(cidx,fvidx_[j][l],fvidx_[j][l]) = one;
                  }
                }
                for (int m=0; m < fp; ++m) {
                  for (int n=0; n < d; ++n) {
                    (*J[n][2])(cidx,fvidx_[j][l],m) = zero;
                  }
                }
              }
            }
          }
        }
      }
    }
    // Combine the Jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }

  void Jacobian_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                   const ROL::TimeStamp<Real> & ts,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                   const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (!useParametricControl_) {
      const Real one(1);
      // GET DIMENSIONS
      int  d = feVel_->gradN()->dimension(3);
      // INITILAIZE JACOBIAN
      pde_->Jacobian_2(jac,uo_coeff,z_coeff,z_param); // Resizes and zeros jac
      std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J;
      fieldHelper_->splitFieldCoeff(J, jac);
      // APPLY DIRICHLET CONDITIONS
      int numSideSets = bdryCellLocIds_.size();
      if (numSideSets > 0) {
        for (int i = 0; i < numSideSets; ++i) {
          if (i==4) {
            int numLocalSideIds = bdryCellLocIds_[i].size();
            for (int j = 0; j < numLocalSideIds; ++j) {
              int numCellsSide = bdryCellLocIds_[i][j].size();
              int numBdryDofs = fvidx_[j].size();
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < numBdryDofs; ++l) {
                  for (int m=0; m < d; ++m) {
                    (*J[m][m])(cidx,fvidx_[j][l],fvidx_[j][l]) = -one;
                  }
                }
              }
            }
          }
        }
      }
      // Combine the jacobians.
      fieldHelper_->combineFieldCoeff(jac, J);
    }
    else {
      throw Exception::Zero(">>> (PDE_NavierStokes::Jacobian_zf): Jacobian is zero.");
    }
  }

  void Jacobian_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & jac,
                   const ROL::TimeStamp<Real> & ts,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                   const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (useParametricControl_) {
      // GET DIMENSIONS
      int  d = feVel_->gradN()->dimension(3);
      // INITILAIZE JACOBIAN
      pde_->Jacobian_3(jac,uo_coeff,z_coeff,z_param); // Resizes and zeros jac
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> J;
      fieldHelper_->splitFieldCoeff(J, jac[0]);
      // APPLY DIRICHLET CONDITIONS
      int numSideSets = bdryCellLocIds_.size();
      if (numSideSets > 0) {
        for (int i = 0; i < numSideSets; ++i) {
          // Apply Dirichlet controls
          if (i==4) {
            Real bv(0);
            int numLocalSideIds = bdryCellLocIds_[i].size();
            for (int j = 0; j < numLocalSideIds; ++j) {
              int numCellsSide = bdryCellLocIds_[i][j].size();
              int numBdryDofs = fvidx_[j].size();
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < numBdryDofs; ++l) {
                  for (int m=0; m < d; ++m) {
                    bv = (*bdryCellDofValues_[i][j])(k,fvidx_[j][l],m);
                    (*J[m])(cidx,fvidx_[j][l]) = -bv;
                  }
                }
              }
            }
          }
        }
      }
      // Combine the jacobians.
      fieldHelper_->combineFieldCoeff(jac[0], J);
    }
    else {
      throw Exception::Zero(">>> (PDE_NavierStokes::Jacobian_zp): Jacobian is zero.");
    }
  }

  void Hessian_uo_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real zero(0), one(1);
    // GET DIMENSIONS
    int  d = feVel_->gradN()->dimension(3);
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew-told;
    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> L;
    ROL::Ptr<Intrepid::FieldContainer<Real>> l0_coeff;
    fieldHelper_->splitFieldCoeff(L, l_coeff);
    // Apply Dirichlet conditions
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==0 || i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                if (useNonPenetratingWalls_) {
                  (*L[1])(cidx,fvidx_[j][l]) = zero;
                }
                else {
                  for (int m=0; m < d; ++m) {
                    (*L[m])(cidx,fvidx_[j][l]) = zero;
                  }
                }
              }
            }
          }
        }
        if (i==3 || i==4) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < d; ++m) {
                  (*L[m])(cidx,fvidx_[j][l]) = zero;
                }
              }
            }
          }
        }
      }
    }
    fieldHelper_->combineFieldCoeff(l0_coeff, L);
    // COMPUTE STEADY RESIDUAL
    pde_->setTime(told);
    pde_->Hessian_11(hess,l0_coeff,uo_coeff,z_coeff,z_param);
    Intrepid::RealSpaceTools<Real>::scale(*hess, (one-theta_)*dt);
  }

  void Hessian_uo_un(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_NavierStokes::Hessian_uo_un): Hessian is zero.");
  }

  void Hessian_uo_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_NavierStokes::Hessian_uo_zf): Hessian is zero.");
  }

  void Hessian_uo_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_NavierStokes::Hessian_uo_zp): Hessian is zero.");
  }

  void Hessian_un_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_NavierStokes::Hessian_un_uo): Hessian is zero.");
  }

  void Hessian_un_un(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real zero(0);
    // GET DIMENSIONS
    int  d = feVel_->gradN()->dimension(3);
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew-told;
    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> L;
    ROL::Ptr<Intrepid::FieldContainer<Real>> l0_coeff;
    fieldHelper_->splitFieldCoeff(L, l_coeff);
    // Apply Dirichlet conditions
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==0 || i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                if (useNonPenetratingWalls_) {
                  (*L[1])(cidx,fvidx_[j][l]) = zero;
                }
                else {
                  for (int m=0; m < d; ++m) {
                    (*L[m])(cidx,fvidx_[j][l]) = zero;
                  }
                }
              }
            }
          }
        }
        if (i==3 || i==4) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < d; ++m) {
                  (*L[m])(cidx,fvidx_[j][l]) = zero;
                }
              }
            }
          }
        }
      }
    }
    fieldHelper_->combineFieldCoeff(l0_coeff, L);
    // COMPUTE STEADY RESIDUAL
    pde_->setTime(told);
    pde_->Hessian_11(hess,l0_coeff,un_coeff,z_coeff,z_param);
    Intrepid::RealSpaceTools<Real>::scale(*hess, theta_*dt);
  }

  void Hessian_un_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_NavierStokes::Hessian_un_zf): Hessian is zero.");
  }

  void Hessian_un_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_NavierStokes::Hessian_un_zp): Hessian is zero.");
  }

  void Hessian_zf_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_NavierStokes::Hessian_zf_uo): Hessian is zero.");
  }

  void Hessian_zf_un(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_NavierStokes::Hessian_zf_un): Hessian is zero.");
  }

  void Hessian_zf_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_NavierStokes::Hessian_zf_zf): Hessian is zero.");
  }

  void Hessian_zf_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_NavierStokes::Hessian_zf_zp): Hessian is zero.");
  }

  void Hessian_zp_uo(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_NavierStokes::Hessian_zp_uo): Hessian is zero.");
  }

  void Hessian_zp_un(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_NavierStokes::Hessian_zp_un): Hessian is zero.");
  }

  void Hessian_zp_zf(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_NavierStokes::Hessian_zp_zf): Hessian is zero.");
  }

  void Hessian_zp_zp(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_NavierStokes::Hessian_zp_zp): Hessian is zero.");
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
    volCellNodes_   = volCellNodes;
    bdryCellNodes_  = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    basisPtrs_      = pde_->getFields();
    feVel_          = pde_->getVelocityFE();
    fePrs_          = pde_->getPressureFE();
    fvidx_          = feVel_->getBoundaryDofs();
    fpidx_          = fePrs_->getBoundaryDofs();
    computeDirichlet();
  }

  void setFieldPattern(const std::vector<std::vector<int>> &fieldPattern) {
    pde_->setFieldPattern(fieldPattern);
    fieldHelper_ = pde_->getFieldHelper();
  }

  const ROL::Ptr<FE<Real>> getVelocityFE(void) const {
    return feVel_;
  }

  const ROL::Ptr<FE<Real>> getPressureFE(void) const {
    return fePrs_;
  }

  const ROL::Ptr<Intrepid::FieldContainer<Real> > getCellNodes(void) const {
    return pde_->getCellNodes();
  }

  const std::vector<ROL::Ptr<FE<Real>>> getVelocityBdryFE(const int sideset = -1) const {
    return pde_->getVelocityBdryFE(sideset);
  }

  const std::vector<std::vector<int>> getBdryCellLocIds(const int sideset = -1) const {
    int side = (sideset < 0 ? 4 : sideset);
    return bdryCellLocIds_[side];
  }

  const ROL::Ptr<FieldHelper<Real>> getFieldHelper(void) const {
    return fieldHelper_;
  }

private:
  Real DirichletFunc(const std::vector<Real> & coords, int sideset, int locSideId, int dir) const {
    const Real one(1), two(2);
    const Real x = coords[0], y = coords[1];
    Real val(0);
    if ((sideset!=4) && (dir==0)) {
      val = (useParabolicInflow_ ? (two + y) * (two - y) : one);
    }
    if (sideset==4) {
      Real tx = y-cy_, ty = cx_-x;
      val = (dir==0 ? tx : ty)/r_;
    }
    return val;
  }

  void computeDirichlet(void) {
    // Compute Dirichlet values at DOFs.
    int d = basisPtrs_[0]->getBaseCellTopology().getDimension();
    int numSidesets = bdryCellLocIds_.size();
    bdryCellDofValues_.resize(numSidesets);
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellDofValues_[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
        int f = basisPtrs_[0]->getCardinality();
        bdryCellDofValues_[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, d);
        ROL::Ptr<Intrepid::FieldContainer<Real> > coords =
          ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, d);
        if (c > 0) {
          feVel_->computeDofCoords(coords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<f; ++l) {
            std::vector<Real> dofpoint(d);
            for (int m=0; m<d; ++m) {
              dofpoint[m] = (*coords)(k, l, m);
            }
            for (int m=0; m<d; ++m) {
              (*bdryCellDofValues_[i][j])(k, l, m) = DirichletFunc(dofpoint, i, j, m);
            }
          }
        }
      }
    }
  }
}; // DynamicPDE_NavierStokes

#endif
