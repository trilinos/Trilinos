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

#ifndef PDEOPT_DYNAMIC_THERMALFLUIDS_HPP
#define PDEOPT_DYNAMIC_THERMALFLUIDS_HPP

#include "../../TOOLS/dynpde.hpp"
#include "../../TOOLS/fe.hpp"
#include "pde_thermal-fluids.hpp"

template <class Real>
class DynamicPDE_ThermalFluids : public DynamicPDE<Real> {
private:
  // Finite element basis information
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrs_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real>> volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real> > feVel_;
  ROL::Ptr<FE<Real> > fePrs_;
  ROL::Ptr<FE<Real> > feThr_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fvidx_;
  std::vector<std::vector<int>> fpidx_;
  std::vector<std::vector<int>> fhidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellVDofValues_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellTDofValues_;
  // Steady PDE without Dirichlet BC
  ROL::Ptr<PDE_ThermalFluids<Real>> pde_;
  Real theta_;   // Time integration factor
  Real cx_, cy_; // Cylinder center
  Real r_;       // Cylinder radius

  ROL::Ptr<FieldHelper<Real>> fieldHelper_;

public:
  DynamicPDE_ThermalFluids(Teuchos::ParameterList &parlist) {
    pde_ = ROL::makePtr<PDE_ThermalFluids<Real>>(parlist);
    // Problem data
    theta_ = parlist.sublist("Time Discretization").get("Theta", 1.0);
    cx_    = parlist.sublist("Problem").get("Cylinder Center X", 0.0);
    cy_    = parlist.sublist("Problem").get("Cylinder Center Y", 0.0);
    r_     = parlist.sublist("Problem").get("Cylinder Radius",   0.5);
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real>> & res,
                const ROL::TimeStamp<Real> & ts,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real one(1);
    // GET DIMENSIONS
    const int  c = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fh = feThr_->gradN()->dimension(1);
    const int  p = feVel_->gradN()->dimension(2);
    const int  d = feVel_->gradN()->dimension(3);
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew-told;
    // INITIALIZE STORAGE
    ROL::Ptr<Intrepid::FieldContainer<Real>> r1, r2;
    Intrepid::FieldContainer<Real> velX_res(c, fv), velY_res(c, fv), thr_res(c, fh);
    ROL::Ptr<Intrepid::FieldContainer<Real>> valVelX_eval, valVelY_eval, valThr_eval;
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
    valThr_eval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fieldHelper_->splitFieldCoeff(U, uo_coeff);
    feVel_->evaluateValue(valVelX_eval, U[0]);
    feVel_->evaluateValue(valVelY_eval, U[1]);
    feThr_->evaluateValue(valThr_eval,  U[3]);
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
    Intrepid::FunctionSpaceTools::integrate<Real>(thr_res,
                                                  *valThr_eval,         // T
                                                  *(feThr_->NdetJ()),   // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::RealSpaceTools<Real>::subtract(*R[0], velX_res);
    Intrepid::RealSpaceTools<Real>::subtract(*R[1], velY_res);
    Intrepid::RealSpaceTools<Real>::subtract(*R[3], thr_res);
    // Integrate Unew * N
    valVelX_eval->initialize();
    valVelY_eval->initialize();
    valThr_eval->initialize();
    fieldHelper_->splitFieldCoeff(U, un_coeff);
    feVel_->evaluateValue(valVelX_eval, U[0]);
    feVel_->evaluateValue(valVelY_eval, U[1]);
    feThr_->evaluateValue(valThr_eval,  U[3]);
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
    Intrepid::FunctionSpaceTools::integrate<Real>(thr_res,
                                                  *valThr_eval,         // T
                                                  *(feThr_->NdetJ()),   // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::RealSpaceTools<Real>::add(*R[0], velX_res);
    Intrepid::RealSpaceTools<Real>::add(*R[1], velY_res);
    Intrepid::RealSpaceTools<Real>::add(*R[3], thr_res);
    // Apply boundary conditions
    // -- Velocity:
    //    -- Free Stream: i=0, i=2 (top and bottom wall)
    //    --      Inflow: i=3 (left wall)
    //    --     Outflow: i=1 (right wall)
    // -- Thermal:
    //      --     Robin: i=0, i=2, i=4 (top and bottom wall and cylinder)
    //      -- Dirichlet: i=3 (left wall)
    //      --   Neumann: i=1 (right wall)
    int numSideSets = bdryCellLocIds_.size();
    Real bv(0);
    if (numSideSets > 0) {
      // APPLY BOUNDARY CONDITIONS
      for (int i = 0; i < numSideSets; ++i) {
        // VELOCITY BOUNDARY CONDITIONS
        if (i==0 || i==2 || i==3 ) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < d; ++m) {
                  bv = (*bdryCellVDofValues_[i][j])(k,fvidx_[j][l],m);
                  (*R[m])(cidx,fvidx_[j][l]) = (*U[m])(cidx,fvidx_[j][l]) - bv;
                }
              }
            }
          }
        }
        // Apply Dirichlet control on cylinder
        if (i==4) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < d; ++m) {
                  bv = (*bdryCellVDofValues_[i][j])(k,fvidx_[j][l],m);
                  (*R[m])(cidx,fvidx_[j][l]) = (*U[m])(cidx,fvidx_[j][l]) - (*z_param)[0] * bv;
                }
              }
            }
          }
        }
        // THERMAL BOUNDARY CONDITIONS
        // Apply Dirichlet conditions on bottom wall
        if (i==3) {
          int numLocalSideIds = bdryCellLocIds_[i].size();   // Number of sides per cell
          for (int j = 0; j < numLocalSideIds; ++j) {        // Loop over sides of cell: Quad = {0, 1, 2, 3}
            int numCellsSide = bdryCellLocIds_[i][j].size(); // Number of cells with side j
            int numHBdryDofs = fhidx_[j].size();             // Number of thermal boundary DOFs
            for (int k = 0; k < numCellsSide; ++k) {         // Loop over cells with side j
              int cidx = bdryCellLocIds_[i][j][k];           // Cell index
              for (int l = 0; l < numHBdryDofs; ++l) {       // Loop over all fields of cell k on side j
                bv = (*bdryCellTDofValues_[i][j])(k,fhidx_[j][l]);
                (*R[d+1])(cidx,fhidx_[j][l]) = (*U[d+1])(cidx,fhidx_[j][l]) - bv;
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
    const int  c = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int fh = feThr_->gradN()->dimension(1);
    const int  d = feVel_->gradN()->dimension(3);
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew-told;
    // INITIALIZE STORAGE
    ROL::Ptr<Intrepid::FieldContainer<Real>> j;
    Intrepid::FieldContainer<Real> velX_jac(c, fv, fv), velY_jac(c, fv, fv), thr_jac(c, fh, fh);
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
    Intrepid::FunctionSpaceTools::integrate<Real>(thr_jac,
                                                  *(feThr_->N()),       // Phi
                                                  *(feThr_->NdetJ()),   // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::RealSpaceTools<Real>::subtract(*J[0][0], velX_jac);
    Intrepid::RealSpaceTools<Real>::subtract(*J[1][1], velY_jac);
    Intrepid::RealSpaceTools<Real>::subtract(*J[3][3], thr_jac);
    // Apply boundary conditions
    // -- Velocity:
    //    -- Free Stream: i=0, i=2 (top and bottom wall)
    //    --      Inflow: i=3 (left wall)
    //    --     Outflow: i=1 (right wall)
    // -- Thermal:
    //      --     Robin: i=0, i=2, i=4 (top and bottom wall and cylinder)
    //      -- Dirichlet: i=3 (left wall)
    //      --   Neumann: i=1 (right wall)
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        // APPLY VELOCITY BOUNDARY CONDITIONS
        if (i==0 || i==2 || i==3 || i==4) {
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
                for (int n=0; n < d; ++n) {
                  for (int m=0; m < fp; ++m) {
                    (*J[n][d])(cidx,fvidx_[j][l],m) = zero;
                  }
                  for (int m=0; m < fh; ++m) {
                    (*J[n][d+1])(cidx,fvidx_[j][l],m) = zero;
                  }
                }
              }
            }
          }
        }
        if (i==3) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numHBdryDofs = fhidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numHBdryDofs; ++l) {
                for (int m = 0; m < fv; ++m) {
                  for (int n = 0; n < d; ++n) {
                    (*J[d+1][n])(cidx,fhidx_[j][l],m) = zero;
                  }
                }
                for (int m = 0; m < fp; ++m) {
                  (*J[d+1][d])(cidx,fhidx_[j][l],m) = zero;
                }
                for (int m = 0; m < fh; ++m) {
                  (*J[d+1][d+1])(cidx,fhidx_[j][l],m) = zero;
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
    const int  c = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int fh = feThr_->gradN()->dimension(1);
    const int  d = feVel_->gradN()->dimension(3);
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew-told;
    // INITIALIZE STORAGE
    ROL::Ptr<Intrepid::FieldContainer<Real>> j;
    Intrepid::FieldContainer<Real> velX_jac(c, fv, fv), velY_jac(c, fv, fv), thr_jac(c, fh, fh);
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
    Intrepid::FunctionSpaceTools::integrate<Real>(thr_jac,
                                                  *(feThr_->N()),       // Phi
                                                  *(feThr_->NdetJ()),   // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::RealSpaceTools<Real>::add(*J[0][0], velX_jac);
    Intrepid::RealSpaceTools<Real>::add(*J[1][1], velY_jac);
    Intrepid::RealSpaceTools<Real>::add(*J[3][3], thr_jac);
    // APPLY BOUNDARY CONDITIONS
    // -- Velocity:
    //    -- Free Stream: i=0, i=2 (top and bottom wall)
    //    --      Inflow: i=3 (left wall)
    //    --     Outflow: i=1 (right wall)
    // -- Thermal:
    //      --     Robin: i=0, i=2, i=4 (top and bottom wall and cylinder)
    //      -- Dirichlet: i=3 (left wall)
    //      --   Neumann: i=1 (right wall)
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==0 || i==2 || i==3 || i==4) {
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
                for (int n=0; n < d; ++n) {
                  for (int m=0; m < fp; ++m) {
                    (*J[n][d])(cidx,fvidx_[j][l],m) = zero;
                  }
                  for (int m=0; m < fh; ++m) {
                    (*J[n][d+1])(cidx,fvidx_[j][l],m) = zero;
                  }
                }
              }
            }
          }
        }
        if (i==3) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numHBdryDofs = fhidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numHBdryDofs; ++l) {
                for (int m = 0; m < fv; ++m) {
                  for (int n = 0; n < d; ++n) {
                    (*J[d+1][n])(cidx,fhidx_[j][l],m) = zero;
                  }
                }
                for (int m = 0; m < fp; ++m) {
                  (*J[d+1][d])(cidx,fhidx_[j][l],m) = zero;
                }
                for (int m = 0; m < fh; ++m) {
                  (*J[d+1][d+1])(cidx,fhidx_[j][l],m) = zero;
                }
                (*J[d+1][d+1])(cidx,fhidx_[j][l],fhidx_[j][l]) = one;
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
    throw Exception::Zero(">>> (PDE_NavierStokes::Jacobian_zf): Jacobian is zero.");
  }

  void Jacobian_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & jac,
                   const ROL::TimeStamp<Real> & ts,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                   const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int  d = feVel_->gradN()->dimension(3);
    // INITILAIZE JACOBIAN
    pde_->Jacobian_3(jac,uo_coeff,z_coeff,z_param); // Resizes and zeros jac
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> J;
    fieldHelper_->splitFieldCoeff(J, jac[0]);
    // APPLY DIRICHLET CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      Real bv(0);
      for (int i = 0; i < numSideSets; ++i) {
        // Apply Dirichlet controls
        if (i==4) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < d; ++m) {
                  bv = (*bdryCellVDofValues_[i][j])(k,fvidx_[j][l],m);
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

  void Hessian_uo_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real zero(0), one(1);
    // GET DIMENSIONS
    const int  d = feVel_->gradN()->dimension(3);
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
        // Velocity boundary conditions
        if (i==0 || i==2 || i==3 || i==4) {
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
        // Thermal boundary conditions
        if (i==3) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numHBdryDofs = fhidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numHBdryDofs; ++l) {
                (*L[d+1])(cidx,fhidx_[j][l]) = zero;
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
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_uo_un): Hessian is zero.");
  }

  void Hessian_uo_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_uo_zf): Hessian is zero.");
  }

  void Hessian_uo_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_uo_zp): Hessian is zero.");
  }

  void Hessian_un_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_un_uo): Hessian is zero.");
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
    const int  d = feVel_->gradN()->dimension(3);
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
        // Velocity boundary conditions
        if (i==0 || i==2 || i==3 || i==4) {
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
        // Thermal boundary conditions
        if (i==3) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numHBdryDofs = fhidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numHBdryDofs; ++l) {
                (*L[d+1])(cidx,fhidx_[j][l]) = zero;
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
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_un_zf): Hessian is zero.");
  }

  void Hessian_un_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_un_zp): Hessian is zero.");
  }

  void Hessian_zf_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_zf_uo): Hessian is zero.");
  }

  void Hessian_zf_un(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_zf_un): Hessian is zero.");
  }

  void Hessian_zf_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_zf_zf: Hessian is zero.");
  }

  void Hessian_zf_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_zf_zp): Hessian is zero.");
  }

  void Hessian_zp_uo(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_zp_uo): Hessian is zero.");
  }

  void Hessian_zp_un(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_zp_un): Hessian is zero.");
  }

  void Hessian_zp_zf(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_zp_zf): Hessian is zero.");
  }

  void Hessian_zp_zp(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_zp_zp): Hessian is zero.");
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
    feThr_          = pde_->getThermalFE();
    fvidx_          = feVel_->getBoundaryDofs();
    fpidx_          = fePrs_->getBoundaryDofs();
    fhidx_          = feThr_->getBoundaryDofs();
    computeDirichlet();
  }

  void setFieldPattern(const std::vector<std::vector<int>> &fieldPattern) {
    pde_->setFieldPattern(fieldPattern);
    fieldHelper_ = pde_->getFieldHelper();
  }

  const ROL::Ptr<Intrepid::FieldContainer<Real> > getCellNodes(void) const {
    return pde_->getCellNodes();
  }

  const ROL::Ptr<FE<Real>> getVelocityFE(void) const {
    return feVel_;
  }

  const ROL::Ptr<FE<Real>> getPressureFE(void) const {
    return fePrs_;
  }

  const ROL::Ptr<FE<Real>> getThermalFE(void) const {
    return feThr_;
  }

  const std::vector<std::vector<int>> getBdryCellLocIds(const int sideset = -1) const {
    int side = (sideset < 0 ? 4 : sideset);
    return bdryCellLocIds_[side];
  }

  const std::vector<ROL::Ptr<FE<Real>>> getVelocityBdryFE(const int sideset = -1) const {
    int side = (sideset < 0 ? 4 : sideset);
    return pde_->getVelocityBdryFE(side);
  }

  const std::vector<ROL::Ptr<FE<Real>>> getPressureBdryFE(const int sideset = -1) const {
    int side = (sideset < 0 ? 4 : sideset);
    return pde_->getPressureBdryFE(side);
  }

  const std::vector<ROL::Ptr<FE<Real>>> getThermalBdryFE(const int sideset = -1) const {
    int side = (sideset < 0 ? 4 : sideset);
    return pde_->getThermalBdryFE(side);
  }

  const ROL::Ptr<FieldHelper<Real>> getFieldHelper(void) const {
    return fieldHelper_;
  }

private:
  Real velocityDirichletFunc(const std::vector<Real> & coords, int sideset, int locSideId, int dir) const {
    const Real one(1), zero(0);
    if (sideset==4) {
      Real x  = coords[0];
      Real y  = coords[1];
      Real tx = y-cy_, ty = cx_-x;
      return (dir==0 ? tx : ty)/r_;
    }
    return (dir==0 ? one : zero);
  }

  Real thermalDirichletFunc(const std::vector<Real> & coords, int sideset, int locSideId) const {
    const Real zero(0);
    return zero;
  }

  void computeDirichlet(void) {
    // Compute Dirichlet values at DOFs.
    const int fv = feVel_->gradN()->dimension(1);
    const int fh = feThr_->gradN()->dimension(1);
    const int  d = feVel_->gradN()->dimension(3);
    int numSidesets = bdryCellLocIds_.size();
    bdryCellVDofValues_.resize(numSidesets);
    bdryCellTDofValues_.resize(numSidesets);
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellVDofValues_[i].resize(numLocSides);
      bdryCellTDofValues_[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
        bdryCellVDofValues_[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, d);
        bdryCellTDofValues_[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);
        ROL::Ptr<Intrepid::FieldContainer<Real> > Vcoords, Tcoords;
        Vcoords = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, d);
        Tcoords = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, d);
        if (c > 0) {
          feVel_->computeDofCoords(Vcoords, bdryCellNodes_[i][j]);
          feThr_->computeDofCoords(Tcoords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<fv; ++l) {
            std::vector<Real> dofpoint(d);
            for (int m=0; m<d; ++m) {
              dofpoint[m] = (*Vcoords)(k, l, m);
            }

            for (int m=0; m<d; ++m) {
              (*bdryCellVDofValues_[i][j])(k, l, m) = velocityDirichletFunc(dofpoint, i, j, m);
            }
          }
          for (int l=0; l<fh; ++l) {
            std::vector<Real> dofpoint(d);
            for (int m=0; m<d; ++m) {
              dofpoint[m] = (*Tcoords)(k, l, m);
            }
            (*bdryCellTDofValues_[i][j])(k, l) = thermalDirichletFunc(dofpoint, i, j);
          }
        }
      }
    }
  }
}; // DynamicPDE_ThermalFluids

#endif
