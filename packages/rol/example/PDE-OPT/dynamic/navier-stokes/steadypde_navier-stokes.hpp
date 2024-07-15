// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde_navier-stokes.hpp
    \brief Implements the local PDE interface for the Navier-Stokes control problem.
*/

#ifndef STEADYPDE_DYNAMIC_NAVIERSTOKES_HPP
#define STEADYPDE_DYNAMIC_NAVIERSTOKES_HPP

#include "pde_navier-stokes.hpp"

template <class Real>
class SteadyPDE_NavierStokes : public PDE<Real> {
private:
  // Cell node information
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellNodes_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellDofValues_;
  // Finite element definition
  ROL::Ptr<FE<Real>> feVel_;
  ROL::Ptr<FE<Real>> fePrs_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fvidx_;
  std::vector<std::vector<int>> fpidx_;
  // Field pattern, offsets, etc.
  ROL::Ptr<FieldHelper<Real>> fieldHelper_;

  ROL::Ptr<PDE_NavierStokes<Real>> pde_;

  Real cx_, cy_; // Cylinder center
  bool useParabolicInflow_;
  bool useNonPenetratingWalls_;

public:
  SteadyPDE_NavierStokes(Teuchos::ParameterList &parlist) {
    pde_ = ROL::makePtr<PDE_NavierStokes<Real>>(parlist);
    useParabolicInflow_     = parlist.sublist("Problem").get("Use Parabolic Inflow", true);
    useNonPenetratingWalls_ = parlist.sublist("Problem").get("Use Non-Penetrating Walls", false);
    cx_ = static_cast<Real>(-3);
    cy_ = static_cast<Real>(0);
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real>> & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    int  d = feVel_->gradN()->dimension(3);
    pde_->residual(res, u_coeff, z_coeff, z_param);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U, R;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(R, res);
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
          Real omega(0);
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < d; ++m) {
                  (*R[m])(cidx,fvidx_[j][l]) = (*U[m])(cidx,fvidx_[j][l]);
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

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real zero(0), one(1);
    int fv = feVel_->gradN()->dimension(1);
    int fp = fePrs_->gradN()->dimension(1);
    int d  = feVel_->gradN()->dimension(3);
    pde_->Jacobian_1(jac, u_coeff, z_coeff, z_param);
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J;
    fieldHelper_->splitFieldCoeff(J, jac);
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


  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (SteadyPDE_NavierStokes::Jacobian_2): Jacobian is zero.");
  }

  void Jacobian_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (SteadyPDE_NavierStokes::Jacobian_3): Jacobian is zero.");
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real zero(0);
    int d  = feVel_->gradN()->dimension(3);
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
    pde_->Hessian_11(hess, l0_coeff, u_coeff, z_coeff, z_param);
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (SteadyPDE_NavierStokes::Hessian_12): Hessian is zero.");
  }

  void Hessian_13(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (SteadyPDE_NavierStokes::Hessian_13): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (SteadyPDE_NavierStokes::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (SteadyPDE_NavierStokes::Hessian_22): Hessian is zero.");
  }

  void Hessian_23(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (SteadyPDE_NavierStokes::Hessian_23): Hessian is zero.");
  }

  void Hessian_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (SteadyPDE_NavierStokes::Hessian_31): Hessian is zero.");
  }

  void Hessian_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (SteadyPDE_NavierStokes::Hessian_32): Hessian is zero.");
  }

  void Hessian_33(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (SteadyPDE_NavierStokes::Hessian_33): Hessian is zero.");
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
    bdryCellNodes_  = bdryCellNodes;
    feVel_ = pde_->getVelocityFE();
    fePrs_ = pde_->getPressureFE();
    fvidx_ = feVel_->getBoundaryDofs();
    fpidx_ = fePrs_->getBoundaryDofs();
    computeDirichlet();
  }

  void setFieldPattern(const std::vector<std::vector<int>> & fieldPattern) {
    pde_->setFieldPattern(fieldPattern);
    fieldHelper_ = pde_->getFieldHelper();
  }

  const ROL::Ptr<FE<Real>> getVelocityFE(void) const {
    return feVel_;
  }

  const ROL::Ptr<FE<Real>> getPressureFE(void) const {
    return fePrs_;
  }

  const std::vector<ROL::Ptr<FE<Real>>> getVelocityBdryFE(void) const {
    return pde_->getVelocityBdryFE();
  }

  const std::vector<std::vector<int>> getBdryCellLocIds(const int sideset = -1) const {
    return pde_->getBdryCellLocIds(sideset);
  }

  const ROL::Ptr<FieldHelper<Real>> getFieldHelper(void) const {
    return fieldHelper_;
  }

private:
  Real DirichletFunc(const std::vector<Real> & coords, int sideset, int locSideId, int dir) const {
    const Real one(1);
    const Real x = coords[0], y = coords[1];
    Real val(0);
    if ((sideset!=4) && (dir==0)) {
      val = (useParabolicInflow_ ? (one + y) * (one - y) : one);
    }
    if (sideset==4) {
      Real tx = y-cy_, ty = cx_-x;
      val = (dir==0 ? tx : ty);
    }
    return val;
  }

  void computeDirichlet(void) {
    // Compute Dirichlet values at DOFs.
    int d = feVel_->gradN()->dimension(3);
    int numSidesets = bdryCellLocIds_.size();
    bdryCellDofValues_.resize(numSidesets);
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellDofValues_[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
        int f = (pde_->getFields())[0]->getCardinality();
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
}; // PDE_NavierStokes

#endif
