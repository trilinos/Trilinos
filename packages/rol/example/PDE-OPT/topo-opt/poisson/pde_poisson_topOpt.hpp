// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde_poisson_topOpt.hpp
    \brief Implements the local PDE interface for the Poisson
           topology optimization problem.
*/

#ifndef PDE_POISSON_TOPOPT_HPP
#define PDE_POISSON_TOPOPT_HPP

#include "../../TOOLS/pde.hpp"
#include "../../TOOLS/fe.hpp"
#include "../../TOOLS/Intrepid_HGRAD_C0_FEM.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real>
class PDE_Poisson_TopOpt : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>> basisPtr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrs_;
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>> basisPtrDens_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrsDens_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real>> cellCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real>> volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real>> fe_vol_, fe_dens_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellDofValues_;
  // Force function evaluated at cubature points
  ROL::Ptr<Intrepid::FieldContainer<Real>> force_eval_;
  ROL::Ptr<Intrepid::FieldContainer<Real>> nforce_eval_;
  // Inputs
  Real minConductivity_;
  Real SIMPpower_;

  bool getFields2called_;

  Real ForceFunc(const std::vector<Real> &x) const {
    return static_cast<Real>(0.01);
  }

  Real dirichletFunc(const std::vector<Real> & coords, const int sideset, const int locSideId) const {
    return static_cast<Real>(0);
  }

public:
  PDE_Poisson_TopOpt(Teuchos::ParameterList &parlist) : getFields2called_(false) {
    // Finite element fields.
    int basisOrder     = parlist.sublist("Problem").get("Order of FE discretization",1);
    int basisOrderDens = parlist.sublist("Problem").get("Density Basis Order",0);
    TEUCHOS_TEST_FOR_EXCEPTION(basisOrder > 2 || basisOrder < 1, std::invalid_argument,
      ">>> PDE-OPT/topopt/poisson/pde_poisson_topOpt.hpp: Basis order is not 1 or 2!");
    TEUCHOS_TEST_FOR_EXCEPTION(basisOrderDens > 1 || basisOrderDens < 0, std::invalid_argument,
      ">>> PDE-OPT/topopt/poisson/pde_poisson_topOpt.hpp: Basis order is not 0 or 1!");
    if (basisOrder == 1)
      basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real,Intrepid::FieldContainer<Real>>>();
    else if (basisOrder == 2)
      basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real,Intrepid::FieldContainer<Real>>>();
    if (basisOrderDens == 1)
      basisPtrDens_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real,Intrepid::FieldContainer<Real>>>();
    else
      basisPtrDens_ = ROL::makePtr<Intrepid::Basis_HGRAD_C0_FEM<Real,Intrepid::FieldContainer<Real>>>();
    basisPtrs_.clear();     basisPtrs_.push_back(basisPtr_);
    basisPtrsDens_.clear(); basisPtrsDens_.push_back(basisPtrDens_); // Density component
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();        // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                       // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree",2);     // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                       // create default cubature

    minConductivity_ = parlist.sublist("Problem").get("Minimum Conductivity",1.e-3);
    SIMPpower_       = parlist.sublist("Problem").get("SIMP Power",3.0);
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real>> & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real one(1);
    // Get dimensions
    const int c = fe_vol_->gradN()->dimension(0);
    const int f = fe_vol_->gradN()->dimension(1);
    const int p = fe_vol_->gradN()->dimension(2);
    const int d = fe_vol_->gradN()->dimension(3);
    // Initialize residual storage
    res = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Temporary storage
    ROL::Ptr<Intrepid::FieldContainer<Real>> valZ_eval, gradU_eval, KgradU;
    valZ_eval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    gradU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    KgradU     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    // Build SIMP density at cubature points
    fe_dens_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valZ_eval)(i,j) = minConductivity_
          + (one - minConductivity_) * std::pow((*valZ_eval)(i,j),SIMPpower_);
      }
    }
    // Build flux function
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*KgradU,
                                                               *valZ_eval,
                                                               *gradU_eval);
    // Integrate stiffness term
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *KgradU,
                                                  *(fe_vol_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // Add force term
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *nforce_eval_,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // Apply Dirichlet boundary conditions
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                (*res)(cidx,fidx_[j][l]) = (*u_coeff)(cidx,fidx_[j][l]) - (*bdryCellDofValues_[i][j])(k,fidx_[j][l]);
              }
            }
          }
        }
      }
    }
  }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real one(1);
    // Get dimensions
    const int c = fe_vol_->gradN()->dimension(0);
    const int f = fe_vol_->gradN()->dimension(1);
    const int p = fe_vol_->gradN()->dimension(2);
    const int d = fe_vol_->gradN()->dimension(3);
    // Initialize Jacobian storage
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    ROL::Ptr<Intrepid::FieldContainer<Real>> valZ_eval, KgradN;
    valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    KgradN    = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p, d);
    // Build density-dependent conductivity function
    fe_dens_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valZ_eval)(i,j) = minConductivity_
          + (one - minConductivity_) * std::pow((*valZ_eval)(i,j),SIMPpower_);
      }
    }
    // Build flux function
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*KgradN,
                                                                *valZ_eval,
                                                                *(fe_vol_->gradN()));
    // Integrate stiffness term
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *KgradN,
                                                  *(fe_vol_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);

    // Apply Dirichlet boundary conditions
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m = 0; m < f; ++m) {
                  (*jac)(cidx,fidx_[j][l],m) = static_cast<Real>(0);
                }
                (*jac)(cidx,fidx_[j][l],fidx_[j][l]) = static_cast<Real>(1);
              }
            }
          }
        }
      }
    }
  }

  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real one(1);
    // Get dimensions
    const int c  = fe_vol_->gradN()->dimension(0);
    const int f  = fe_vol_->gradN()->dimension(1);
    const int fd = fe_dens_->gradN()->dimension(1);
    const int p  = fe_vol_->gradN()->dimension(2);
    const int d  = fe_vol_->gradN()->dimension(3);
    // Initialize Jacobian storage
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, fd);
    ROL::Ptr<Intrepid::FieldContainer<Real>> valZ_eval, gradU_eval, dKN, gradUgradN;
    // Temporary storage
    valZ_eval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    gradU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    dKN        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fd, p);
    gradUgradN = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    // Build density-dependent conductivity function
    fe_dens_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valZ_eval)(i,j) = SIMPpower_*(one - minConductivity_)
          * std::pow((*valZ_eval)(i,j),SIMPpower_-one);
      }
    }
    // Build derivative of conductivity function
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*dKN,
                                                                *valZ_eval,
                                                                *(fe_dens_->N()));
    // Integrate stiffness term
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*gradUgradN,
                                                             *gradU_eval,
                                                             (*fe_vol_->gradNdetJ()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *gradUgradN,
                                                  *dKN,
                                                  Intrepid::COMP_CPP, false);

    // Apply Dirichlet boundary conditions
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < f; ++m) {
                  (*jac)(cidx,fidx_[j][l],m) = static_cast<Real>(0);
                }
              }
            }
          }
	}
      }
    }
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Poisson_TopOpt::Hessian_11): Zero Hessian.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real one(1);
    // Get dimensions
    const int c  = fe_vol_->gradN()->dimension(0);
    const int f  = fe_vol_->gradN()->dimension(1);
    const int fd = fe_dens_->gradN()->dimension(1);
    const int p  = fe_vol_->gradN()->dimension(2);
    const int d  = fe_vol_->gradN()->dimension(3);
    // Initialize Hessian storage
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fd, f);
    // Temporary storage
    ROL::Ptr<Intrepid::FieldContainer<Real>> valZ_eval, l0_coeff, dKN, gradL_eval, gradLgradN;
    valZ_eval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    l0_coeff   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    dKN        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fd, p);
    gradL_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    gradLgradN = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    // Build density-dependent conductivity function
    fe_dens_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valZ_eval)(i,j) = SIMPpower_*(one - minConductivity_)
          * std::pow((*valZ_eval)(i,j),SIMPpower_-one);
      }
    }
    // Apply Dirichlet conditions to the multipliers.
    *l0_coeff = *l_coeff;
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                (*l0_coeff)(cidx,fidx_[j][l]) = static_cast<Real>(0);
              }
            }
          }
        }
      }
    }
    // Build derivative of conductivity function
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*dKN,
                                                                *valZ_eval,
                                                                *(fe_dens_->N()));
    // Integrate stiffness term
    fe_vol_->evaluateGradient(gradL_eval, l0_coeff);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*gradLgradN,
                                                             *gradL_eval,
                                                             (*fe_vol_->gradNdetJ()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *dKN,
                                                  *gradLgradN,
                                                  Intrepid::COMP_CPP, false);
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real one(1);
    // Get dimensions
    const int c  = fe_vol_->gradN()->dimension(0);
    const int f  = fe_vol_->gradN()->dimension(1);
    const int fd = fe_dens_->gradN()->dimension(1);
    const int p  = fe_vol_->gradN()->dimension(2);
    const int d  = fe_vol_->gradN()->dimension(3);
    // Initialize Hessian storage
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, fd);
    // Temporary storage
    ROL::Ptr<Intrepid::FieldContainer<Real>> valZ_eval, l0_coeff, dKN, gradL_eval, gradLgradN;
    valZ_eval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    l0_coeff   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    dKN        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fd, p);
    gradL_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    gradLgradN = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    // Build density-dependent conductivity function
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valZ_eval)(i,j) = SIMPpower_*(one - minConductivity_)
          * std::pow((*valZ_eval)(i,j),SIMPpower_-one);
      }
    }
    // Apply Dirichlet conditions to the multipliers.
    *l0_coeff = *l_coeff;
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                (*l0_coeff)(cidx,fidx_[j][l]) = static_cast<Real>(0);
              }
            }
          }
        }
      }
    }
    // Build derivative of conductivity function
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*dKN,
                                                                *valZ_eval,
                                                                *(fe_dens_->N()));
    // Integrate stiffness term
    fe_vol_->evaluateGradient(gradL_eval, l0_coeff);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*gradLgradN,
                                                             *gradL_eval,
                                                             (*fe_vol_->gradNdetJ()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *gradLgradN,
                                                  *dKN,
                                                  Intrepid::COMP_CPP, false);
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real one(1), two(2);
    // Get dimensions
    const int c  = fe_vol_->gradN()->dimension(0);
    const int f  = fe_vol_->gradN()->dimension(1);
    const int fd = fe_dens_->gradN()->dimension(1);
    const int p  = fe_vol_->gradN()->dimension(2);
    const int d  = fe_vol_->gradN()->dimension(3);
    // Initialize Hessian storage
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fd, fd);
    // Temporary storage
    ROL::Ptr<Intrepid::FieldContainer<Real>> valZ_eval, l0_coeff, dKN, gradU_eval, gradL_eval, gradUgradL, NgradUgradL;
    valZ_eval   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    l0_coeff    = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    dKN         = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fd, p);
    gradU_eval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    gradL_eval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    gradUgradL  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    NgradUgradL = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fd, p);
    // Build density-dependent conductivity function
    fe_dens_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valZ_eval)(i,j) = SIMPpower_ * (SIMPpower_ - one)*(one - minConductivity_)
          * std::pow((*valZ_eval)(i,j),SIMPpower_-two);
      }
    }
    // Apply Dirichlet conditions to the multipliers.
    *l0_coeff = *l_coeff;
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                (*l0_coeff)(cidx,fidx_[j][l]) = static_cast<Real>(0);
              }
            }
          }
        }
      }
    }
    // Build derivative of conductivity function
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*dKN,
                                                                *valZ_eval,
                                                                *(fe_dens_->N()));
    // Integrate stiffness term
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    fe_vol_->evaluateGradient(gradL_eval, l0_coeff);
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(*gradUgradL,
                                                            *gradU_eval,
                                                            *gradL_eval);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*NgradUgradL,
                                                                *gradUgradL,
                                                                *(fe_dens_->NdetJ()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *dKN,
                                                  *NgradUgradL,
                                                  Intrepid::COMP_CPP, false);
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
    int c = fe_dens_->N()->dimension(0);
    int f = fe_dens_->N()->dimension(1);
    // INITIALIZE RIESZ
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    *riesz = *fe_dens_->massMat();
  }

  // This must be called before getFields2
  void setDensityFields(const std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> &basisPtrs) {
    TEUCHOS_TEST_FOR_EXCEPTION(getFields2called_, std::invalid_argument,
      ">>> PDE-OPT/topo-opt/elasticity/src/pde_elasticity.hpp: Must call before getFields2!");

    basisPtrDens_  = basisPtrs[0];
    basisPtrsDens_ = basisPtrs;
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields() {
    return basisPtrs_;
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields2() {
    getFields2called_ = true;
    return basisPtrsDens_;
  }

  void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real>> &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds) {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fe_vol_  = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtr_,cellCub_);
    fe_dens_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrDens_,cellCub_,false);
    fidx_ = fe_vol_->getBoundaryDofs();
    computeForce();
    computeDirichlet();
  }

  void computeForce(void) {
    int c = fe_vol_->cubPts()->dimension(0);
    int p = fe_vol_->cubPts()->dimension(1);
    int d = fe_vol_->cubPts()->dimension(2);
    std::vector<Real> pt(d,0);
    force_eval_  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    nforce_eval_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (*fe_vol_->cubPts())(i,j,k);
        }
        (*force_eval_)(i,j)  = ForceFunc(pt);
        (*nforce_eval_)(i,j) = -(*force_eval_)(i,j);
      }
    }
  }

  void computeDirichlet(void) {
    // Compute Dirichlet values at DOFs.
    int d = basisPtr_->getBaseCellTopology().getDimension();
    int numSidesets = bdryCellLocIds_.size();
    bdryCellDofValues_.resize(numSidesets);
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellDofValues_[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
        int f = basisPtr_->getCardinality();
        bdryCellDofValues_[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
        ROL::Ptr<Intrepid::FieldContainer<Real>> coords =
          ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, d);
        if (c > 0) {
          fe_vol_->computeDofCoords(coords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<f; ++l) {
            std::vector<Real> dofpoint(d);
            for (int m=0; m<d; ++m) {
              dofpoint[m] = (*coords)(k, l, m);
            }
            (*bdryCellDofValues_[i][j])(k, l) = dirichletFunc(dofpoint, i, j);
          }
        }
      }
    }
  }

  const ROL::Ptr<FE<Real>> getFE(void) const {
    return fe_vol_;
  }

  const ROL::Ptr<FE<Real>> getDensityFE(void) const {
    return fe_dens_;
  }

  const ROL::Ptr<Intrepid::FieldContainer<Real>> getForce(void) const {
    return force_eval_;
  }

}; // PDE_Poisson

#endif
