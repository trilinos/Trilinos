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

/*! \file  pde_poisson_topOpt.hpp
    \brief Implements the local PDE interface for the Poisson
           topology optimization problem.
*/

#ifndef PDE_POISSON_TOPOPT_HPP
#define PDE_POISSON_TOPOPT_HPP

#include "../../TOOLS/pde.hpp"
#include "../../TOOLS/fe.hpp"

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
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real> > cellCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real> > volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;
  std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real> > fe_vol_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int> > fidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellDofValues_;
  // Force function evaluated at cubature points
  ROL::Ptr<Intrepid::FieldContainer<Real> > force_eval_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > nforce_eval_;
  // Inputs
  Real minConductivity_;
  Real SIMPpower_;

  Real ForceFunc(const std::vector<Real> &x) const {
    return static_cast<Real>(0.01);
  }

  Real dirichletFunc(const std::vector<Real> & coords, const int sideset, const int locSideId) const {
    return static_cast<Real>(0);
  }

public:
  PDE_Poisson_TopOpt(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("Problem").get("Order of FE discretization",1);
    if (basisOrder == 1) {
      basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >>();
    }
    else if (basisOrder == 2) {
      basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >>();
    }
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();        // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                       // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree",2);     // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                       // create default cubature

    minConductivity_ = parlist.sublist("Problem").get("Minimum Conductivity",1.e-3);
    SIMPpower_       = parlist.sublist("Problem").get("SIMP Power",3.0);
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real> > & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    const Real one(1);
    // Get dimensions
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    // Initialize residual storage
    res = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Evaluate density at cubature points
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    // Build SIMP density at cubature points
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valZ_eval)(i,j) = minConductivity_
          + (one - minConductivity_) * std::pow((*valZ_eval)(i,j),SIMPpower_);
      }
    }
    // Build flux function
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    Intrepid::FieldContainer<Real> KgradU(c,p,d);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(KgradU,
                                                               *valZ_eval,
                                                               *gradU_eval);
    // Integrate stiffness term
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  KgradU,
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
        if ((i==2)) {
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

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    const Real one(1);
    // Get dimensions
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    // Initialize Jacobian storage
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // Build density-dependent conductivity function
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valZ_eval)(i,j) = minConductivity_
          + (one - minConductivity_) * std::pow((*valZ_eval)(i,j),SIMPpower_);
      }
    }
    // Build flux function
    Intrepid::FieldContainer<Real> KgradN(c,f,p,d);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(KgradN,
                                                                *valZ_eval,
                                                                *(fe_vol_->gradN()));
    // Integrate stiffness term
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  KgradN,
                                                  *(fe_vol_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);

    // Apply Dirichlet boundary conditions
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if ((i==2)) {
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

  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    const Real one(1);
    // Get dimensions
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    // Initialize Jacobian storage
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // Build density-dependent conductivity function
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valZ_eval)(i,j) = SIMPpower_*(one - minConductivity_)
          * std::pow((*valZ_eval)(i,j),SIMPpower_-one);
      }
    }
    // Build derivative of conductivity function
    Intrepid::FieldContainer<Real> dKN(c,f,p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(dKN,
                                                                *valZ_eval,
                                                                *(fe_vol_->N()));
    // Integrate stiffness term
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    Intrepid::FieldContainer<Real> gradUgradN(c,f,p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(gradUgradN,
                                                             *gradU_eval,
                                                             (*fe_vol_->gradNdetJ()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  gradUgradN,
                                                  dKN,
                                                  Intrepid::COMP_CPP, false);

    // Apply Dirichlet boundary conditions
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if ((i==2)) {
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

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Poisson_TopOpt::Hessian_11): Zero Hessian.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    const Real one(1);
    // Get dimensions
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    // Initialize Hessian storage
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // Build density-dependent conductivity function
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valZ_eval)(i,j) = SIMPpower_*(one - minConductivity_)
          * std::pow((*valZ_eval)(i,j),SIMPpower_-one);
      }
    }
    // Apply Dirichlet conditions to the multipliers.
    ROL::Ptr<Intrepid::FieldContainer<Real> > l0_coeff
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    *l0_coeff = *l_coeff;
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if ((i==2)) {
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
    Intrepid::FieldContainer<Real> dKN(c,f,p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(dKN,
                                                                *valZ_eval,
                                                                *(fe_vol_->N()));
    // Integrate stiffness term
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradL_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_vol_->evaluateGradient(gradL_eval, l0_coeff);
    Intrepid::FieldContainer<Real> gradLgradN(c,f,p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(gradLgradN,
                                                             *gradL_eval,
                                                             (*fe_vol_->gradNdetJ()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  dKN,
                                                  gradLgradN,
                                                  Intrepid::COMP_CPP, false);
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    const Real one(1);
    // Get dimensions
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    // Initialize Hessian storage
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // Build density-dependent conductivity function
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valZ_eval)(i,j) = SIMPpower_*(one - minConductivity_)
          * std::pow((*valZ_eval)(i,j),SIMPpower_-one);
      }
    }
    // Apply Dirichlet conditions to the multipliers.
    ROL::Ptr<Intrepid::FieldContainer<Real> > l0_coeff
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    *l0_coeff = *l_coeff;
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if ((i==2)) {
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
    Intrepid::FieldContainer<Real> dKN(c,f,p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(dKN,
                                                                *valZ_eval,
                                                                *(fe_vol_->N()));
    // Integrate stiffness term
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradL_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_vol_->evaluateGradient(gradL_eval, l0_coeff);
    Intrepid::FieldContainer<Real> gradLgradN(c,f,p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(gradLgradN,
                                                             *gradL_eval,
                                                             (*fe_vol_->gradNdetJ()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  gradLgradN,
                                                  dKN,
                                                  Intrepid::COMP_CPP, false);
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    const Real one(1), two(2);
    // Get dimensions
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    // Initialize Hessian storage
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // Build density-dependent conductivity function
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valZ_eval)(i,j) = SIMPpower_ * (SIMPpower_ - one)*(one - minConductivity_)
          * std::pow((*valZ_eval)(i,j),SIMPpower_-two);
      }
    }
    // Apply Dirichlet conditions to the multipliers.
    ROL::Ptr<Intrepid::FieldContainer<Real> > l0_coeff
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    *l0_coeff = *l_coeff;
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if ((i==2)) {
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
    Intrepid::FieldContainer<Real> dKN(c,f,p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(dKN,
                                                                *valZ_eval,
                                                                *(fe_vol_->N()));
    // Integrate stiffness term
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradL_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_vol_->evaluateGradient(gradL_eval, l0_coeff);
    Intrepid::FieldContainer<Real> gradUgradL(c,p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(gradUgradL,
                                                            *gradU_eval,
                                                            *gradL_eval);
    Intrepid::FieldContainer<Real> NgradUgradL(c,f,p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(NgradUgradL,
                                                                gradUgradL,
                                                                *(fe_vol_->NdetJ()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  dKN,
                                                  NgradUgradL,
                                                  Intrepid::COMP_CPP, false);
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    // GET DIMENSIONS
    int c = fe_vol_->N()->dimension(0);
    int f = fe_vol_->N()->dimension(1);
    // INITIALIZE RIESZ
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    *riesz = *fe_vol_->stiffMat();
    Intrepid::RealSpaceTools<Real>::add(*riesz,*(fe_vol_->massMat()));
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    // GET DIMENSIONS
    int c = fe_vol_->N()->dimension(0);
    int f = fe_vol_->N()->dimension(1);
    // INITIALIZE RIESZ
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    *riesz = *fe_vol_->massMat();
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > getFields() {
    return basisPtrs_;
  }

  void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real> > &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds) {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fe_vol_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtr_,cellCub_);
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
        ROL::Ptr<Intrepid::FieldContainer<Real> > coords =
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

  const ROL::Ptr<FE<Real> > getFE(void) const {
    return fe_vol_;
  }

  const ROL::Ptr<Intrepid::FieldContainer<Real> > getForce(void) const {
    return force_eval_;
  }

}; // PDE_Poisson



template <class Real>
class PDE_Filter : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real> > cellCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real> > volCellNodes_;
  // Finite element definition
  ROL::Ptr<FE<Real> > fe_;
  // Problem parameters.
  Real lengthScale_;

public:
  PDE_Filter(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >>();
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();            // get the cell type from the basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                           // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 2);        // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                           // create default cubature

    basisPtrs_.clear();
    basisPtrs_.push_back(basisPtr_);  // Filter components; there is only one, but we need d because of the infrastructure.

    // Other problem parameters.
    Real filterRadius = parlist.sublist("Problem").get("Filter Radius",  0.1);
    lengthScale_ = std::pow(filterRadius, 2)/12.0;
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real> > & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
 
    // Initialize residuals.
    res = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f);

    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real> > valU_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradU_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    fe_->evaluateValue(valZ_eval, z_coeff);
    fe_->evaluateGradient(gradU_eval, u_coeff);

    Intrepid::RealSpaceTools<Real>::scale(*gradU_eval, lengthScale_);
    Intrepid::RealSpaceTools<Real>::scale(*valZ_eval,  static_cast<Real>(-1));

    /*** Evaluate weak form of the residual. ***/
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *gradU_eval,           // R*gradU
                                                  *fe_->gradNdetJ(),     // gradN
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *valU_eval,            // U
                                                  *fe_->NdetJ(),         // N
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *valZ_eval,            // -Z
                                                  *fe_->NdetJ(),         // N
                                                  Intrepid::COMP_CPP,
                                                  true);
  }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
 
    // Initialize Jacobians.
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f);

    /*** Evaluate weak form of the Jacobian. ***/
    *jac = *(fe_->stiffMat());
    Intrepid::RealSpaceTools<Real>::scale(*jac, lengthScale_);    // ls*gradN1 . gradN2
    Intrepid::RealSpaceTools<Real>::add(*jac,*(fe_->massMat()));  // + N1 * N2
  }


  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
 
    // Initialize Jacobians.
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f);

    /*** Evaluate weak form of the Jacobian. ***/
    *jac = *(fe_->massMat());
    Intrepid::RealSpaceTools<Real>::scale(*jac, static_cast<Real>(-1));  // -N1 * N2
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Filter::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Filter::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Filter::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Filter::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    //throw Exception::NotImplemented(">>> (PDE_Filter::RieszMap_1): Not implemented.");
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
 
    // Initialize Riesz map.
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f);

    *riesz = *(fe_->stiffMat());
    Intrepid::RealSpaceTools<Real>::add(*riesz,*(fe_->massMat()));
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    //throw Exception::NotImplemented(">>> (PDE_Filter::RieszMap_2): Not implemented.");
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
 
    // Initialize Riesz map.
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f);

    *riesz = *(fe_->massMat());
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > getFields() {
    return basisPtrs_;
  }

  void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real> > &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds) {
    volCellNodes_ = volCellNodes;
    // Finite element definition.
    fe_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtr_,cellCub_);
  }

  void setFieldPattern(const std::vector<std::vector<int> > & fieldPattern) {}

}; // PDE_Filter

#endif
