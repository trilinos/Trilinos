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

/*! \file  pde.hpp
    \brief Implements the local PDE interface for the Poisson control problem.
*/

#ifndef PDE_ADV_DIFF_SUR_HPP
#define PDE_ADV_DIFF_SUR_HPP

#include "../../TOOLS/pde.hpp"
#include "../../TOOLS/fe.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C2_FEM.hpp"
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
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>> basisPtr2_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrs2_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real>> cellCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real>> volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real>> fe_vol_, fe_ctrl_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellDofValues_;

  int nx_, ny_;
  Real XL_, XU_, YL_, YU_;
  bool usePC_;

public:
  PDE_adv_diff(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("Problem").get("Order of FE Discretization",1);
    int cubDegree  = parlist.sublist("Problem").get("Cubature Degree",4);
    int probDim    = parlist.sublist("Problem").get("Problem Dimension",2);
    if (probDim == 2) {
      if (basisOrder == 1) {
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
      else if (basisOrder == 2) {
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
    }
    else if (probDim == 3) {
      if (basisOrder == 1) {
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_HEX_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
      else if (basisOrder == 2) {
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_HEX_C2_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
    }
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();        // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                       // create cubature factory
    cellCub_ = cubFactory.create(cellType, cubDegree);                       // create default cubature

    nx_ = parlist.sublist("Problem").get("Number of X-Cells", 4);
    ny_ = parlist.sublist("Problem").get("Number of Y-Cells", 2);
    XL_ = parlist.sublist("Geometry").get("X0", 0.0);
    YL_ = parlist.sublist("Geometry").get("Y0", 0.0);
    XU_ = XL_ + parlist.sublist("Geometry").get("Width",  2.0);
    YU_ = YL_ + parlist.sublist("Geometry").get("Height", 1.0);
    usePC_ = parlist.sublist("Problem").get("Piecewise Constant Controls", true);
    if (!usePC_) {
      if (probDim == 2) {
        basisPtr2_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
      else if (probDim == 3) {
        basisPtr2_ = ROL::makePtr<Intrepid::Basis_HGRAD_HEX_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
      basisPtrs2_.clear(); basisPtrs2_.push_back(basisPtr2_);
    }
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
    // INITIALIZE RESIDUAL AND STORAGE
    res = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    ROL::Ptr<Intrepid::FieldContainer<Real>> kappa, V;
    ROL::Ptr<Intrepid::FieldContainer<Real>> gradU_eval, kappa_gradU, V_gradU, valZ_eval;
    kappa       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    V           = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    gradU_eval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    kappa_gradU = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    V_gradU     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    valZ_eval   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    // COMPUTE PDE COEFFICIENTS
    computeCoefficients(kappa,V);
    // EVALUE GRADIENT OF U AT QUADRATURE POINTS
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    // COMPUTE DIFFUSION TERM
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(*kappa_gradU,
                                                               *kappa,
                                                               *gradU_eval);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *kappa_gradU,
                                                  *fe_vol_->gradNdetJ(),
                                                  Intrepid::COMP_CPP, false);
    // ADD ADVECTION TERM TO RESIDUAL
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(*V_gradU,
                                                            *V,
                                                            *gradU_eval);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *V_gradU,
                                                  *fe_vol_->NdetJ(),
                                                  Intrepid::COMP_CPP, true);

    // ADD CONTROL TERM TO RESIDUAL
    if (z_coeff != ROL::nullPtr) {
      fe_ctrl_->evaluateValue(valZ_eval, z_coeff);
      Intrepid::RealSpaceTools<Real>::scale(*valZ_eval,static_cast<Real>(-1));
      Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                    *valZ_eval,
                                                    *fe_vol_->NdetJ(),
                                                    Intrepid::COMP_CPP, true);
    }
    else {
      addControlOperator(valZ_eval,z_param);
      Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                    *valZ_eval,
                                                    *fe_vol_->NdetJ(),
                                                    Intrepid::COMP_CPP, true);
    }
    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
        int cidx = bdryCellLocIds_[0][j][k];
        for (int l = 0; l < numBdryDofs; ++l) {
          (*res)(cidx,fidx_[j][l])
            = (*u_coeff)(cidx,fidx_[j][l]) - (*bdryCellDofValues_[0][j])(k,fidx_[j][l]);
        }
      }
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
    // INITILAIZE JACOBIAN AND STORAGE
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    ROL::Ptr<Intrepid::FieldContainer<Real>> kappa, V, kappa_gradN, V_gradN;
    kappa       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    V           = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    kappa_gradN = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p, d);
    V_gradN     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    // COMPUTE PDE COEFFICIENTS
    computeCoefficients(kappa,V);
    // COMPUTE DIFFUSION TERM
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(*kappa_gradN,
                                                                *kappa,
                                                                *fe_vol_->gradN());
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *kappa_gradN,
                                                  *fe_vol_->gradNdetJ(),
                                                  Intrepid::COMP_CPP, false);
    // ADD ADVECTION TERM TO JACOBIAN
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*V_gradN,
                                                             *V,
                                                             *fe_vol_->gradN());
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *fe_vol_->NdetJ(),
                                                  *V_gradN,
                                                  Intrepid::COMP_CPP, true);
    // APPLY DIRICHLET CONDITIONS
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

  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_coeff != ROL::nullPtr) {
      // GET DIMENSIONS
      int c = fe_vol_->N()->dimension(0);
      int f1 = fe_vol_->N()->dimension(1);
      int f2 = fe_ctrl_->N()->dimension(1);
      // INITIALIZE RIESZ
      jac  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f1, f2);
      Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                    *fe_vol_->N(),
                                                    *fe_ctrl_->NdetJ(),
                                                    Intrepid::COMP_CPP, false);
      //*jac = *fe_vol_->massMat();
      Intrepid::RealSpaceTools<Real>::scale(*jac,static_cast<Real>(-1));
      // APPLY DIRICHLET CONDITIONS
      int numLocalSideIds = bdryCellLocIds_[0].size();
      for (int j = 0; j < numLocalSideIds; ++j) {
        int numCellsSide = bdryCellLocIds_[0][j].size();
        int numBdryDofs = fidx_[j].size();
        for (int k = 0; k < numCellsSide; ++k) {
          int cidx = bdryCellLocIds_[0][j][k];
          for (int l = 0; l < numBdryDofs; ++l) {
            //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
            for (int m = 0; m < f2; ++m) {
              (*jac)(cidx,fidx_[j][l],m) = static_cast<Real>(0);
            }
          }
        }
      }
    }
    else {
      throw Exception::Zero(">>> (PDE_stoch_adv_diff::Jacobian_2): Jacobian is zero.");
    }
  }

  void Jacobian_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_coeff != ROL::nullPtr) {
      throw Exception::Zero(">>> (PDE_stoch_adv_diff::Jacobian_3): Jacobian is zero.");
    }
    else {
      // GET DIMENSIONS
      int c = fe_vol_->gradN()->dimension(0);
      int f = fe_vol_->gradN()->dimension(1);
      int p = fe_vol_->gradN()->dimension(2);
      // ADD CONTROL TERM TO RESIDUAL
      ROL::Ptr<Intrepid::FieldContainer<Real>> ctrl
        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      ROL::Ptr<Intrepid::FieldContainer<Real>> B
        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      for (int i = 0; i < nx_; ++i) {
        for (int j = 0; j < ny_; ++j) {
          int ind = i + j*nx_;
          jac[ind] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
          addControlJaobian(B,i,j);
          Intrepid::FunctionSpaceTools::integrate<Real>(*jac[ind],
                                                        *B,
                                                        *fe_vol_->NdetJ(),
                                                        Intrepid::COMP_CPP, false);
          // APPLY DIRICHLET CONDITIONS
          int numLocalSideIds = bdryCellLocIds_[0].size();
          for (int k = 0; k < numLocalSideIds; ++k) {
            int numCellsSide = bdryCellLocIds_[0][k].size();
            int numBdryDofs = fidx_[k].size();
            for (int l = 0; l < numCellsSide; ++l) {
              int cidx = bdryCellLocIds_[0][k][l];
              for (int m = 0; m < numBdryDofs; ++m) {
                (*(jac[ind]))(cidx,fidx_[k][m]) = static_cast<Real>(0);
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
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_12): Hessian is zero.");
  }

  void Hessian_13(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_13): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_22): Hessian is zero.");
  }

  void Hessian_23(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_23): Hessian is zero.");
  }

  void Hessian_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_31): Hessian is zero.");
  }

  void Hessian_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_32): Hessian is zero.");
  }

  void Hessian_33(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_33): Hessian is zero.");
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
    int c = fe_ctrl_->N()->dimension(0);
    int f = fe_ctrl_->N()->dimension(1);
    // INITIALIZE RIESZ
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    *riesz = *fe_ctrl_->massMat();
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields() {
    return basisPtrs_;
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields2() {
    return basisPtrs2_;
  }

  void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real>> &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds) {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fe_vol_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtr_,cellCub_);
    if (!usePC_) {
      fe_ctrl_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtr2_,cellCub_);
    }
    // Set local boundary DOFs.
    fidx_ = fe_vol_->getBoundaryDofs();
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
            (*bdryCellDofValues_[i][j])(k, l) = evaluateDirichlet(dofpoint, i, j);
          }
        }
      }
    }
  }

  const ROL::Ptr<FE<Real>> getFE(void) const {
    return fe_vol_;
  }

  const ROL::Ptr<FE<Real>> getFE2(void) const {
    return fe_ctrl_;
  }

  void print(void) const {
    std::ofstream xfile, yfile;
    xfile.open("X.txt");
    yfile.open("Y.txt");
    for (int i = 0; i < nx_; ++i) {
      for (int j = 0; j < ny_; ++j) {
        xfile << i << std::endl;
        yfile << j << std::endl;
      }
    }
    xfile.close();
    yfile.close();
  }

private:

  Real evaluateDirichlet(const std::vector<Real> & coords, int sideset, int locSideId) const {
    return static_cast<Real>(0);
  }

  Real evaluateDiffusivity(const std::vector<Real> &x) const {
    return static_cast<Real>(0.01);
  }

  void evaluateVelocity(std::vector<Real> &adv, const std::vector<Real> &x) const {
    adv.assign(x.size(),static_cast<Real>(0));
    adv[0] = static_cast<Real>(1);
  }

  void computeCoefficients(ROL::Ptr<Intrepid::FieldContainer<Real>> &kappa,
                           ROL::Ptr<Intrepid::FieldContainer<Real>> &V) const {
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
      }
    }
  }

  void addControlOperator(ROL::Ptr<Intrepid::FieldContainer<Real>> &Bz,
                          const ROL::Ptr<const std::vector<Real>> &z) const {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    Real xl(0), xu(0), yl(0), yu(0);
    Bz->initialize();
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*fe_vol_->cubPts())(i,j,k);
        }
        for (int l = 0; l < nx_; ++l) {
          xl = XL_ + static_cast<Real>(l)*(XU_-XL_)/static_cast<Real>(nx_);
          xu = XL_ + static_cast<Real>(l+1)*(XU_-XL_)/static_cast<Real>(nx_);
          if ( pt[0] < xu && pt[0] >= xl ) {
            for (int m = 0; m < ny_; ++m) {
              int ind = l + m*nx_;
              yl = YL_ + static_cast<Real>(m)*(YU_-YL_)/static_cast<Real>(ny_);
              yu = YL_ + static_cast<Real>(m+1)*(YU_-YL_)/static_cast<Real>(ny_);
              if ( pt[1] < yu && pt[1] >= yl ) {
                (*Bz)(i,j) -= (*z)[ind];
              }
            }
          }
        }
      }
    }
  }

  void addControlJaobian(ROL::Ptr<Intrepid::FieldContainer<Real>> &B,
                         int l, int m) const {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    Real xl(0), xu(0), yl(0), yu(0);
    xl = XL_ + static_cast<Real>(l)*(XU_-XL_)/static_cast<Real>(nx_);
    xu = XL_ + static_cast<Real>(l+1)*(XU_-XL_)/static_cast<Real>(nx_);
    yl = YL_ + static_cast<Real>(m)*(YU_-YL_)/static_cast<Real>(ny_);
    yu = YL_ + static_cast<Real>(m+1)*(YU_-YL_)/static_cast<Real>(ny_);
    int ind = l + m*nx_;
    B->initialize();
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*fe_vol_->cubPts())(i,j,k);
        }
        if ( pt[0] < xu && pt[0] >= xl && pt[1] < yu && pt[1] >= yl ) {
          (*B)(i,j) = static_cast<Real>(-1);
        }
      }
    }
  }

}; // PDE_stoch_adv_diff

#endif
