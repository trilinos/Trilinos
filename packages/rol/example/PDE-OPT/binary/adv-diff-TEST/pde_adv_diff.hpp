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
#include "hilbert.hpp"

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
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real>> cellCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real>> volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real>> fe_vol_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellDofValues_;

  int order_;
  Real XL_, XU_, YL_, YU_;
  std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> ctrl_wts_;
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

    order_ = parlist.sublist("Problem").get("Hilbert Curve Order", 2);
    XL_ = parlist.sublist("Geometry").get("X0", 0.0);
    YL_ = parlist.sublist("Geometry").get("Y0", 0.0);
    XU_ = XL_ + parlist.sublist("Geometry").get("Width",  2.0);
    YU_ = YL_ + parlist.sublist("Geometry").get("Height", 1.0);
    usePC_ = parlist.sublist("Problem").get("Piecewise Constant Controls", true);
  }

//  PDE_adv_diff(const PDE_adv_diff &rpde)
//    : basisPtr_(rpde.basisPtr_), basisPtrs_(rpde.basisPtrs_),
//      cellCub_(rpde.cellCub_), volCellNodes_(rpde.volCellNodes_),
//      bdryCellNodes_(rpde.bdryCellNodes_),
//      bdryCellLocIds_(rpde.bdryCellLocIds_),
//      fe_vol_(rpde.fe_vol_), fidx_(rpde.fidx_),
//      bdryCellDofValues_(rpde.bdryCellDofValues_),
//      order_(rpde.order_), XL_(rpde.XL_), XU_(rpde.XU_), YL(rpde.YL_),
//      YU_(rpde.YU_), ctrl_wts_(rpde.ctrl_wts_), usePC_(rpde.usePC_) {}

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
      fe_vol_->evaluateValue(valZ_eval, z_coeff);
      Intrepid::RealSpaceTools<Real>::scale(*valZ_eval,static_cast<Real>(-1));
      Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                    *valZ_eval,
                                                    *fe_vol_->NdetJ(),
                                                    Intrepid::COMP_CPP, true);
    }
    else {
      int n = std::pow(2,order_);
      for (int i = 0; i < n*n; ++i) {
        *valZ_eval = *ctrl_wts_[i];
        Intrepid::RealSpaceTools<Real>::scale(*valZ_eval,-(*z_param)[i]);
        Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                      *valZ_eval,
                                                      *fe_vol_->NdetJ(),
                                                      Intrepid::COMP_CPP, true);
      }
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
      int f = fe_vol_->N()->dimension(1);
      // INITIALIZE RIESZ
      jac  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
      *jac = *fe_vol_->massMat();
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
            for (int m = 0; m < f; ++m) {
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
      int n = std::pow(2,order_);
      ROL::Ptr<Intrepid::FieldContainer<Real>> ctrl
        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      for (int i = 0; i < n*n; ++i) {
        jac[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
        Intrepid::FunctionSpaceTools::integrate<Real>(*jac[i],
                                                      *ctrl_wts_[i],
                                                      *fe_vol_->NdetJ(),
                                                      Intrepid::COMP_CPP, false);
        Intrepid::RealSpaceTools<Real>::scale(*jac[i],static_cast<Real>(-1));
        // APPLY DIRICHLET CONDITIONS
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
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fe_vol_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtr_,cellCub_);
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
    if (usePC_) {
      computeControlWeights();
    }
  }

  const ROL::Ptr<FE<Real>> getFE(void) const {
    return fe_vol_;
  }

  void print(void) const {
    int n = std::pow(2,order_);
    int x(0), y(0);
    std::ofstream xfile, yfile;
    xfile.open("X.txt");
    yfile.open("Y.txt");
    for (int i = 0; i < n*n; ++i) {
      hilbert::d2xy(order_, i, x, y);
      xfile << x << std::endl;
      yfile << y << std::endl;
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

  void computeControlWeights(void) {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    Real xl(0), xu(0), yl(0), yu(0);
    int n = std::pow(2,order_), x(0), y(0), D(0);
    ctrl_wts_.clear(); ctrl_wts_.resize(n*n);
    for (int l = 0; l < n*n; ++l) {
      ctrl_wts_[l] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
      ctrl_wts_[l]->initialize();
    }
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*fe_vol_->cubPts())(i,j,k);
        }
        for (int l = 0; l < n; ++l) {
          for (int m = 0; m < n; ++m) {
            hilbert::xy2d(order_,l,m,D);
            D = l + m*n;
            xl = XL_ + static_cast<Real>(l)*(XU_-XL_)/static_cast<Real>(n);
            xu = XL_ + static_cast<Real>(l+1)*(XU_-XL_)/static_cast<Real>(n);
            yl = YL_ + static_cast<Real>(m)*(YU_-YL_)/static_cast<Real>(n);
            yu = YL_ + static_cast<Real>(m+1)*(YU_-YL_)/static_cast<Real>(n);
            if ( (pt[0] < xu && pt[0] >= xl) && (pt[1] < yu && pt[1] >= yl) ) {
              (*ctrl_wts_[D])(i,j) = static_cast<Real>(1);
            }
          }
        }
      }
    }
  }

}; // PDE_stoch_adv_diff

#endif
