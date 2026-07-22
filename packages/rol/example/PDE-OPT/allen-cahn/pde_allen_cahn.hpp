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

#ifndef PDE_ALLEN_CAHN_HPP
#define PDE_ALLEN_CAHN_HPP

#include "../TOOLS/pde.hpp"
#include "../TOOLS/fe.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real>
class PDE_Allen_Cahn : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>> basisPtr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real>> cellCub_;
  ROL::Ptr<Intrepid::Cubature<Real>> bdryCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real>> volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real>> fe_vol_;
  std::vector<std::vector<ROL::Ptr<FE<Real>>>> fe_bdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fidx_;

  Real uScale_, vScale_;
  Real robinCoeff_;

public:
  PDE_Allen_Cahn(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("Problem").get("Order of FE discretization",1);
    if (basisOrder == 1) {
      basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
    }
    else if (basisOrder == 2) {
      basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real>>>();
    }
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();                  // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                                 // create cubature factory
    int cubDegree = parlist.sublist("PDE Poisson Boltzmann").get("Cubature Degree",2); // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                                 // create default cubature

    int d = cellType.getDimension();
    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2); // set cubature degree, e.g., 2
    bdryCub_ = cubFactory.create(bdryCellType, bdryCubDegree);

    uScale_     = parlist.sublist("Problem").get("Scale Third Order Term",1.0);
    vScale_     = parlist.sublist("Problem").get("Scale First Order Term",1.0);
    robinCoeff_ = parlist.sublist("Problem").get("Robin Coefficient",1e0);
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real>> & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN()->dimension(0);
    const int f = fe_vol_->gradN()->dimension(1);
    const int p = fe_vol_->gradN()->dimension(2);
    const int d = fe_vol_->gradN()->dimension(3);
    // INITIALIZE RESIDUAL
    res = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // STORAGE
    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval, gradU_eval, lambda, lambda_gradU_eval, phi_valU_eval;
    valU_eval         = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    gradU_eval        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    lambda            = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    lambda_gradU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    phi_valU_eval     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    // COMPUTE PDE COEFFICIENTS
    computeCoefficients(lambda);
    // EVALUATE STATE AT QUADRATURE POINTS
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    // ADD STIFFNESS TERM TO RESIDUAL
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(*lambda_gradU_eval,
                                                               *lambda,
                                                               *gradU_eval);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *lambda_gradU_eval,
                                                  *(fe_vol_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // ADD NONLINEAR TERM TO RESIDUAL
    computeNonlinearity(phi_valU_eval, valU_eval,0);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *phi_valU_eval,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // APPLY ROBIN CONDITIONS
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellDofValues;
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      int numLocalSideIds = bdryCellLocIds_[i].size();
      for (int j = 0; j < numLocalSideIds; ++j) {
        int numCellsSide = bdryCellLocIds_[i][j].size();
        const int numCubPerSide = bdryCub_->getNumPoints();
        if (numCellsSide) {
          // Get U and Z coefficients on Robin boundary
          ROL::Ptr<Intrepid::FieldContainer<Real >> u_coeff_bdry, z_coeff_bdry;
          u_coeff_bdry = getBoundaryCoeff(*u_coeff, i, j);
          z_coeff_bdry = getBoundaryCoeff(*z_coeff, i, j);
          // Evaluate U and Z on FE basis
          ROL::Ptr<Intrepid::FieldContainer<Real >> valU_eval_bdry, valZ_eval_bdry;
          valU_eval_bdry = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          valZ_eval_bdry = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          fe_bdry_[i][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
          fe_bdry_[i][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
          // Compute Robin residual
          Intrepid::FieldContainer<Real> robinRes(numCellsSide, f);
          ROL::Ptr< Intrepid::FieldContainer<Real>> robinVal
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          computeRobin(robinVal,valU_eval_bdry,valZ_eval_bdry,i,j,0);
          Intrepid::FunctionSpaceTools::integrate<Real>(robinRes,
                                                        *robinVal,
                                                        *(fe_bdry_[i][j]->NdetJ()),
                                                        Intrepid::COMP_CPP, false);
          // Add Robin residual to volume residual
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
            for (int l = 0; l < f; ++l) { 
              (*res)(cidx,l) += robinRes(k,l);
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
    // GET DIMENSIONS
    const int c = fe_vol_->gradN()->dimension(0);
    const int f = fe_vol_->gradN()->dimension(1);
    const int p = fe_vol_->gradN()->dimension(2);
    const int d = fe_vol_->gradN()->dimension(3);
    // INITIALIZE JACOBIAN
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // STORAGE
    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval, lambda, lambda_gradN_eval, dphi_valU_eval, NdphiU;
    valU_eval         = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    lambda            = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    lambda_gradN_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p, d);
    dphi_valU_eval    = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    NdphiU            = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    // COMPUTE PDE COEFFICIENTS
    computeCoefficients(lambda);
    // EVALUATE STATE AT QUADRATURE POINTS
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    // ADD STIFFNESS TERM TO JACOBIAN
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(*lambda_gradN_eval,
                                                                *lambda,
                                                                *(fe_vol_->gradN()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *lambda_gradN_eval,
                                                  *(fe_vol_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // ADD NONLINEAR TERM
    computeNonlinearity(dphi_valU_eval, valU_eval, 1);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*NdphiU,
                                                                *dphi_valU_eval,
                                                                *(fe_vol_->N()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *NdphiU,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // APPLY ROBIN CONDITIONS
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      int numLocalSideIds = bdryCellLocIds_[i].size();
      for (int j = 0; j < numLocalSideIds; ++j) {
        int numCellsSide = bdryCellLocIds_[i][j].size();
        const int numCubPerSide = bdryCub_->getNumPoints();
        if (numCellsSide) {
          // Get U and Z coefficients on Robin boundary
          ROL::Ptr<Intrepid::FieldContainer<Real >> u_coeff_bdry, z_coeff_bdry;
          u_coeff_bdry = getBoundaryCoeff(*u_coeff, i, j);
          z_coeff_bdry = getBoundaryCoeff(*z_coeff, i, j);
          // Evaluate U and Z on FE basis
          ROL::Ptr<Intrepid::FieldContainer<Real >> valU_eval_bdry, valZ_eval_bdry;
          valU_eval_bdry = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          valZ_eval_bdry = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          fe_bdry_[i][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
          fe_bdry_[i][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
          // Compute Robin residual
          Intrepid::FieldContainer<Real> robinVal_N(numCellsSide, f, numCubPerSide);
          Intrepid::FieldContainer<Real> robinJac(numCellsSide, f, f);
          ROL::Ptr< Intrepid::FieldContainer<Real>> robinVal
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          computeRobin(robinVal,valU_eval_bdry,valZ_eval_bdry,i,j,1,1);
          Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(robinVal_N,
                                                                      *robinVal,
                                                                      *(fe_bdry_[i][j]->N()));
          Intrepid::FunctionSpaceTools::integrate<Real>(robinJac,
                                                        robinVal_N,
                                                        *(fe_bdry_[i][j]->NdetJ()),
                                                        Intrepid::COMP_CPP, false);
          // Add Robin residual to volume residual
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
            for (int l = 0; l < f; ++l) { 
              for (int m = 0; m < f; ++m) { 
                (*jac)(cidx,l,m) += robinJac(k,l,m);
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
    // GET DIMENSIONS
    const int c = fe_vol_->gradN()->dimension(0);
    const int f = fe_vol_->gradN()->dimension(1);
    // INITIALIZE JACOBIAN
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // APPLY ROBIN CONDITIONS
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      int numLocalSideIds = bdryCellLocIds_[i].size();
      for (int j = 0; j < numLocalSideIds; ++j) {
        int numCellsSide = bdryCellLocIds_[i][j].size();
        const int numCubPerSide = bdryCub_->getNumPoints();
        if (numCellsSide) {
          // Get U and Z coefficients on Robin boundary
          ROL::Ptr<Intrepid::FieldContainer<Real >> u_coeff_bdry, z_coeff_bdry;
          u_coeff_bdry = getBoundaryCoeff(*u_coeff, i, j);
          z_coeff_bdry = getBoundaryCoeff(*z_coeff, i, j);
          // Evaluate U and Z on FE basis
          ROL::Ptr<Intrepid::FieldContainer<Real >> valU_eval_bdry, valZ_eval_bdry;
          valU_eval_bdry = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          valZ_eval_bdry = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          fe_bdry_[i][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
          fe_bdry_[i][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
          // Compute Robin residual
          Intrepid::FieldContainer<Real> robinVal_N(numCellsSide, f, numCubPerSide);
          Intrepid::FieldContainer<Real> robinJac(numCellsSide, f, f);
          ROL::Ptr< Intrepid::FieldContainer<Real>> robinVal
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          computeRobin(robinVal,valU_eval_bdry,valZ_eval_bdry,i,j,1,2);
          Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(robinVal_N,
                                                                      *robinVal,
                                                                      *(fe_bdry_[i][j]->N()));
          Intrepid::FunctionSpaceTools::integrate<Real>(robinJac,
                                                        robinVal_N,
                                                        *(fe_bdry_[i][j]->NdetJ()),
                                                        Intrepid::COMP_CPP, false);
          // Add Robin residual to volume residual
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
            for (int l = 0; l < f; ++l) { 
              for (int m = 0; m < f; ++m) { 
                (*jac)(cidx,l,m) += robinJac(k,l,m);
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
    // GET DIMENSIONS
    const int c = fe_vol_->gradN()->dimension(0);
    const int f = fe_vol_->gradN()->dimension(1);
    const int p = fe_vol_->gradN()->dimension(2);
    // INITIALIZE HESSIAN
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // STORAGE
    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval, valL_eval, d2phi_valU_eval, Ld2phiU, NLd2phiU;
    valU_eval       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    valL_eval       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    d2phi_valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Ld2phiU         = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    NLd2phiU        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    // EVALUATE STATE AT QUADRATURE POINTS
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    fe_vol_->evaluateValue(valL_eval, l_coeff);
    // COMPUTE NONLINEAR TERM
    computeNonlinearity(d2phi_valU_eval, valU_eval, 2);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*Ld2phiU,
                                                               *d2phi_valU_eval,
                                                               *valL_eval);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*NLd2phiU,
                                                                *Ld2phiU,
                                                                *(fe_vol_->N()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *NLd2phiU,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Allen_Cahn:Hessian_12: Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Allen_Cahn:Hessian_21: Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Allen_Cahn:Hessian_22: Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    // GET DIMENSIONS
    const int c = fe_vol_->N()->dimension(0);
    const int f = fe_vol_->N()->dimension(1);
    // INITIALIZE RIESZ
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    *riesz = *fe_vol_->stiffMat();
    Intrepid::RealSpaceTools<Real>::add(*riesz,*(fe_vol_->massMat()));
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    // GET DIMENSIONS
    const int c = fe_vol_->N()->dimension(0);
    const int f = fe_vol_->N()->dimension(1);
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
    // Construct boundary FEs
    const int numSidesets = bdryCellNodes.size();
    fe_bdry_.resize(numSidesets);
    for(int i = 0; i < numSidesets; ++i) {
      int numLocSides = bdryCellNodes[i].size();
      fe_bdry_[i].resize(numLocSides);
      for (int j = 0; j < numLocSides; ++j) {
        if (bdryCellNodes_[i][j] != ROL::nullPtr) {
          fe_bdry_[i][j] = ROL::makePtr<FE<Real>>(bdryCellNodes_[i][j],basisPtr_,bdryCub_,j);
        }
      }
    }
  }

  const ROL::Ptr<FE<Real>> getFE(void) const {
    return fe_vol_;
  }

  const std::vector<std::vector<ROL::Ptr<FE<Real>>>> getBdryFE(void) const {
    return fe_bdry_;
  }

  const ROL::Ptr<Intrepid::FieldContainer<Real>> getCellNodes(void) const {
    return volCellNodes_;
  }

  const std::vector<std::vector<std::vector<int>>> getBdryCellLocIds(void) const {
    return bdryCellLocIds_;
  }

private:

  Real evaluateLambda(const std::vector<Real> &x) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    Real val(1);
    if (param.size()) {
      val += param[0];
    }
    return val;
  }

  void computeCoefficients(ROL::Ptr<Intrepid::FieldContainer<Real>> &lambda) const {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN()->dimension(0);
    const int p = fe_vol_->gradN()->dimension(2);
    const int d = fe_vol_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*fe_vol_->cubPts())(i,j,k);
        }
        (*lambda)(i,j) = evaluateLambda(pt);
      }
    }
  }

  Real evaluateDelta(const std::vector<Real> &x) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    Real val(1);
    if (param.size()) {
      val += param[1];
    }
    return val;
  }

  Real evaluateNonlinearity(const std::vector<Real> &x, const Real u, const int deriv) const {
    Real val(0);
    if (deriv == 0) {
      val = uScale_*std::pow(u,3) - vScale_*u;
    }
    if (deriv == 1) {
      val = uScale_*static_cast<Real>(3)*std::pow(u,2) - vScale_;
    }
    if (deriv == 2) {
      val = uScale_*static_cast<Real>(6)*u;
    }
    return evaluateDelta(x) * val;
  }

  void computeNonlinearity(ROL::Ptr<Intrepid::FieldContainer<Real>> &val,
                           const ROL::Ptr<Intrepid::FieldContainer<Real>> &u,
                           const int deriv = 0) const {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN()->dimension(0);
    const int p = fe_vol_->gradN()->dimension(2);
    const int d = fe_vol_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*fe_vol_->cubPts())(i,j,k);
        }
        (*val)(i,j) = evaluateNonlinearity(pt, (*u)(i,j), deriv);
      }
    }
  }

  Real evaluateRobin(const Real u, const Real z, const std::vector<Real> &x,
                     const int sideset, const int locSideId,
                     const int deriv = 0, const int component = 1) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    Real h(robinCoeff_);
    if (param.size()) {
      h += param[2];
    }
    if ( deriv == 1 ) {
      return (component==1) ? h : -h;
    }
    if ( deriv > 1 ) {
      return static_cast<Real>(0);
    }
    return h * (u - z);
  }

  void computeRobin(ROL::Ptr<Intrepid::FieldContainer<Real>> &robin,
                    const ROL::Ptr<Intrepid::FieldContainer<Real>> &u,
                    const ROL::Ptr<Intrepid::FieldContainer<Real>> &z,
                    const int sideset,
                    const int locSideId,
                    const int deriv = 0,
                    const int component = 1) const {
    const int c = u->dimension(0);
    const int p = u->dimension(1);
    const int d = fe_vol_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (*fe_bdry_[sideset][locSideId]->cubPts())(i,j,k);
        }
        (*robin)(i,j) = evaluateRobin((*u)(i,j),(*z)(i,j),pt,sideset,locSideId,deriv,component);
      }
    }
  }

  ROL::Ptr<Intrepid::FieldContainer<Real>> getBoundaryCoeff(
      const Intrepid::FieldContainer<Real> & cell_coeff,
      int sideSet, int cell) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][cell];
    const int numCellsSide = bdryCellLocId.size();
    const int f = basisPtr_->getCardinality();
    
    ROL::Ptr<Intrepid::FieldContainer<Real >> bdry_coeff = 
      ROL::makePtr<Intrepid::FieldContainer<Real >>(numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        (*bdry_coeff)(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

}; // PDE_Allen_Cahn

#endif
