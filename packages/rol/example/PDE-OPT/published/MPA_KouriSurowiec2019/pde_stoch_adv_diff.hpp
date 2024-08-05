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

#ifndef PDE_STOCH_ADV_DIFF_HPP
#define PDE_STOCH_ADV_DIFF_HPP

#include "../../TOOLS/pde.hpp"
#include "../../TOOLS/fe.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real>
class PDE_stoch_adv_diff : public PDE<Real> {
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

public:
  PDE_stoch_adv_diff(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("PDE Poisson").get("Basis Order",1);
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
    int cubDegree = parlist.sublist("PDE Poisson").get("Cubature Degree",2); // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                       // create default cubature
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real> > & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    // INITIALIZE RESIDUAL
    res = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // COMPUTE PDE COEFFICIENTS
    ROL::Ptr<Intrepid::FieldContainer<Real> > kappa
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > V
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    ROL::Ptr<Intrepid::FieldContainer<Real> > rhs
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    computeCoefficients(kappa,V,rhs);
    // COMPUTE DIFFUSION TERM
    // Compute grad(U)
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    // Multiply kappa * grad(U)
    Intrepid::FieldContainer<Real> kappa_gradU(c, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(kappa_gradU,
                                                               *kappa,
                                                               *gradU_eval);
    // Integrate (kappa * grad(U)) . grad(N)
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  kappa_gradU,
                                                  *(fe_vol_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // ADD ADVECTION TERM TO RESIDUAL
    // Multiply V . grad(U)
    Intrepid::FieldContainer<Real> V_gradU(c, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(V_gradU,
                                                            *V,
                                                            *gradU_eval);
    // Integrate (V . grad(U)) * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  V_gradU,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // ADD RHS TO RESIDUAL
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *rhs,
                                                  (*fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);

    // ADD CONTROL TERM TO RESIDUAL
    int size = z_param->size();
    ROL::Ptr<Intrepid::FieldContainer<Real> > ctrl
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < size; ++i) {
      computeControlOperator(ctrl,(*z_param)[i],i);
      Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                    *ctrl,
                                                    *(fe_vol_->NdetJ()),
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

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    // INITILAIZE JACOBIAN
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // COMPUTE PDE COEFFICIENTS
    ROL::Ptr<Intrepid::FieldContainer<Real> > kappa
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > V
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    ROL::Ptr<Intrepid::FieldContainer<Real> > rhs
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    computeCoefficients(kappa,V,rhs);
    // COMPUTE DIFFUSION TERM
    // Multiply kappa * grad(N)
    Intrepid::FieldContainer<Real> kappa_gradN(c, f, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(kappa_gradN,
                                                                *kappa,
                                                                *(fe_vol_->gradN()));
    // Integrate (kappa * grad(N)) . grad(N)
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  kappa_gradN,
                                                  *(fe_vol_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // ADD ADVECTION TERM TO JACOBIAN
    // Multiply V . grad(N)
    Intrepid::FieldContainer<Real> V_gradN(c, f, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(V_gradN,
                                                             *V,
                                                             *(fe_vol_->gradN()));
    // Integrate (V . grad(U)) * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *(fe_vol_->NdetJ()),
                                                  V_gradN,
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

  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Jacobian_2): Jacobian is zero.");
  }

  void Jacobian_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    // ADD CONTROL TERM TO RESIDUAL
    int size = z_param->size();
    ROL::Ptr<Intrepid::FieldContainer<Real> > ctrl
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < size; ++i) {
      jac[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      computeControlOperator(ctrl,static_cast<Real>(1),i);
      Intrepid::FunctionSpaceTools::integrate<Real>(*(jac[i]),
                                                    *ctrl,
                                                    *(fe_vol_->NdetJ()),
                                                    Intrepid::COMP_CPP, false);
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

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_12): Hessian is zero.");
  }

  void Hessian_13(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_13): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_22): Hessian is zero.");
  }

  void Hessian_23(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_23): Hessian is zero.");
  }

  void Hessian_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_31): Hessian is zero.");
  }

  void Hessian_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_32): Hessian is zero.");
  }

  void Hessian_33(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_stoch_adv_diff::Hessian_33): Hessian is zero.");
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
            (*bdryCellDofValues_[i][j])(k, l) = evaluateDirichlet(dofpoint, i, j);
          }
        }
      }
    }
  }

  const ROL::Ptr<FE<Real> > getFE(void) const {
    return fe_vol_;
  }

private:

  Real evaluateDirichlet(const std::vector<Real> & coords, int sideset, int locSideId) const {
    return static_cast<Real>(0);
  }

  Real evaluateDiffusivity(const std::vector<Real> &x) const {
    // random diffusion coefficient from i. babuska, f. nobile, r. tempone 2010.
    // simplified model for random stratified media.
    const int ns = 10;
    const Real one(1), two(2), three(3), eight(8), sixteen(16), half(0.5);
    const Real lc = one/sixteen, sqrtpi = std::sqrt(M_PI);
    const Real xi = std::sqrt(sqrtpi*lc), sqrt3 = std::sqrt(three);
    const std::vector<Real> param = PDE<Real>::getParameter();
    Real f(0), phi(0);                                          
    Real p25 = (param.size() > 25 ? param[25] : static_cast<Real>(0)); 
    Real val = one + sqrt3*p25*std::sqrt(sqrtpi*lc*half);
    Real arg = one + sqrt3*std::sqrt(sqrtpi*lc*half);
    for (int i = 1; i < ns; ++i) {
      f = floor(half*static_cast<Real>(i+1));     
      phi = ((i+1)%2 ? std::sin(f*M_PI*x[0]) : std::cos(f*M_PI*x[0]));
      Real pi25 = (static_cast<int>(param.size()) > i+25 ? param[i+25] : static_cast<Real>(0));
      val += xi*std::exp(-std::pow(f*M_PI*lc,two)/eight)*phi*sqrt3*pi25;
      arg += xi*std::exp(-std::pow(f*M_PI*lc,two)/eight)*std::abs(phi)*sqrt3;
    }
    return half + two*std::exp(val)/std::exp(arg);
  }

  void evaluateVelocity(std::vector<Real> &adv, const std::vector<Real> &x) const {
    const Real half(0.5), one(1), five(5), ten(10);
    const std::vector<Real> param = PDE<Real>::getParameter();
    Real p35 = (param.size() > 35 ? param[35] : static_cast<Real>(0));
    Real p36 = (param.size() > 36 ? param[36] : static_cast<Real>(0));
    const Real a = five*half*(p36+one);
    const Real b = five + (ten-five)*half*(p35+one);
    adv[0] = b - a*x[0];
    adv[1] =     a*x[1];
  }

  Real evaluateRHS(const std::vector<Real> &x) const {
    const int ns = 5;             
    const Real half(0.5), one(1), two(2);
    const std::vector<Real> param = PDE<Real>::getParameter();
    Real source(0), arg1(0), arg2(0), mag(0), x0(0), y0(0), sx(0), sy(0);
    // Upper and lower bounds on source magintudes
    const std::vector<Real> ml = {1.5, 1.2, 1.5, 1.2, 1.1};
    const std::vector<Real> mu = {2.5, 1.8, 1.9, 2.6, 1.5};
    // Upper and lower bounds on source locations
    const std::vector<Real> xl = {0.45, 0.75, 0.40, 0.05, 0.85};
    const std::vector<Real> xu = {0.55, 0.85, 0.60, 0.35, 0.95};
    const std::vector<Real> yl = {0.25, 0.55, 0.50, 0.45, 0.45};
    const std::vector<Real> yu = {0.35, 0.65, 0.70, 0.65, 0.55};
    // Upper and lower bounds on source widths
    const std::vector<Real> sxl = {0.03, 0.02, 0.01, 0.02, 0.015};
    const std::vector<Real> sxu = {0.07, 0.04, 0.05, 0.04, 0.025};
    const std::vector<Real> syl = {0.04, 0.01, 0.02, 0.02, 0.01};
    const std::vector<Real> syu = {0.12, 0.05, 0.04, 0.04, 0.03};
    for (int i=0; i<ns; ++i) {
      Real pi  = (static_cast<int>(param.size()) > i      ? param[i]      : static_cast<Real>(0));
      Real pi1 = (static_cast<int>(param.size()) > i+1*ns ? param[i+1*ns] : static_cast<Real>(0));
      Real pi2 = (static_cast<int>(param.size()) > i+2*ns ? param[i+2*ns] : static_cast<Real>(0));
      Real pi3 = (static_cast<int>(param.size()) > i+3*ns ? param[i+3*ns] : static_cast<Real>(0));
      Real pi4 = (static_cast<int>(param.size()) > i+4*ns ? param[i+4*ns] : static_cast<Real>(0));

      mag  = ml[i] + (mu[i]-ml[i])*half*(pi+one);
      x0   = xl[i] + (xu[i]-xl[i])*half*(pi1+one);
      y0   = yl[i] + (yu[i]-yl[i])*half*(pi3+one);
      sx   = sxl[i] + (sxu[i]-sxl[i])*half*(pi2+one);
      sy   = syl[i] + (syu[i]-syl[i])*half*(pi4+one);
      arg1 = std::pow((x[0]-x0)/sx, two);
      arg2 = std::pow((x[1]-y0)/sy, two);
      source += mag*std::exp(-half*(arg1+arg2));
    }
    return source;
  }

  void computeCoefficients(ROL::Ptr<Intrepid::FieldContainer<Real> > &kappa,
                           ROL::Ptr<Intrepid::FieldContainer<Real> > &V,
                           ROL::Ptr<Intrepid::FieldContainer<Real> > &rhs) const {
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
        // Compute forcing term f
        (*rhs)(i,j) = -evaluateRHS(pt);
      }
    }
  }

  Real evaluateControlOperator(const std::vector<Real> &x, const int i) const {
    const Real sx(0.05), sy(0.05), half(0.5);
    const std::vector<Real> xl = {0.25, 0.50, 0.75, 0.25, 0.50, 0.75, 0.25, 0.50, 0.75};
    const std::vector<Real> yl = {0.25, 0.25, 0.25, 0.50, 0.50, 0.50, 0.75, 0.75, 0.75};
    return -std::exp(- half*(x[0]-xl[i])*(x[0]-xl[i]) / (sx*sx)
                     - half*(x[1]-yl[i])*(x[1]-yl[i]) / (sy*sy));
  }
  
  void computeControlOperator(ROL::Ptr<Intrepid::FieldContainer<Real> > &ctrl,
                              const Real z, const int I) const {
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
        // Compute control operator
        (*ctrl)(i,j) = -z*evaluateControlOperator(pt,I);
      }
    }
  }

}; // PDE_stoch_adv_diff

#endif
