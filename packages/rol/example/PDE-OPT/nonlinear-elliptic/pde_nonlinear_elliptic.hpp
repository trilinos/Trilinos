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

#ifndef PDE_NONLINEAR_ELLIPTIC_HPP
#define PDE_NONLINEAR_ELLIPTIC_HPP

#include "../TOOLS/pde.hpp"
#include "../TOOLS/fe.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real>
class PDE_Nonlinear_Elliptic : public PDE<Real> {
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
  ROL::Ptr<FE<Real> > fe_;

  void computeBetaGradU(const ROL::Ptr<Intrepid::FieldContainer<Real> > &BgradU,
                        const ROL::Ptr<const Intrepid::FieldContainer<Real> > &gradU) const {
    // GET DIMENSIONS
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    std::vector<Real> U(d), BU(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          U[k] = (*gradU)(i,j,k);
        }
        beta_value(BU,U);
        for (int k = 0; k < d; ++k) {
          (*BgradU)(i,j,k) = BU[k];
        }
      }
    }
  }

  void computeBetaJacobianGradU(const ROL::Ptr<Intrepid::FieldContainer<Real> > &BgradU,
                                const ROL::Ptr<const Intrepid::FieldContainer<Real> > &gradU) const {
    // GET DIMENSIONS
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    std::vector<Real> U(d);
    std::vector<std::vector<Real> > BU(d, U);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          U[k] = (*gradU)(i,j,k);
        }
        beta_jacobian(BU,U);
        for (int k = 0; k < d; ++k) {
          for (int m = 0; m < d; ++m) {
            (*BgradU)(i,j,k,m) = BU[k][m];
          }
        }
      }
    }
  }

  void computeBetaHessianGradUGradL(const ROL::Ptr<Intrepid::FieldContainer<Real> > &BgradUgradL,
                                    const ROL::Ptr<const Intrepid::FieldContainer<Real> > &gradU,
                                    const ROL::Ptr<const Intrepid::FieldContainer<Real> > &gradL) const {
    // GET DIMENSIONS
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    std::vector<Real> U(d), L(d);
    std::vector<std::vector<Real> > BU(d, U);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          U[k] = (*gradU)(i,j,k);
          L[k] = (*gradL)(i,j,k);
        }
        beta_hessian(BU,U,L);
        for (int k = 0; k < d; ++k) {
          for (int m = 0; m < d; ++m) {
            (*BgradUgradL)(i,j,k,m) = BU[k][m];
          }
        }
      }
    }
  }

  void computeDiffusivity(ROL::Ptr<Intrepid::FieldContainer<Real> > &kappa) const {
    // GET DIMENSIONS
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*fe_->cubPts())(i,j,k);
        }
        // Compute diffusivity kappa
        (*kappa)(i,j) = kappa_value(pt);
      }
    }
  }

  void computeRHS(ROL::Ptr<Intrepid::FieldContainer<Real> > &rhs) const {
    // GET DIMENSIONS
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*fe_->cubPts())(i,j,k);
        }
        // Compute diffusivity kappa
        (*rhs)(i,j) = -rhs_value(pt);
      }
    }
  }

public:
  PDE_Nonlinear_Elliptic(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("Problem").get("Basis Order",1);
    if (basisOrder == 1) {
      basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >>();
    }
    else if (basisOrder == 2) {
      basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >>();
    }
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();    // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                   // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree",2); // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                   // create default cubature
  }

  virtual void beta_value(std::vector<Real> &BU,
                          const std::vector<Real> &U) const {
    const int size = U.size();
    Real magU(0);
    for (int i = 0; i < size; ++i) {
      magU += U[i] * U[i];
    }
    for (int i = 0; i < size; ++i) {
      BU[i] = magU * U[i];
    }
  }

  virtual void beta_jacobian(std::vector<std::vector<Real> > &BU,
                             const std::vector<Real> &U) const {
    const int size = U.size();
    Real magU(0);
    for (int i = 0; i < size; ++i) {
      magU += U[i] * U[i];
    }
    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < size; ++j) {
        BU[i][j] = static_cast<Real>(2) * U[j] * U[i];
      }
      BU[i][i] += magU;
    }
  }

  virtual void beta_hessian(std::vector<std::vector<Real> > &BU,
                            const std::vector<Real> &U,
                            const std::vector<Real> &L) const {
    const int size = U.size();
    Real UL(0);
    for (int i = 0; i < size; ++i) {
      UL += static_cast<Real>(2) * U[i] * L[i];
    }
    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < size; ++j) {
        BU[i][j] = static_cast<Real>(2) * (L[j] * U[i] + L[i] * U[j]);
      }
      BU[i][i] += UL;
    }
  }

  Real kappa_value(const std::vector<Real> &x) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    const int ns = param.size(), d = x.size();
    if (ns > d+1) { 
      // random diffusion coefficient from i. babuska, f. nobile, r. tempone 2010.
      // simplified model for random stratified media.
      const Real one(1), two(2), three(3), eight(8), sixteen(16), half(0.5);
      const Real lc = one/sixteen, sqrtpi = std::sqrt(M_PI);
      const Real xi = std::sqrt(sqrtpi*lc), sqrt3 = std::sqrt(three);
      Real f(0), phi(0);                                          
      Real val = one + sqrt3*param[0]*std::sqrt(sqrtpi*lc*half);
      Real arg = one + sqrt3*std::sqrt(sqrtpi*lc*half);
      for (int i = 1; i < ns-(d+1); ++i) {
        f = floor(half*static_cast<Real>(i+1));     
        phi = ((i+1)%2 ? std::sin(f*M_PI*x[0]) : std::cos(f*M_PI*x[0]));
        val += xi*std::exp(-std::pow(f*M_PI*lc,two)/eight)*phi*sqrt3*param[i];
        arg += xi*std::exp(-std::pow(f*M_PI*lc,two)/eight)*std::abs(phi)*sqrt3;
      }
      return half + two*std::exp(val)/std::exp(arg);
    }
    return static_cast<Real>(1);
  }

  Real rhs_value(const std::vector<Real> &x) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    const int ns = param.size(), d = x.size();
    if (ns) {
      const int nk = ns - (d+1);
      const Real ten(10), two(2), pi(M_PI);
      Real val = std::pow(ten,two*param[nk]-two);
      for (int i = 0; i < d; ++i) {
        val *= std::sin(two*param[nk+i+1]*pi*x[i]);
      }
      return val;
    }
    return static_cast<Real>(0);
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real> > & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    // INITIALIZE RESIDUAL
    res = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // COMPUTE STIFFNESS TERM
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_->evaluateGradient(gradU_eval, u_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real> > kappa =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    computeDiffusivity(kappa);
    ROL::Ptr<Intrepid::FieldContainer<Real> > KgradU =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*KgradU,
                                                               *kappa,
                                                               *gradU_eval);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *KgradU,
                                                  *(fe_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // ADD MASS TERM
    ROL::Ptr<Intrepid::FieldContainer<Real> > valU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *valU_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // ADD NONLINEAR TERM
    ROL::Ptr<Intrepid::FieldContainer<Real> > BgradU
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    computeBetaGradU(BgradU,gradU_eval);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *BgradU,
                                                  *(fe_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // ADD RHS
    ROL::Ptr<Intrepid::FieldContainer<Real> > rhs =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    computeRHS(rhs);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *rhs,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // ADD CONTROL TERM
    if ( z_coeff != ROL::nullPtr ) {
      ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval =
        ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      fe_->evaluateValue(valZ_eval, z_coeff);
      Intrepid::RealSpaceTools<Real>::scale(*valZ_eval,static_cast<Real>(-1));
      Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                    *valZ_eval,
                                                    *(fe_->NdetJ()),
                                                    Intrepid::COMP_CPP, true);
    }
  }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    // INITIALIZE JACOBIAN
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // COMPUTE STIFFNESS TERM
    ROL::Ptr<Intrepid::FieldContainer<Real> > kappa =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    computeDiffusivity(kappa);
    ROL::Ptr<Intrepid::FieldContainer<Real> > KgradN =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p, d);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*KgradN,
                                                                *kappa,
                                                                *(fe_->gradN()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *KgradN,
                                                  *(fe_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // ADD MASS TERM
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *(fe_->N()),
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // ADD NONLINEAR TERM
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradU_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_->evaluateGradient(gradU_eval, u_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real> > BgradU
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d, d);
    computeBetaJacobianGradU(BgradU,gradU_eval);
    ROL::Ptr<Intrepid::FieldContainer<Real> > BgradUgradN
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(*BgradUgradN,
                                                                *BgradU,
                                                                *(fe_->gradN()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *BgradUgradN,
                                                  *(fe_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, true);
  }

  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    if ( z_coeff != ROL::nullPtr ) {
      // GET DIMENSIONS
      int c = fe_->N()->dimension(0);
      int f = fe_->N()->dimension(1);
      // INITIALIZE JACOBIAN
      jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
      // ADD CONTROL TERM
      Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                    *(fe_->N()),
                                                    *(fe_->NdetJ()),
                                                    Intrepid::COMP_CPP, false);
      Intrepid::RealSpaceTools<Real>::scale(*jac,static_cast<Real>(-1));
    }
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    // INITIALIZE JACOBIAN
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // COMPUTE NONLINEAR TERM
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradU_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_->evaluateGradient(gradU_eval, u_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradL_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_->evaluateGradient(gradL_eval, l_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real> > BgradUgradL
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d, d);
    computeBetaHessianGradUGradL(BgradUgradL,gradU_eval,gradL_eval);
    ROL::Ptr<Intrepid::FieldContainer<Real> > BgradUgradLgradN
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(*BgradUgradLgradN,
                                                                *BgradUgradL,
                                                                *(fe_->gradN()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *BgradUgradLgradN,
                                                  *(fe_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Nonlinear_Elliptic:Hessian_12: Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Nonlinear_Elliptic:Hessian_21: Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Nonlinear_Elliptic:Hessian_22: Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    // GET DIMENSIONS
    int c = fe_->N()->dimension(0);
    int f = fe_->N()->dimension(1);
    // INITIALIZE RIESZ
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    *riesz = *fe_->stiffMat();
    Intrepid::RealSpaceTools<Real>::add(*riesz,*(fe_->massMat()));
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    // GET DIMENSIONS
    int c = fe_->N()->dimension(0);
    int f = fe_->N()->dimension(1);
    // INITIALIZE RIESZ
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    *riesz = *fe_->massMat();
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
    fe_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtr_,cellCub_);
  }

  const ROL::Ptr<FE<Real> > getFE(void) const {
    return fe_;
  }

}; // PDE_Nonlinear_Elliptic

#endif
