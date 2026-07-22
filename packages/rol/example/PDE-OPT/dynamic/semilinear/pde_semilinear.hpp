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

#ifndef PDE_PARABOLIC_SEMILINEAR_HPP
#define PDE_PARABOLIC_SEMILINEAR_HPP

#include "../../TOOLS/pde.hpp"
#include "../../TOOLS/fe.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real>
class PDE_Semilinear : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>> basisPtr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real>> cellCub_;
  // Finite element definition
  ROL::Ptr<FE<Real>> fe_vol_;

  int type_;

public:
  PDE_Semilinear(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("Semilinear").get("Order of FE discretization",1);
    if (basisOrder == 1) {
      basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
    }
    else if (basisOrder == 2) {
      basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real>>>();
    }
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();                    // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                                   // create cubature factory
    int cubDegree = parlist.sublist("Semilinear").get("Cubature Degree",2);              // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                                   // create default cubature

    type_ = parlist.sublist("Semilinear").get("Nonlinearity Type",0);
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
    // COMPUTE PDE COEFFICIENTS
    Intrepid::FieldContainer<Real> diff(c, p), adv(c, p, d), rhs(c, p);
    computeCoefficients(diff,adv,rhs);
    // Evaluate PDE solution
    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real>> gradU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    // Evaluate control
    ROL::Ptr<Intrepid::FieldContainer<Real>> valZ_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    // Multiply diff * grad(U)
    Intrepid::FieldContainer<Real> diff_gradU(c, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(diff_gradU,
                                                               diff,
                                                               *gradU_eval);
    // Integrate (diff * grad(U)) . grad(N)
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  diff_gradU,
                                                  *(fe_vol_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // Multiply adv . grad(U)
    Intrepid::FieldContainer<Real> adv_gradU(c, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(adv_gradU,
                                                            adv,
                                                            *gradU_eval);
    // Integrate (adv . grad(U)) * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  adv_gradU,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // Integrate L * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  rhs,
                                                  (*fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // Compute F(U)
    Intrepid::FieldContainer<Real> phi_valU_eval(c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        phi_valU_eval(i,j) = evaluateReaction((*valU_eval)(i,j),0);
      }
    }
    // Integrate F(U) * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  phi_valU_eval,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // Integrate -Z *N
    Intrepid::FieldContainer<Real> valZ_scal(c, p);
    Intrepid::RealSpaceTools<Real>::scale(valZ_scal,*valZ_eval,static_cast<Real>(-1));
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  valZ_scal,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
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
    // COMPUTE PDE COEFFICIENTS
    Intrepid::FieldContainer<Real> diff(c, p), adv(c, p, d), rhs(c, p);
    computeCoefficients(diff,adv,rhs);
    // Evaluate PDE solution
    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    // Evaluate control
    ROL::Ptr<Intrepid::FieldContainer<Real>> valZ_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    // Multiply diff * grad(N)
    Intrepid::FieldContainer<Real> diff_gradN(c, f, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(diff_gradN,
                                                                diff,
                                                                *(fe_vol_->gradN()));
    // Integrate (diff * grad(N)) . grad(N)
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  diff_gradN,
                                                  *(fe_vol_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // Multiply adv . grad(N)
    Intrepid::FieldContainer<Real> adv_gradN(c, f, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(adv_gradN,
                                                             adv,
                                                             *(fe_vol_->gradN()));
    // Integrate (adv . grad(N)) * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  adv_gradN,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // Compute F'(U) * N
    Intrepid::FieldContainer<Real> dphi_valU_eval(c, p);
    Intrepid::FieldContainer<Real> NdphiU(c, f, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        dphi_valU_eval(i,j) = evaluateReaction((*valU_eval)(i,j),1);
      }
    }
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(NdphiU,
                                                                dphi_valU_eval,
                                                                *(fe_vol_->N()));
    // Integrate F'(U) * N * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  NdphiU,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
  }

  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN()->dimension(0);
    const int f = fe_vol_->gradN()->dimension(1);
    const int p = fe_vol_->gradN()->dimension(2);
    // INITIALIZE JACOBIAN
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // ADD CONTROL TERM
    Intrepid::FieldContainer<Real> valN_scal(c, f, p);
    Intrepid::RealSpaceTools<Real>::scale(valN_scal,*(fe_vol_->N()),static_cast<Real>(-1));
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  valN_scal,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
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
    // COMPUTE NONLINEAR TERM
    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real>> valL_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(valL_eval, l_coeff);
    Intrepid::FieldContainer<Real> d2phi_valU_eval(c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        d2phi_valU_eval(i,j) = (*valL_eval)(i,j)*evaluateReaction((*valU_eval)(i,j),2);
      }
    }
    Intrepid::FieldContainer<Real> NLd2phiU(c,f,p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(NLd2phiU,
                                                                d2phi_valU_eval,
                                                                *(fe_vol_->N()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  NLd2phiU,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Semilinear:Hessian_12: Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Semilinear:Hessian_21: Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Semilinear:Hessian_22: Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    // GET DIMENSIONS
    const int c = fe_vol_->N()->dimension(0);
    const int f = fe_vol_->N()->dimension(1);
    // INITIALIZE RIESZ
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    //*riesz = *fe_vol_->stiffMat();
    //Intrepid::RealSpaceTools<Real>::add(*riesz,*(fe_vol_->massMat()));
    *riesz = *fe_vol_->massMat();
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    // GET DIMENSIONS
    const int c = fe_vol_->N()->dimension(0);
    const int f = fe_vol_->N()->dimension(1);
    // INITIALIZE RIESZ
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    //*riesz = *fe_vol_->stiffMat();
    //Intrepid::RealSpaceTools<Real>::add(*riesz,*(fe_vol_->massMat()));
    *riesz = *fe_vol_->massMat();
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields() {
    return basisPtrs_;
  }

  void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real>> &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds) {
    // Finite element definition.
    fe_vol_ = ROL::makePtr<FE<Real>>(volCellNodes,basisPtr_,cellCub_);
  }

  const ROL::Ptr<FE<Real>> getFE(void) const {
    return fe_vol_;
  }

private:

  Real evaluateDiffusivity(const std::vector<Real> &x) const {
    return static_cast<Real>(1.e-1);
  }

  void evaluateAdvection(std::vector<Real> &v, const std::vector<Real> &x) const {
    const Real a(2.5), b(7.5);
    v[0] = b - a*x[0];
    v[1] =     a*x[1];
  }

  Real evaluateLoad(const std::vector<Real> &x) const {
    const int dim = x.size();
    const Real zero(0), center(0.1), one(1), rad(0.07);
    Real norm(0);
    for (int i = 0; i < dim; ++i) {
      norm += std::pow(x[i]-center,2);
    }
    norm = std::sqrt(norm);
    return (norm < rad ? one : zero);
  }

  void computeCoefficients(Intrepid::FieldContainer<Real> &diff,
                           Intrepid::FieldContainer<Real> &adv,
                           Intrepid::FieldContainer<Real> &rhs) const {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN()->dimension(0);
    const int p = fe_vol_->gradN()->dimension(2);
    const int d = fe_vol_->gradN()->dimension(3);
    std::vector<Real> pt(d), vel(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*fe_vol_->cubPts())(i,j,k);
        }
        // Compute scaled Debye length
        diff(i,j) = evaluateDiffusivity(pt);
        // Compute advection velocity field V
        evaluateAdvection(vel,pt);
        for (int k = 0; k < d; ++k) {
          adv(i,j,k) = adv[k];
        }
        // Compute forcing term f
        rhs(i,j) = -evaluateLoad(pt);
      }
    }
  }

  Real evaluateReaction(const Real u, const int deriv = 0) const {
    const Real zero(0), one(1), three(3), six(6);
    Real val(0);
    if (type_ == 0) {
      val = (deriv == 0 ? u : (deriv == 1 ? one : zero));
    }
    else if (type_ == 1) {
      val = (deriv == 0 ? std::pow(u, 3) - u
          : (deriv == 1 ? three * std::pow(u, 2) - one
          : six * u));
    }
    else if (type_ == 2) {
      val = (deriv == 0 ? std::exp(u) - std::exp(-u)
          : (deriv == 1 ? std::exp(u) + std::exp(-u)
          : std::exp(u) - std::exp(-u)));
    }
    else {
      val = (deriv == 0 ? u : (deriv == 1 ? one : zero));
    }
    return val;
  }
}; // PDE_Semilinear

#endif
