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

#ifndef PDE_POISSON_HPP
#define PDE_POISSON_HPP

#include "../TOOLS/pde.hpp"
#include "../TOOLS/fe.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "Teuchos_RCP.hpp"

template <class Real>
class PDE_Poisson_TopOpt : public PDE<Real> {
private:
  // Finite element basis information
  Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtr_;
  std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;
  // Cell cubature information
  Teuchos::RCP<Intrepid::Cubature<Real> > cellCub_;
  // Cell node information
  Teuchos::RCP<Intrepid::FieldContainer<Real> > volCellNodes_;
  std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;
  std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
  // Finite element definition
  Teuchos::RCP<FE<Real> > fe_vol_;
  // Force function evaluated at cubature points
  Teuchos::RCP<Intrepid::FieldContainer<Real> > force_eval_;
  // Inputs
  Real z0_;
  Real p_;

  Real ForceFunc(const std::vector<Real> &x) const {
    int dim = x.size();
    std::vector<Real> l1(dim,0), u1(dim,0);
    std::vector<Real> l2(dim,0), u2(dim,0);
    if ( dim > 0 ) { l1[0] = 0.4; u1[0] = 0.6; l2[0] = 0.4; u2[0] = 0.6; }
    if ( dim > 1 ) { l1[1] = 0.0; u1[1] = 0.1; l2[1] = 0.9; u2[1] = 1.0; }
    if ( dim > 2 ) { l1[2] = 0.4; u1[2] = 0.6; l2[2] = 0.4; u2[2] = 0.6; }
    bool flag1 = true, flag2 = true;
    for (int i = 0; i < dim; ++i) {
      if (x[i] < l1[i] || x[i] > u1[i]) {
        flag1 *= false;
      }
      if (x[i] < l2[i] || x[i] > u2[i]) {
        flag2 *= false;
      }
    }
    return (flag1 ? 1 : 0) + (flag2 ? -1 : 0);
  }

public:
  PDE_Poisson_TopOpt(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("PDE Poisson TopOpt").get("Basis Order",1);
    if (basisOrder == 1) {
      basisPtr_ = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >);
    }
    else if (basisOrder == 2) {
      basisPtr_ = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >);
    }
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();        // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                       // create cubature factory
    int cubDegree = parlist.sublist("PDE Poisson").get("Cubature Degree",2); // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                       // create default cubature

    z0_ = parlist.sublist("PDE Poisson TopOpt").get("Minimum Conductivity",1.e-4);
    p_  = parlist.sublist("PDE Poisson TopOpt").get("SIMP Parameter",3);
  }

  void residual(Teuchos::RCP<Intrepid::FieldContainer<Real> > & res,
                const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff) {
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    res = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
    // Build density-dependent conductivity function
    const Real one(1);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valZ_eval)(i,j) = z0_ + (one - z0_)*std::pow((*valZ_eval)(i,j),p_);
      }
    }
    // Build flux function
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
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
    // Integrate reaction term
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    Intrepid::FieldContainer<Real> KU(c,p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(KU,
                                                               *valZ_eval,
                                                               *valU_eval);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  KU,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // Add force term
    Intrepid::FieldContainer<Real> mforce_eval(c,p);
    Intrepid::RealSpaceTools<Real>::scale(mforce_eval,*force_eval_,static_cast<Real>(-1));
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  mforce_eval,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
  }

  void Jacobian_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & jac,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff) {
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    jac = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    // Build density-dependent conductivity function
    const Real one(1);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valZ_eval)(i,j) = z0_ + (one - z0_)*std::pow((*valZ_eval)(i,j),p_);
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
    // Integrate reaction term
    Intrepid::FieldContainer<Real> KN(c,f,p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(KN,
                                                               *valZ_eval,
                                                               *(fe_vol_->N()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  KN,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
  }

  void Jacobian_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & jac,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff) {
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    jac = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    // Build density-dependent conductivity function
    const Real one(1);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valZ_eval)(i,j) = p_*(one - z0_)*std::pow((*valZ_eval)(i,j),p_-one);
      }
    }
    // Build derivative of conductivity function
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    Intrepid::FieldContainer<Real> dKN(c,f,p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(dKN,
                                                                *valZ_eval,
                                                                *(fe_vol_->N()));
    // Integrate stiffness term
    Intrepid::FieldContainer<Real> gradUgradN(c,f,p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(gradUgradN,
                                                             *gradU_eval,
                                                             (*fe_vol_->gradNdetJ()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  gradUgradN,
                                                  dKN,
                                                  Intrepid::COMP_CPP, false);
    // Integrate reaction term
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    Intrepid::FieldContainer<Real> UN(c,f,p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(UN,
                                                                *valU_eval,
                                                                (*fe_vol_->NdetJ()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  UN,
                                                  dKN,
                                                  Intrepid::COMP_CPP, true);
  }

  void Hessian_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff) {
    throw Exception::Zero(">>> (PDE_Poisson_TopOpt::Hessian_11): Zero Hessian.");
  }

  void Hessian_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff) {
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    hess = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    // Build density-dependent conductivity function
    const Real one(1);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valZ_eval)(i,j) = p_*(one - z0_)*std::pow((*valZ_eval)(i,j),p_-one);
      }
    }
    // Build derivative of conductivity function
    Intrepid::FieldContainer<Real> dKN(c,f,p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(dKN,
                                                                *valZ_eval,
                                                                *(fe_vol_->N()));
    // Integrate stiffness term
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradL_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    fe_vol_->evaluateGradient(gradL_eval, l_coeff);
    Intrepid::FieldContainer<Real> gradLgradN(c,f,p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(gradLgradN,
                                                             *gradL_eval,
                                                             (*fe_vol_->gradNdetJ()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  dKN,
                                                  gradLgradN,
                                                  Intrepid::COMP_CPP, false);
    // Integrate reaction term
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valL_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_vol_->evaluateValue(valL_eval, l_coeff);
    Intrepid::FieldContainer<Real> LN(c,f,p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(LN,
                                                                *valL_eval,
                                                                (*fe_vol_->NdetJ()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  dKN,
                                                  LN,
                                                  Intrepid::COMP_CPP, true);
  }

  void Hessian_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff) {
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    hess = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    // Build density-dependent conductivity function
    const Real one(1);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valZ_eval)(i,j) = p_*(one - z0_)*std::pow((*valZ_eval)(i,j),p_-one);
      }
    }
    // Build derivative of conductivity function
    Intrepid::FieldContainer<Real> dKN(c,f,p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(dKN,
                                                                *valZ_eval,
                                                                *(fe_vol_->N()));
    // Integrate stiffness term
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradL_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    fe_vol_->evaluateGradient(gradL_eval, l_coeff);
    Intrepid::FieldContainer<Real> gradLgradN(c,f,p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(gradLgradN,
                                                             *gradL_eval,
                                                             (*fe_vol_->gradNdetJ()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  gradLgradN,
                                                  dKN,
                                                  Intrepid::COMP_CPP, false);
    // Integrate reaction term
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valL_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_vol_->evaluateValue(valL_eval, l_coeff);
    Intrepid::FieldContainer<Real> LN(c,f,p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(LN,
                                                                *valL_eval,
                                                                (*fe_vol_->NdetJ()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  LN,
                                                  dKN,
                                                  Intrepid::COMP_CPP, true);
  }

  void Hessian_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff) {
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    hess = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    // Build density-dependent conductivity function
    const Real one(1), two(2);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valZ_eval)(i,j) = p_*(p_-one)*(one - z0_)*std::pow((*valZ_eval)(i,j),p_-two);
      }
    }
    // Build derivative of conductivity function
    Intrepid::FieldContainer<Real> dKN(c,f,p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(dKN,
                                                                *valZ_eval,
                                                                *(fe_vol_->N()));
    // Integrate stiffness term
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradL_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    fe_vol_->evaluateGradient(gradL_eval, l_coeff);
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
    // Integrate reaction term
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valL_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_vol_->evaluateValue(valL_eval, l_coeff);
    Intrepid::FieldContainer<Real> UL(c,p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(UL,
                                                               *valU_eval,
                                                               *valL_eval);
    Intrepid::FieldContainer<Real> ULN(c,f,p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(ULN,
                                                                UL,
                                                                *(fe_vol_->NdetJ()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  dKN,
                                                  ULN,
                                                  Intrepid::COMP_CPP, true);
  }

  std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > getFields() {
    return basisPtrs_;
  }

  void setCellNodes(const Teuchos::RCP<Intrepid::FieldContainer<Real> > &volCellNodes,
                    const std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds) {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fe_vol_ = Teuchos::rcp(new FE<Real>(volCellNodes_,basisPtr_,cellCub_));
    computeForce();
  }

  void computeForce(void) {
    int c = fe_vol_->cubPts()->dimension(0);
    int p = fe_vol_->cubPts()->dimension(1);
    int d = fe_vol_->cubPts()->dimension(2);
    std::vector<Real> pt(d,0);
    force_eval_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,p));
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (*fe_vol_->cubPts())(i,j,k);
        }
        (*force_eval_)(i,j) = ForceFunc(pt);
      }
    }
  }

  const Teuchos::RCP<FE<Real> > getFE(void) const {
    return fe_vol_;
  }

  const Teuchos::RCP<Intrepid::FieldContainer<Real> > getForce(void) const {
    return force_eval_;
  }

}; // PDE_Poisson

#endif
