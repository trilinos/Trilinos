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

/*! \file  pde_helmholtz.hpp
    \brief Implements the local PDE interface for the optimal control of
           Helmholtz.
*/

#ifndef PDE_HELMHOLTZ_REAL_HPP
#define PDE_HELMHOLTZ_REAL_HPP

#include "../../TOOLS/pde.hpp"
#include "../../TOOLS/fe.hpp"
#include "../../TOOLS/fieldhelper.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"


template <class Real>
class PDE_Helmholtz_Real : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>> basisPtr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real>> cellCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real> > volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;
  std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real>> fe_;

  Real innerAnnulusRadius_;
  Real outerAnnulusRadius_;
  int example_;
  Real waveNumber_;

  ROL::Ptr<Intrepid::FieldContainer<Real>> ctrlWeight_;
  ROL::Ptr<Intrepid::FieldContainer<Real>> ctrlJac_;

  bool insideControlDomain(const std::vector<Real> &x) const {
    bool val = true;
    if (example_==1) {
      const Real eps = std::sqrt(ROL::ROL_EPSILON<Real>());
      Real xnorm(0);
      const int d = x.size();
      for (int i = 0; i < d; ++i) {
        xnorm += x[i]*x[i];
      }
      xnorm = std::sqrt(xnorm);
      val = (xnorm <= outerAnnulusRadius_+eps && xnorm >= innerAnnulusRadius_-eps);
    }
    return val;
  }

  void computeControlWeight(void) {
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
   
    ctrlWeight_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);

    const Real zero(0), one(1);
    bool inside(false);
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      inside = false;
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (*fe_->cubPts())(i,j,k);
        }
        if ( insideControlDomain(x) ) {
          inside = true;
          break;
        }
      }
      for (int j = 0; j < p; ++j) {
        (*ctrlWeight_)(i,j) = (inside ? one : zero);
      }
    }
  }

  void buildControlJacobian(void) {
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);

    // Build force/control term
    ROL::Ptr<Intrepid::FieldContainer<Real>> F
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*F, *ctrlWeight_, *fe_->N());
    ctrlJac_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*ctrlJac_,
                                                  *F,
                                                  *fe_->NdetJ(),
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::RealSpaceTools<Real>::scale(*ctrlJac_,static_cast<Real>(-1));
  }

public:
  PDE_Helmholtz_Real(Teuchos::ParameterList &parlist)
    : innerAnnulusRadius_(2.5), outerAnnulusRadius_(2.6),
      example_(parlist.sublist("Problem").get("Example",1)),
      waveNumber_(parlist.sublist("Problem").get("Wave Number",10.0)) {
    // Finite element fields.
    basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();            // get the cell type from the basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                           // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 2);        // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                           // create default cubature
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real>> & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real one(1), two(2);
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    // Initialize residuals.
    res = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f);
    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> valz_eval, valu_eval, gradu_eval;
    valz_eval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    valu_eval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    gradu_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_->evaluateValue(valz_eval, z_coeff);
    fe_->evaluateValue(valu_eval, u_coeff);
    fe_->evaluateGradient(gradu_eval, u_coeff);
    // Integrate PDE term
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *valu_eval,    // U
                                                  *fe_->NdetJ(), // N
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::RealSpaceTools<Real>::scale(*res, -std::pow(waveNumber_,two));
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *gradu_eval,       // grad U
                                                  *fe_->gradNdetJ(), // grad N
                                                  Intrepid::COMP_CPP,
                                                  true);
    // Build control term
    ROL::Ptr<Intrepid::FieldContainer<Real>> F;
    F = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*F, *ctrlWeight_, *valz_eval);
    Intrepid::RealSpaceTools<Real>::scale(*F,-one);
    // Volumetric intregration
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *F,              // F
                                                  *fe_->NdetJ(),   // N
                                                  Intrepid::COMP_CPP,
                                                  true);

  }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real two(2);
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    // Initialize Jacobian.
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f);
    // Add PDE terms
    *jac = *fe_->massMat();
    Intrepid::RealSpaceTools<Real>::scale(*jac, -std::pow(waveNumber_,two));
    Intrepid::RealSpaceTools<Real>::add(*jac, *fe_->stiffMat());
  }


  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    // Initialize Jacobian.
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f);
    // Add control term
    *jac = *ctrlJac_;
    Intrepid::RealSpaceTools<Real>::scale(*jac, static_cast<Real>(-1));
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Helmholtz::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Helmholtz::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Helmholtz::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Helmholtz::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    // Initialize Jacobian.
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f);
    *riesz = *fe_->stiffMat();
    Intrepid::RealSpaceTools<Real>::add(*riesz, *fe_->massMat());
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    // Initialize Jacobian.
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f);
    *riesz = *fe_->massMat();
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
    fe_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtr_,cellCub_);
    // Compute control weight
    computeControlWeight();
    buildControlJacobian();
  }

  const ROL::Ptr<FE<Real>> getFE(void) const {
    return fe_;
  }

}; // PDE_Helmholtz


#endif
