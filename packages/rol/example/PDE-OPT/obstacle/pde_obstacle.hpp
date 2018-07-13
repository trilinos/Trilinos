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
    \brief Implements the local PDE interface for the obstacle problem.
*/

#ifndef PDE_POISSON_HPP
#define PDE_POISSON_HPP

#include "../TOOLS/pde.hpp"
#include "../TOOLS/fe.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real>
class PDE_Obstacle : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real> > cellCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real> > volCellNodes_;
  // Finite element definition
  ROL::Ptr<FE<Real> > fe_vol_;
  // Load term
  ROL::Ptr<Intrepid::FieldContainer<Real> > load_;
  Real coeff_;
  // Use Riesz map for state variables?
  bool useStateRiesz_;

  Real evaluateLoad(const std::vector<Real> &coord) const {
    return static_cast<Real>(2)*coeff_;
  }

public:
  PDE_Obstacle(Teuchos::ParameterList &parlist) {
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
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();                  // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                                 // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 2);              // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                                 // create default cubature
    coeff_ = parlist.sublist("Problem").get("Load Magnitude", 1.0);                    // Pointwise load magnitude
    useStateRiesz_ = parlist.sublist("Problem").get("Use State Riesz Map", true);      // use Riesz map for state variables?
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
    // COMPUTE STIFFNESS TERM
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradU_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *gradU_eval,
                                                  *(fe_vol_->gradNdetJ()),
                                                  Intrepid::COMP_CPP,
                                                  false);
//    // ADD MASS TERM
//    ROL::Ptr<Intrepid::FieldContainer<Real> > valU_eval
//      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
//    fe_vol_->evaluateValue(valU_eval, u_coeff);
//    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
//                                                  *valU_eval,
//                                                  *(fe_vol_->NdetJ()),
//                                                  Intrepid::COMP_CPP,
//                                                  true);
    // ADD LOAD TERM
    Intrepid::FunctionSpaceTools::integrate<Real>(*res, 
                                                  *load_,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP,
                                                  true);
  }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    // INITIALIZE JACOBIAN
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // COMPUTE STIFFNESS TERM
    *jac = *fe_vol_->stiffMat();
//    // ADD MASS TERM
//    Intrepid::RealSpaceTools<Real>::add(*jac, *fe_vol_->massMat());
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Obstacle::Hessian_11): Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    // Optionally disable Riesz map ...
    if (!useStateRiesz_) {
      throw Exception::NotImplemented(">>> (PDE_Obstacle::RieszMap_1): Not implemented.");
    }

    // ...otherwise ...

    int c = fe_vol_->N()->dimension(0);
    int f = fe_vol_->N()->dimension(1);
    // INITIALIZE RIESZ MAP
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    *riesz = *fe_vol_->stiffMat();
    Intrepid::RealSpaceTools<Real>::add(*riesz,*(fe_vol_->massMat()));
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > getFields() {
    return basisPtrs_;
  }

  void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real> > &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds) {
    volCellNodes_ = volCellNodes;
    // Finite element definition.
    fe_vol_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtr_,cellCub_);
    // Compute load function
    computeLoad();
  }

  void computeLoad(void) {
    int c = fe_vol_->gradN()->dimension(0);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    std::vector<Real> coord(d);
    load_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          coord[k] = (*fe_vol_->cubPts())(i,j,k);
        }
        (*load_)(i,j) = -evaluateLoad(coord);
      }
    }
  }

  const ROL::Ptr<FE<Real> > getFE(void) const {
    return fe_vol_;
  }

  const ROL::Ptr<Intrepid::FieldContainer<Real> > getCellNodes(void) const {
    return volCellNodes_;
  }

}; // PDE_Obstacle

#endif
