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
#include "Intrepid_CellTools.hpp"

#include "Teuchos_RCP.hpp"

template <class Real>
class PDE_Poisson : public PDE<Real> {
private:
  Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtr_;
  std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;
  Teuchos::RCP<Intrepid::Cubature<Real> > cellCub_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > volCellNodes_;
  Teuchos::RCP<FE<Real> > fe_vol_;

public:
  PDE_Poisson(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("PDE Poisson").get("Basis Order",1);
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
  }

  void residual(Teuchos::RCP<Intrepid::FieldContainer<Real> > & res,
                Teuchos::RCP<Intrepid::FieldContainer<Real> > & u_coeff,
                Teuchos::RCP<Intrepid::FieldContainer<Real> > & z_coeff) {
    int c = u_coeff->dimension(0);
    int p = cellCub_->getNumPoints();
    int f = basisPtr_->getCardinality();
    int d = cellCub_->getDimension();
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    res = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res, *gradU_eval, *(fe_vol_->gradNdetJ()), Intrepid::COMP_CPP, false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res, *valU_eval, *(fe_vol_->NdetJ()), Intrepid::COMP_CPP, true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res, *valZ_eval, *(fe_vol_->NdetJ()), Intrepid::COMP_CPP, true);
  }

  void Jacobian_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & jac,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & u_coeff,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & z_coeff) {
    int c = u_coeff->dimension(0);
    int f = basisPtr_->getCardinality();
    jac = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac, *(fe_vol_->gradN()), *(fe_vol_->gradNdetJ()), Intrepid::COMP_CPP, false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac, *(fe_vol_->N()), *(fe_vol_->NdetJ()), Intrepid::COMP_CPP, true);
  }

  void Jacobian_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & jac,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & u_coeff,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & z_coeff) {
    int c = u_coeff->dimension(0);
    int f = basisPtr_->getCardinality();
    jac = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac, *(fe_vol_->N()), *(fe_vol_->NdetJ()), Intrepid::COMP_CPP, false);
  }

  void Hessian_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & u_coeff,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & z_coeff,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & l_coeff) {
    throw Exception::Zero(">>> Zero Hessian.");
  }

  void Hessian_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & u_coeff,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & z_coeff,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & l_coeff) {
    throw Exception::Zero(">>> Zero Hessian.");
  }

  void Hessian_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & u_coeff,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & z_coeff,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & l_coeff) {
    throw Exception::Zero(">>> Zero Hessian.");
  }

  void Hessian_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & u_coeff,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & z_coeff,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & l_coeff) {
    throw Exception::Zero(">>> Zero Hessian.");
  }

  std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > getFields() {
    return basisPtrs_;
  }

  void setCellNodes(const Teuchos::RCP<Intrepid::FieldContainer<Real> > &volCellNodes,
                    const std::vector<Teuchos::RCP<BoundaryCells<Real> > > &bdryCells) {
    volCellNodes_ = volCellNodes;
    // Finite element definition.
    fe_vol_ = Teuchos::rcp(new FE<Real>(volCellNodes_,basisPtr_,cellCub_));
  }

}; // PDE_Poisson

#endif
