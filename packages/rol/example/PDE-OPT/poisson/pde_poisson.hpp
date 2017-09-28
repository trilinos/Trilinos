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

#include "Intrepid_HGRAD_LINE_Cn_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "Teuchos_RCP.hpp"

template <class Real>
class PDE_Poisson : public PDE<Real> {
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
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int> > fidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > bdryCellDofValues_;

  bool useStateRiesz_;
  bool useControlRiesz_;
  Real alpha_;

  Real DirichletFunc(const std::vector<Real> & coords, int sideset, int locSideId) const {
    return 0;
  }

  Real evaluateRHS(const std::vector<Real> &x) const {
    const Real pi(M_PI), eight(8);
    Real s1(1), s2(1);
    int dim = x.size();
    for (int i=0; i<dim; ++i) {
      s1 *= std::sin(eight*pi*x[i]);
      s2 *= std::sin(pi*x[i]);
    }
    Real coeff1(64), coeff2(dim);
    return s1/(alpha_*coeff1*coeff2*pi*pi) + coeff2*pi*pi*s2;
  }

  void computeRHS(Teuchos::RCP<Intrepid::FieldContainer<Real> > &rhs) const {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*fe_vol_->cubPts())(i,j,k);
        }
        // Compute forcing term f
        (*rhs)(i,j) = -evaluateRHS(pt);
      }
    }
  }

public:
  PDE_Poisson(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("Problem").get("Basis Order",1);
    int cubDegree  = parlist.sublist("Problem").get("Cubature Degree",4);
    int probDim    = parlist.sublist("Problem").get("Problem Dimension",2);
    if (probDim > 3 || probDim < 1) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> PDE-OPT/poisson/pde_poisson.hpp: Problem dimension is not 1, 2 or 3!");
    }
    if (basisOrder > 2 || basisOrder < 1) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> PDE-OPT/poisson/pde_poisson.hpp: Basis order is not 1 or 2!");
    }
    if (probDim == 1) {
      basisPtr_ = Teuchos::rcp(new Intrepid::Basis_HGRAD_LINE_Cn_FEM<Real, Intrepid::FieldContainer<Real> >(basisOrder,Intrepid::POINTTYPE_EQUISPACED));
    } else if (probDim == 2) {
      if (basisOrder == 1) {
        basisPtr_ = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >);
      }
      else if (basisOrder == 2) {
        basisPtr_ = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >);
      }
    } else if (probDim == 3) {
      if (basisOrder == 1) {
        basisPtr_ = Teuchos::rcp(new Intrepid::Basis_HGRAD_HEX_C1_FEM<Real, Intrepid::FieldContainer<Real> >);
      }
      else if (basisOrder == 2) {
        basisPtr_ = Teuchos::rcp(new Intrepid::Basis_HGRAD_HEX_C2_FEM<Real, Intrepid::FieldContainer<Real> >);
      }
    }
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();
    Intrepid::DefaultCubatureFactory<Real> cubFactory;
    cellCub_ = cubFactory.create(cellType, cubDegree);
    // Problem data.
    useStateRiesz_   = parlist.sublist("Problem").get("Use State Riesz Map", true);      // use Riesz map for state variables?
    useControlRiesz_ = parlist.sublist("Problem").get("Use Control Riesz Map", true);    // use Riesz map for control variables?
    alpha_           = parlist.sublist("Problem").get("Control penalty parameter",1e-2);
  }

  void residual(Teuchos::RCP<Intrepid::FieldContainer<Real> > & res,
                const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // GET DIMENSIONS
    int c = u_coeff->dimension(0);
    int p = cellCub_->getNumPoints();
    int f = basisPtr_->getCardinality();
    int d = cellCub_->getDimension();
    // INITIALIZE RESIDUAL
    res = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
    // COMPUTE STIFFNESS TERM
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradU_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res, *gradU_eval, *(fe_vol_->gradNdetJ()), Intrepid::COMP_CPP, false);
    // COMPUTE RHS
    Teuchos::RCP<Intrepid::FieldContainer<Real> > rhs
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    computeRHS(rhs);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res, *rhs, *(fe_vol_->NdetJ()), Intrepid::COMP_CPP, true);
    // ADD CONTROL TERM TO RESIDUAL
    if ( z_coeff != Teuchos::null ) {
      Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval =
        Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      fe_vol_->evaluateValue(valZ_eval, z_coeff);
      Intrepid::RealSpaceTools<Real>::scale(*valZ_eval,static_cast<Real>(-1));
      Intrepid::FunctionSpaceTools::integrate<Real>(*res, *valZ_eval, *(fe_vol_->NdetJ()), Intrepid::COMP_CPP, true);
    }
    // APPLY DIRICHLET CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        int numLocalSideIds = bdryCellLocIds_[i].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          int numCellsSide = bdryCellLocIds_[i][j].size();
          int numBdryDofs = fidx_[j].size();
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
            for (int l = 0; l < numBdryDofs; ++l) {
              //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx_[j][l];
              (*res)(cidx,fidx_[j][l]) = (*u_coeff)(cidx,fidx_[j][l]) - (*bdryCellDofValues_[i][j])(k,fidx_[j][l]);
            }
          }
        }
      }
    }
  }

  void Jacobian_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & jac,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // GET DIMENSIONS
    int c = u_coeff->dimension(0);
    int f = basisPtr_->getCardinality();
    // INITILAIZE JACOBIAN
    jac = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    // COMPUTE STIFFNESS TERM
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac, *(fe_vol_->gradN()), *(fe_vol_->gradNdetJ()), Intrepid::COMP_CPP, false);
    // APPLY DIRICHLET CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        int numLocalSideIds = bdryCellLocIds_[i].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          int numCellsSide = bdryCellLocIds_[i][j].size();
          int numBdryDofs = fidx_[j].size();
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
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
    }
  }

  void Jacobian_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & jac,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    if ( z_coeff != Teuchos::null ) {
      // GET DIMENSIONS
      int c = u_coeff->dimension(0);
      int f = basisPtr_->getCardinality();
      // INITIALIZE JACOBIAN
      jac = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
      // ADD CONTROL TERM
      Intrepid::FunctionSpaceTools::integrate<Real>(*jac, *(fe_vol_->N()), *(fe_vol_->NdetJ()), Intrepid::COMP_CPP, false);
      Intrepid::RealSpaceTools<Real>::scale(*jac,static_cast<Real>(-1));
      // APPLY DIRICHLET CONDITIONS
      int numSideSets = bdryCellLocIds_.size();
      if (numSideSets > 0) {
        for (int i = 0; i < numSideSets; ++i) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
                for (int m = 0; m < f; ++m) {
                  (*jac)(cidx,fidx_[j][l],m) = static_cast<Real>(0);
                }
              }
            }
          }
        }
      }
    }
    else {
      throw Exception::Zero(">>> (PDE_Poisson::Jacobian_2): Jacobian is zero.");
    }
  }

  void Hessian_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Poisson::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Poisson::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Poisson::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Poisson::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & riesz) {
    // Optionally disable Riesz map ...
    if (!useStateRiesz_) {
      throw Exception::NotImplemented(">>> (PDE_Poisson::RieszMap_1): Not implemented.");
    }

    // ...otherwise ...

    int c = fe_vol_->N()->dimension(0);
    int f = fe_vol_->N()->dimension(1);
    // INITILAIZE JACOBIAN
    riesz = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    *riesz = *fe_vol_->stiffMat();
    Intrepid::RealSpaceTools<Real>::add(*riesz,*(fe_vol_->massMat()));
  }

  void RieszMap_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & riesz) {
    // Optionally disable Riesz map ...
    if (!useControlRiesz_) {
      throw Exception::NotImplemented(">>> (PDE_Poisson::RieszMap_2): Not implemented.");
    }

    // ...otherwise ...

    int c = fe_vol_->N()->dimension(0);
    int f = fe_vol_->N()->dimension(1);
    // INITILAIZE JACOBIAN
    riesz = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    *riesz = *fe_vol_->massMat();
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
        bdryCellDofValues_[i][j] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
        Teuchos::RCP<Intrepid::FieldContainer<Real> > coords =
          Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, d));
        if (c > 0) {
          fe_vol_->computeDofCoords(coords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<f; ++l) {
            std::vector<Real> dofpoint(d);
            for (int m=0; m<d; ++m) {
              dofpoint[m] = (*coords)(k, l, m);
            }
            (*bdryCellDofValues_[i][j])(k, l) = DirichletFunc(dofpoint, i, j);
          }
        }
      }
    }
  }

  const Teuchos::RCP<FE<Real> > getFE(void) const {
    return fe_vol_;
  }

}; // PDE_Poisson

#endif
