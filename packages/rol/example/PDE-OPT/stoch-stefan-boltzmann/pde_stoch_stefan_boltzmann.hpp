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
    \brief Implements the local PDE interface for the Stefan_Boltzmann control problem.
*/

#ifndef PDEOPT_PDE_STOCH_STEFAN_BOLTZMANN_HPP
#define PDEOPT_PDE_STOCH_STEFAN_BOLTZMANN_HPP

#include "../TOOLS/pde.hpp"
#include "../TOOLS/fe.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "Teuchos_RCP.hpp"

template <class Real>
class StochasticStefanBoltzmannPDE : public PDE<Real> {
private:

  // Finite element basis information
  Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtr_;
  std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;
  // Cell cubature information
  Teuchos::RCP<Intrepid::Cubature<Real> > cellCub_;
  Teuchos::RCP<Intrepid::Cubature<Real> > bdryCub_;
  // Cell node information
  Teuchos::RCP<Intrepid::FieldContainer<Real> > volCellNodes_;
  std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;
  std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
  // Finite element definition
  Teuchos::RCP<FE<Real> > fe_vol_;
  std::vector<std::vector<Teuchos::RCP<FE<Real> > > > fe_bdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int> > fidx_;

  const Real scale_;
  Real xmid_;
  Real engTemp_;
  Real airTemp_;


public:

  StochasticStefanBoltzmannPDE(Teuchos::ParameterList &parlist) : scale_(1) {
    xmid_ = parlist.sublist("Geometry").get<Real>("Step height");
    engTemp_ = parlist.sublist("Problem").get("Engine: Ambient Temperature",450.0);
    airTemp_ = parlist.sublist("Problem").get("Air: Ambient Temperature",293.0);
    // Finite element fields.
    int basisOrder = parlist.sublist("Problem").get("Basis Order",1);
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
    int cubDegree = parlist.sublist("Problem").get("Volume Cubature Degree",2); // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                       // create default cubature

    int d = cellType.getDimension();
    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2); // set cubature degree, e.g., 2
    bdryCub_ = cubFactory.create(bdryCellType, bdryCubDegree);
  }
  
  void residual(Teuchos::RCP<Intrepid::FieldContainer<Real> > & res,
                const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    // INITIALIZE RESIDUAL
    res = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
    // EVALUATE STATE ON FE BASIS
    Teuchos::RCP<Intrepid::FieldContainer<Real> > U_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_vol_->evaluateValue(U_eval, u_coeff);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    // COMPUTE CONSTANT PDE COEFFICIENTS
    Teuchos::RCP<Intrepid::FieldContainer<Real> > V
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > rhs
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    computeCoefficients(V,rhs);
    // COMPUTE DIFFUSIVITY
    Teuchos::RCP<Intrepid::FieldContainer<Real> > kappa
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    computeDiffusivity(kappa,U_eval,0);
    // MULTIPLY kappa * grad(U)
    Intrepid::FieldContainer<Real> kappa_gradU(c, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(kappa_gradU,
                                                               *kappa,
                                                               *gradU_eval);
    // INTEGRATE (kappa * grad(U)) . grad(N)
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  kappa_gradU,
                                                  *(fe_vol_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // MULTIPLY V . grad(U)
    Intrepid::FieldContainer<Real> V_gradU(c, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(V_gradU,
                                                            *V,
                                                            *gradU_eval);
    // INTEGRATE (V . grad(U)) * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  V_gradU,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // ADD RHS TO RESIDUAL
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *rhs,
                                                  (*fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);



    // APPLY NEUMANN CONDITIONS: Sideset 1, 3, 6
    // ---> Nothing to do
    int numLocalSideIds(0);
    // APPLY STEFAN-BOLTZMANN CONDITIONS: Sideset 2, 4, 5
    std::vector<int> sidesets = {2, 4, 5};
    for (int i = 0; i < 3; ++i) {
      numLocalSideIds = bdryCellLocIds_[sidesets[i]].size();
      const int numCubPerSide = bdryCub_->getNumPoints();
      for (int j = 0; j < numLocalSideIds; ++j) {
        int numCellsSide = bdryCellLocIds_[sidesets[i]][j].size();
        if (numCellsSide) {
          // Get U coefficients on Stefan-Boltzmann boundary
          Teuchos::RCP<Intrepid::FieldContainer<Real > > u_coeff_bdry
            = getBoundaryCoeff(*u_coeff, sidesets[i], j);
          // Evaluate U on FE basis
          Teuchos::RCP<Intrepid::FieldContainer<Real > > valU_eval_bdry
            = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
          fe_bdry_[sidesets[i]][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
          // Compute Stefan-Boltzmann residual
          Teuchos::RCP< Intrepid::FieldContainer<Real> > sb_valU
            = Teuchos::rcp( new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
          computeStefanBoltzmann(sb_valU,valU_eval_bdry,sidesets[i],j,0);
          Intrepid::FieldContainer<Real> sbRes(numCellsSide, f);
          Intrepid::FunctionSpaceTools::integrate<Real>(sbRes,
                                                        *sb_valU,
                                                        *(fe_bdry_[sidesets[i]][j]->NdetJ()),
                                                        Intrepid::COMP_CPP, false);
          // Add Stefan-Boltzmann residual to volume residual
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[sidesets[i]][j][k];
            for (int l = 0; l < f; ++l) { 
              (*res)(cidx,l) += sbRes(k,l);
            }
          }
        }
      }
    }
    // APPLY DIRICHLET CONTROLS: Sideset 0
    int sideset = 0;
    numLocalSideIds = bdryCellLocIds_[sideset].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[sideset][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
        int cidx = bdryCellLocIds_[sideset][j][k];
        for (int l = 0; l < numBdryDofs; ++l) {
          (*res)(cidx,fidx_[j][l])
            = (*u_coeff)(cidx,fidx_[j][l]) - (*z_coeff)(cidx,fidx_[j][l]);
        }
      }
    }
  }
  
  void Jacobian_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & jac,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    // INITILAIZE JACOBIAN
    jac = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    // EVALUATE STATE ON FE BASIS
    Teuchos::RCP<Intrepid::FieldContainer<Real> > U_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_vol_->evaluateValue(U_eval, u_coeff);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    // COMPUTE CONSTNAT PDE COEFFICIENTS
    Teuchos::RCP<Intrepid::FieldContainer<Real> > V
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > rhs
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    computeCoefficients(V,rhs);
    // COMPUTE DIFFUSIVITY
    Teuchos::RCP<Intrepid::FieldContainer<Real> > kappa
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    computeDiffusivity(kappa,U_eval,0);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > d_kappa
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    computeDiffusivity(d_kappa,U_eval,1);
    // MULTIPLY kappa * grad(N)
    Intrepid::FieldContainer<Real> kappa_gradN(c, f, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(kappa_gradN,
                                                                *kappa,
                                                                *(fe_vol_->gradN()));
    // INTEGRATE (kappa * grad(N)) . grad(N)
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  kappa_gradN,
                                                  *(fe_vol_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // MULTIPLY d_kappa * grad(U)
    Intrepid::FieldContainer<Real> d_kappa_gradU(c, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(d_kappa_gradU,
                                                               *d_kappa,
                                                               *gradU_eval);
    // MULTIPLY (d_kappa * grad(U)) . grad(N)
    Intrepid::FieldContainer<Real> d_kappa_gradU_gradN(c, f, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(d_kappa_gradU_gradN,
                                                             d_kappa_gradU,
                                                             *(fe_vol_->gradNdetJ()));
    // INTEGRATE (d_kappa * grad(U)) . grad(N) * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  d_kappa_gradU_gradN,
                                                  *(fe_vol_->N()),
                                                  Intrepid::COMP_CPP, true);
    // MULTIPLY V . grad(N)
    Intrepid::FieldContainer<Real> V_gradN(c, f, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(V_gradN,
                                                             *V,
                                                             *(fe_vol_->gradN()));
    // INTEGRATE (V . grad(U)) * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *(fe_vol_->NdetJ()),
                                                  V_gradN,
                                                  Intrepid::COMP_CPP, true);

    // APPLY NEUMANN CONDITIONS: Sideset 1, 3, 6
    // ---> Nothing to do
    int numLocalSideIds(0);
    // APPLY STEFAN-BOLTZMANN CONDITIONS: Sideset 2, 4, 5
    std::vector<int> sidesets = {2, 4, 5};
    for (int i = 0; i < 3; ++i) {
      numLocalSideIds = bdryCellLocIds_[sidesets[i]].size();
      const int numCubPerSide = bdryCub_->getNumPoints();
      for (int j = 0; j < numLocalSideIds; ++j) {
        int numCellsSide = bdryCellLocIds_[sidesets[i]][j].size();
        if (numCellsSide) {
          // Get U coefficients on Stefan-Boltzmann boundary
          Teuchos::RCP<Intrepid::FieldContainer<Real > > u_coeff_bdry
            = getBoundaryCoeff(*u_coeff, sidesets[i], j);
          // Evaluate U on FE basis
          Teuchos::RCP<Intrepid::FieldContainer<Real > > valU_eval_bdry
            = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
          fe_bdry_[sidesets[i]][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
          // Compute Stefan-Boltzmann residual
          Teuchos::RCP< Intrepid::FieldContainer<Real> > sb_derivU
            = Teuchos::rcp( new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
          computeStefanBoltzmann(sb_derivU,valU_eval_bdry,sidesets[i],j,1);
          Intrepid::FieldContainer<Real> sb_derivU_N(numCellsSide, f, numCubPerSide);
          Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(sb_derivU_N,
                                                                      *sb_derivU,
                                                                      *(fe_bdry_[sidesets[i]][j]->N()));
          Intrepid::FieldContainer<Real> sbJac(numCellsSide, f, f);
          Intrepid::FunctionSpaceTools::integrate<Real>(sbJac,
                                                        sb_derivU_N,
                                                        *(fe_bdry_[sidesets[i]][j]->NdetJ()),
                                                        Intrepid::COMP_CPP, false);
          // Add Stefan-Boltzmann residual to volume residual
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[sidesets[i]][j][k];
            for (int l = 0; l < f; ++l) { 
              for (int m = 0; m < f; ++m) { 
                (*jac)(cidx,l,m) += sbJac(k,l,m);
              }
            }
          }
        }
      }
    }
    // APPLY DIRICHLET CONDITIONS: Sideset 0
    int sideset = 0;
    numLocalSideIds = bdryCellLocIds_[sideset].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[sideset][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
        int cidx = bdryCellLocIds_[sideset][j][k];
        for (int l = 0; l < numBdryDofs; ++l) {
          for (int m = 0; m < f; ++m) {
            (*jac)(cidx,fidx_[j][l],m) = static_cast<Real>(0);
          }
          (*jac)(cidx,fidx_[j][l],fidx_[j][l]) = static_cast<Real>(1);
        }
      }
    }
  }
  
  void Jacobian_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & jac,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    // INITILAIZE JACOBIAN
    jac = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    // APPLY DIRICHLET CONTROLS: Sideset 0
    int sideset = 0;
    int numLocalSideIds = bdryCellLocIds_[sideset].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[sideset][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
        int cidx = bdryCellLocIds_[sideset][j][k];
        for (int l = 0; l < numBdryDofs; ++l) {
          (*jac)(cidx,fidx_[j][l],fidx_[j][l]) = static_cast<Real>(-1);
        }
      }
    }
  }

  void Hessian_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    // INITILAIZE JACOBIAN
    hess = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    // APPLY DIRICHLET CONDITIONS TO LAGRANGE MULTIPLIERS: Sideset 0, 4
    Teuchos::RCP<Intrepid::FieldContainer<Real> > l_coeff_dbc
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(*l_coeff));
    int sideset = 0;
    int numLocalSideIds = bdryCellLocIds_[sideset].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[sideset][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
        int cidx = bdryCellLocIds_[sideset][j][k];
        for (int l = 0; l < numBdryDofs; ++l) {
          (*l_coeff_dbc)(cidx,fidx_[j][l]) = static_cast<Real>(0);
        }
      }
    }
    // EVALUATE STATE ON FE BASIS
    Teuchos::RCP<Intrepid::FieldContainer<Real> > U_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_vol_->evaluateValue(U_eval, u_coeff);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradL_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    fe_vol_->evaluateGradient(gradL_eval, l_coeff_dbc);
    // COMPUTE DIFFUSIVITY
    Teuchos::RCP<Intrepid::FieldContainer<Real> > d1_kappa
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    computeDiffusivity(d1_kappa,U_eval,1);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > d2_kappa
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    computeDiffusivity(d2_kappa,U_eval,2);
    // MULTIPLY d1_kappa * grad(L)
    Intrepid::FieldContainer<Real> d1_kappa_gradL(c, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(d1_kappa_gradL,
                                                               *d1_kappa,
                                                               *gradL_eval);
    // MULTIPLY (d_1kappa * grad(L)) . grad(N)det(J)
    Intrepid::FieldContainer<Real> d1_kappa_gradL_gradNdetJ(c, f, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(d1_kappa_gradL_gradNdetJ,
                                                             d1_kappa_gradL,
                                                             *(fe_vol_->gradNdetJ()));
    // INTEGRATE (d1_kappa * grad(L)) . grad(N) * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *(fe_vol_->N()),
                                                  d1_kappa_gradL_gradNdetJ,
                                                  Intrepid::COMP_CPP, false);
    // MULTIPLY d1_kappa * grad(L) . grad(N)
    Intrepid::FieldContainer<Real> d1_kappa_gradL_gradN(c, f, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(d1_kappa_gradL_gradN,
                                                             d1_kappa_gradL,
                                                             *(fe_vol_->gradN()));
    // INTEGRATE (d1_kappa * grad(L)) . grad(N) * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  d1_kappa_gradL_gradN,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // MULTIPLY grad(U) . grad(L)
    Intrepid::FieldContainer<Real> gradU_gradL(c, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(gradU_gradL,
                                                            *gradU_eval,
                                                            *gradL_eval);
    // MULTIPLY d2_kappa * grad(U) . grad(L)
    Intrepid::FieldContainer<Real> d2_kappa_gradU_gradL(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(d2_kappa_gradU_gradL,
                                                               *d2_kappa,
                                                               gradU_gradL);
    // MULTIPLY d2_kappa * grad(U) . grad(L) * N
    Intrepid::FieldContainer<Real> d2_kappa_gradU_gradL_N(c, f, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(d2_kappa_gradU_gradL_N,
                                                                d2_kappa_gradU_gradL,
                                                                *(fe_vol_->N()));
    // INTEGRATE (d2_kappa * grad(U) . grad(L) ) * N * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  d2_kappa_gradU_gradL_N,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // APPLY NEUMANN CONDITIONS: Sideset 1, 3, 6
    // ---> Nothing to do
    // APPLY STEFAN-BOLTZMANN CONDITIONS: Sideset 2, 4, 5
    std::vector<int> sidesets = {2, 4, 5};
    for (int i = 0; i < 3; ++i) {
      int numLocalSideIds = bdryCellLocIds_[sidesets[i]].size();
      const int numCubPerSide = bdryCub_->getNumPoints();
      for (int j = 0; j < numLocalSideIds; ++j) {
        int numCellsSide = bdryCellLocIds_[sidesets[i]][j].size();
        if (numCellsSide) {
          // Get U coefficients on Stefan-Boltzmann boundary
          Teuchos::RCP<Intrepid::FieldContainer<Real > > u_coeff_bdry
            = getBoundaryCoeff(*u_coeff, sidesets[i], j);
          Teuchos::RCP<Intrepid::FieldContainer<Real > > l_coeff_bdry
            = getBoundaryCoeff(*l_coeff_dbc, sidesets[i], j);
          // Evaluate U on FE basis
          Teuchos::RCP<Intrepid::FieldContainer<Real > > valU_eval_bdry
            = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
          fe_bdry_[sidesets[i]][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
          Teuchos::RCP<Intrepid::FieldContainer<Real > > valL_eval_bdry
            = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
          fe_bdry_[sidesets[i]][j]->evaluateValue(valL_eval_bdry, l_coeff_bdry);
          // Compute Stefan-Boltzmann residual
          Teuchos::RCP< Intrepid::FieldContainer<Real> > sb_derivU
            = Teuchos::rcp( new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
          computeStefanBoltzmann(sb_derivU,valU_eval_bdry,sidesets[i],j,2);
          Intrepid::FieldContainer<Real> sb_derivU_L(numCellsSide, numCubPerSide);
          Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(sb_derivU_L,
                                                                     *sb_derivU,
                                                                     *valL_eval_bdry);
          Intrepid::FieldContainer<Real> sb_derivU_L_N(numCellsSide, f, numCubPerSide);
          Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(sb_derivU_L_N,
                                                                      sb_derivU_L,
                                                                      *(fe_bdry_[sidesets[i]][j]->N()));
          Intrepid::FieldContainer<Real> sbHess(numCellsSide, f, f);
          Intrepid::FunctionSpaceTools::integrate<Real>(sbHess,
                                                        sb_derivU_L_N,
                                                        *(fe_bdry_[sidesets[i]][j]->NdetJ()),
                                                        Intrepid::COMP_CPP, false);
          // Add Stefan-Boltzmann residual to volume residual
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[sidesets[i]][j][k];
            for (int l = 0; l < f; ++l) { 
              for (int m = 0; m < f; ++m) { 
                (*hess)(cidx,l,m) += sbHess(k,l,m);
              }
            }
          }
        }
      }
    }
  }

  void Hessian_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (StochasticStefanBoltzmannPDE::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (StochasticStefanBoltzmannPDE::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (StochasticStefanBoltzmannPDE::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & riesz) {
    riesz = fe_vol_->stiffMat();
    Intrepid::RealSpaceTools<Real>::add(*riesz,*(fe_vol_->massMat()));
  }

  void RieszMap_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & riesz) {
    riesz = fe_vol_->stiffMat();
    Intrepid::RealSpaceTools<Real>::add(*riesz,*(fe_vol_->massMat()));
  }
 
  std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > getFields(void) {
    return basisPtrs_;
  }

  void setCellNodes(const Teuchos::RCP<Intrepid::FieldContainer<Real> > &volCellNodes,
                    const std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > &bdryCellNodes, 
                    const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds ) {
    volCellNodes_   = volCellNodes;
    bdryCellNodes_  = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition
    fe_vol_ = Teuchos::rcp(new FE<Real>(volCellNodes_,basisPtr_,cellCub_));
    // Set local boundary DOFs
    fidx_ = fe_vol_->getBoundaryDofs();
    // Construct boundary FEs
    const int numSidesets = bdryCellNodes.size();
    fe_bdry_.resize(numSidesets);
    for(int i = 0; i < numSidesets; ++i) {
      int numLocSides = bdryCellNodes[i].size();
      fe_bdry_[i].resize(numLocSides);
      for (int j = 0; j < numLocSides; ++j) {
        if (bdryCellNodes_[i][j] != Teuchos::null) {
          fe_bdry_[i][j] = Teuchos::rcp(new FE<Real>(bdryCellNodes_[i][j],basisPtr_,bdryCub_,j));
        }
      }
    }
  }

  const Teuchos::RCP<FE<Real> > getFE(void) const {
    return fe_vol_;
  }
 
private:
    
  /***************************************************************************/
  /************** EVALUATE PDE COEFFICIENTS AT DOF COORDINATES ***************/
  /***************************************************************************/
  Real evaluateDiffusivity(const Real u, const std::vector<Real> & x, const int deriv = 0) const {
    if ( x[1] < xmid_ ) {
      const std::vector<Real> param = PDE<Real>::getParameter();
      if ( deriv > 0 ) {
        return static_cast<Real>(0);
      }
      return static_cast<Real>(1.43e-6) + static_cast<Real>(2.e-7)*param[0];
    }
    else {
      const std::vector<Real> param = PDE<Real>::getParameter();
      const Real c1 = static_cast<Real>(8.9e-5) + static_cast<Real>(1.e-6)*param[1];
      const Real c2 = static_cast<Real>(1.0e-4) + static_cast<Real>(1.e-5)*param[2];
      const Real c3 = static_cast<Real>(8.2e-3) + static_cast<Real>(3.e-4)*param[3];
      if ( deriv > 0 ) {
        return c2 * std::pow(-c3,deriv) * std::exp(-c3 * u);
      }
      return c2 * std::exp(-c3 * u) + c1;
    }
  }

  void evaluateVelocity(std::vector<Real> &adv, const std::vector<Real> &x) const {
    if ( x[1] < xmid_ ) {
      const std::vector<Real> param = PDE<Real>::getParameter();
      const Real quarter(0.25), half(0.5);
      const Real x1  = xmid_ * (quarter * param[4] + half);
      const Real det = (xmid_ - x1) * x1;
      const Real a   = static_cast<Real>(-1) / det;
      const Real b   = xmid_ / det;
      const Real m   = static_cast<Real>(1.e-6) + static_cast<Real>(1.e-7)*param[5];
      const Real mag = m * (a * std::pow(x[1],2) + b * x[1]);
      adv[0] = -mag;
      adv[1] = static_cast<Real>(0);
    }
    else {
      adv[0] = static_cast<Real>(0);
      adv[1] = static_cast<Real>(0);
    }
  }

  Real evaluateRHS(const std::vector<Real> &x) const {
    return static_cast<Real>(0);
  }

  Real evaluateStefanBoltzmann(const Real u, const std::vector<Real> &x,
                               const int sideset, const int locSideId,
                               const int deriv = 0) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    Real c1(0), c2(0), c3(0);
    if ( sideset == 2 ) {
      c1 = static_cast<Real>(1.0e-8) + static_cast<Real>(5.e-9)        * param[6];
      c2 = airTemp_                  + static_cast<Real>(0.1*airTemp_) * param[7];
      c3 = static_cast<Real>(1.4e-6) + static_cast<Real>(2.e-7)        * param[8];
    }
    else if ( sideset == 4 || sideset == 5 ) {
      c1 = static_cast<Real>(1.0e-8) + static_cast<Real>(5.e-9)        * param[9];
      c2 = engTemp_                  + static_cast<Real>(0.1*engTemp_) * param[10];
      c3 = static_cast<Real>(1.0e-4) + static_cast<Real>(1.e-5)        * param[11];
    }
    if ( deriv == 1 ) {
      return c1 * static_cast<Real>(4) * std::pow(u,3) + c3;
    }
    if ( deriv == 2 ) {
      return c1 * static_cast<Real>(4) * static_cast<Real>(3) * std::pow(u,2);
    }
    return c1 * (std::pow(u,4) - std::pow(c2,4)) + c3 * (u - c2);
  }
 
  /***************************************************************************/
  /************** COMPUTE PDE COEFFICIENTS AT DOFS ***************************/
  /***************************************************************************/
  void computeDiffusivity(Teuchos::RCP<Intrepid::FieldContainer<Real> > &kappa,
                          const Teuchos::RCP<Intrepid::FieldContainer<Real> > &u,
                          const int deriv = 0 ) const {
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
        // Compute diffusivity
        (*kappa)(i,j) = scale_*evaluateDiffusivity((*u)(i,j),pt,deriv);
      }
    }
  }

  void computeCoefficients(Teuchos::RCP<Intrepid::FieldContainer<Real> > &V,
                           Teuchos::RCP<Intrepid::FieldContainer<Real> > &rhs) const {
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
        // Compute advection velocity field V
        evaluateVelocity(adv,pt);
        for (int k = 0; k < d; ++k) {
          (*V)(i,j,k) = scale_*adv[k];
        }
        // Compute forcing term f
        (*rhs)(i,j) = -scale_*evaluateRHS(pt);
      }
    }
  }

  void computeStefanBoltzmann(Teuchos::RCP<Intrepid::FieldContainer<Real> > &sb,
                              const Teuchos::RCP<Intrepid::FieldContainer<Real> > &u,
                              const int sideset,
                              const int locSideId,
                              const int deriv = 0) const {
    const int c = u->dimension(0);
    const int p = u->dimension(1);
    const int d = fe_vol_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (*fe_bdry_[sideset][locSideId]->cubPts())(i,j,k);
        }
        (*sb)(i,j) = evaluateStefanBoltzmann((*u)(i,j),pt,sideset,locSideId,deriv);
      }
    }
  }

  /***************************************************************************/
  /************** EXTRACT COEFFICIENTS ON BOUNDARY ***************************/
  /***************************************************************************/
  Teuchos::RCP<Intrepid::FieldContainer<Real> > getBoundaryCoeff(
      const Intrepid::FieldContainer<Real> & cell_coeff,
      int sideSet, int cell) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][cell];
    const int numCellsSide = bdryCellLocId.size();
    const int f = basisPtr_->getCardinality();
    
    Teuchos::RCP<Intrepid::FieldContainer<Real > > bdry_coeff = 
      Teuchos::rcp(new Intrepid::FieldContainer<Real > (numCellsSide, f));
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        (*bdry_coeff)(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

}; // PDE_stefan_boltzmann

#endif
