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
    \brief Implements the local PDE interface for the Poisson-Boltzmann control problem.
*/

#ifndef PDE_POISSON_BOLTZMANN_EX02_HPP
#define PDE_POISSON_BOLTZMANN_EX02_HPP

#include "../TOOLS/pde.hpp"
#include "../TOOLS/fe.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "Teuchos_RCP.hpp"

template <class Real>
class PDE_Poisson_Boltzmann_ex02 : public PDE<Real> {
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

  Real wx_, wy_, dx_;
  Real uScale_, vScale_;
  bool useRobin_;
  Real robinCoeff_;

public:
  PDE_Poisson_Boltzmann_ex02(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("Problem").get("Order of FE discretization",1);
    if (basisOrder == 1) {
      basisPtr_ = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >);
    }
    else if (basisOrder == 2) {
      basisPtr_ = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >);
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

    wx_      = parlist.sublist("Geometry").get("Width",0.6);
    wy_      = parlist.sublist("Geometry").get("Height",0.2);
    dx_      = wx_/static_cast<Real>(6);
    uScale_  = parlist.sublist("Problem").get("Electron Density",1.0);
    vScale_  = parlist.sublist("Problem").get("Hole Density",1.0);

    useRobin_ = parlist.sublist("Problem").get("Use Robin Conditions",false);
    robinCoeff_ = parlist.sublist("Problem").get("Robin Coefficient",1e2);
  }

  void residual(Teuchos::RCP<Intrepid::FieldContainer<Real> > & res,
                const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN()->dimension(0);
    const int f = fe_vol_->gradN()->dimension(1);
    const int p = fe_vol_->gradN()->dimension(2);
    const int d = fe_vol_->gradN()->dimension(3);
    // INITIALIZE RESIDUAL
    res = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
    // COMPUTE PDE COEFFICIENTS
    Intrepid::FieldContainer<Real> lambda2(c, p), delta2(c, p), scale(c, p);
    computeCoefficients(lambda2,delta2,scale);
    // ADD STIFFNESS TERM TO RESIDUAL
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    Intrepid::FieldContainer<Real> lambda2_gradU_eval(c, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(lambda2_gradU_eval,
                                                               lambda2,
                                                               *gradU_eval);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  lambda2_gradU_eval,
                                                  *(fe_vol_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // ADD NONLINEAR TERM TO RESIDUAL
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    Intrepid::FieldContainer<Real> phi_valU_eval(c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        phi_valU_eval(i,j) = delta2(i,j)*(uScale_*std::exp((*valU_eval)(i,j)) - vScale_*std::exp(-(*valU_eval)(i,j)));
      }
    }
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  phi_valU_eval,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // ADD CONTROL TERM TO RESIDUAL
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_scal =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*valZ_scal,scale,*valZ_eval);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *valZ_scal,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // APPLY DIRICHLET CONDITIONS
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > bdryCellDofValues;
    computeDirichlet(bdryCellDofValues);
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      if ( i==1 || i==2 || i==3 ) {
        int numLocalSideIds = bdryCellLocIds_[i].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          int numCellsSide = bdryCellLocIds_[i][j].size();
          if (!useRobin_) {
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                (*res)(cidx,fidx_[j][l])
                  = (*u_coeff)(cidx,fidx_[j][l]) - (*bdryCellDofValues[i][j])(k,l);
              }
            }
          }
          else {
            const int numCubPerSide = bdryCub_->getNumPoints();
            if (numCellsSide) {
              // Get U and Z coefficients on Robin boundary
              Teuchos::RCP<Intrepid::FieldContainer<Real > > u_coeff_bdry, z_coeff_bdry;
              u_coeff_bdry = getBoundaryCoeff(*u_coeff, i, j);
              z_coeff_bdry = getBoundaryCoeff(*z_coeff, i, j);
              // Evaluate U and Z on FE basis
              Teuchos::RCP<Intrepid::FieldContainer<Real > > valU_eval_bdry, valZ_eval_bdry;
              valU_eval_bdry = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
              valZ_eval_bdry = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
              fe_bdry_[i][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
              fe_bdry_[i][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
              // Compute Robin residual
              Intrepid::FieldContainer<Real> robinRes(numCellsSide, f);
              Teuchos::RCP< Intrepid::FieldContainer<Real> > robinVal
                = Teuchos::rcp( new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
              computeRobin(robinVal,valU_eval_bdry,valZ_eval_bdry,i,j,0);
              Intrepid::FunctionSpaceTools::integrate<Real>(robinRes,
                                                            *robinVal,
                                                            *(fe_bdry_[i][j]->NdetJ()),
                                                            Intrepid::COMP_CPP, false);
              // Add Stefan-Boltzmann residual to volume residual
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
    }
  }

  void Jacobian_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & jac,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN()->dimension(0);
    const int f = fe_vol_->gradN()->dimension(1);
    const int p = fe_vol_->gradN()->dimension(2);
    const int d = fe_vol_->gradN()->dimension(3);
    // INITIALIZE JACOBIAN
    jac = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    // COMPUTE PDE COEFFICIENTS
    Intrepid::FieldContainer<Real> lambda2(c, p), delta2(c, p), scale(c, p);
    computeCoefficients(lambda2,delta2,scale);
    // ADD STIFFNESS TERM TO JACOBIAN
    Intrepid::FieldContainer<Real> lambda2_gradN_eval(c, f, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(lambda2_gradN_eval,
                                                                lambda2,
                                                                *(fe_vol_->gradN()));
    // COMPUTE STIFFNESS TERM
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  lambda2_gradN_eval,
                                                  *(fe_vol_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // ADD NONLINEAR TERM
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    Intrepid::FieldContainer<Real> dphi_valU_eval(c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        dphi_valU_eval(i,j) = delta2(i,j)*(uScale_*std::exp((*valU_eval)(i,j)) + vScale_*std::exp(-(*valU_eval)(i,j)));
      }
    }
    Intrepid::FieldContainer<Real> NexpU(c, f, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(NexpU,
                                                                dphi_valU_eval,
                                                                *(fe_vol_->N()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  NexpU,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // APPLY DIRICHLET CONDITIONS
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      if ( i==1 || i==2 || i==3 ) {
        int numLocalSideIds = bdryCellLocIds_[i].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          int numCellsSide = bdryCellLocIds_[i][j].size();
          if (!useRobin_) {
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
          else {
            const int numCubPerSide = bdryCub_->getNumPoints();
            if (numCellsSide) {
              // Get U and Z coefficients on Robin boundary
              Teuchos::RCP<Intrepid::FieldContainer<Real > > u_coeff_bdry, z_coeff_bdry;
              u_coeff_bdry = getBoundaryCoeff(*u_coeff, i, j);
              z_coeff_bdry = getBoundaryCoeff(*z_coeff, i, j);
              // Evaluate U and Z on FE basis
              Teuchos::RCP<Intrepid::FieldContainer<Real > > valU_eval_bdry, valZ_eval_bdry;
              valU_eval_bdry = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
              valZ_eval_bdry = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
              fe_bdry_[i][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
              fe_bdry_[i][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
              // Compute Robin residual
              Intrepid::FieldContainer<Real> robinVal_N(numCellsSide, f, numCubPerSide);
              Intrepid::FieldContainer<Real> robinJac(numCellsSide, f, f);
              Teuchos::RCP< Intrepid::FieldContainer<Real> > robinVal
                = Teuchos::rcp( new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
              computeRobin(robinVal,valU_eval_bdry,valZ_eval_bdry,i,j,1,1);
              Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(robinVal_N,
                                                                          *robinVal,
                                                                          *(fe_bdry_[i][j]->N()));
              Intrepid::FunctionSpaceTools::integrate<Real>(robinJac,
                                                            robinVal_N,
                                                            *(fe_bdry_[i][j]->NdetJ()),
                                                            Intrepid::COMP_CPP, false);
              // Add Stefan-Boltzmann residual to volume residual
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
    }
  }

  void Jacobian_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & jac,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN()->dimension(0);
    const int f = fe_vol_->gradN()->dimension(1);
    const int p = fe_vol_->gradN()->dimension(2);
    // INITIALIZE JACOBIAN
    jac = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    // COMPUTE PDE COEFFICIENTS
    Intrepid::FieldContainer<Real> lambda2(c, p), delta2(c, p), scale(c, p);
    computeCoefficients(lambda2,delta2,scale);
    // ADD CONTROL TERM
    Intrepid::FieldContainer<Real> valN_scal(c, f, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(valN_scal,scale,*(fe_vol_->N()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  valN_scal,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // APPLY DIRICHLET CONDITIONS
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      if ( i==1 || i==2 || i==3 ) {
        int numLocalSideIds = bdryCellLocIds_[i].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          int numCellsSide = bdryCellLocIds_[i][j].size();
          if (!useRobin_) {
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
          else {
            const int numCubPerSide = bdryCub_->getNumPoints();
            if (numCellsSide) {
              // Get U and Z coefficients on Robin boundary
              Teuchos::RCP<Intrepid::FieldContainer<Real > > u_coeff_bdry, z_coeff_bdry;
              u_coeff_bdry = getBoundaryCoeff(*u_coeff, i, j);
              z_coeff_bdry = getBoundaryCoeff(*z_coeff, i, j);
              // Evaluate U and Z on FE basis
              Teuchos::RCP<Intrepid::FieldContainer<Real > > valU_eval_bdry, valZ_eval_bdry;
              valU_eval_bdry = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
              valZ_eval_bdry = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
              fe_bdry_[i][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
              fe_bdry_[i][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
              // Compute Robin residual
              Intrepid::FieldContainer<Real> robinVal_N(numCellsSide, f, numCubPerSide);
              Intrepid::FieldContainer<Real> robinJac(numCellsSide, f, f);
              Teuchos::RCP< Intrepid::FieldContainer<Real> > robinVal
                = Teuchos::rcp( new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
              computeRobin(robinVal,valU_eval_bdry,valZ_eval_bdry,i,j,1,2);
              Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(robinVal_N,
                                                                          *robinVal,
                                                                          *(fe_bdry_[i][j]->N()));
              Intrepid::FunctionSpaceTools::integrate<Real>(robinJac,
                                                            robinVal_N,
                                                            *(fe_bdry_[i][j]->NdetJ()),
                                                            Intrepid::COMP_CPP, false);
              // Add Stefan-Boltzmann residual to volume residual
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
    }
  }

  void Hessian_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN()->dimension(0);
    const int f = fe_vol_->gradN()->dimension(1);
    const int p = fe_vol_->gradN()->dimension(2);
    // INITIALIZE HESSIAN
    hess = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    // COMPUTE PDE COEFFICIENTS
    Intrepid::FieldContainer<Real> lambda2(c, p), delta2(c, p), scale(c, p);
    computeCoefficients(lambda2,delta2,scale);
    // APPLY DIRICHLET CONDITIONS
    Teuchos::RCP<Intrepid::FieldContainer<Real> > l_dbc
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    *l_dbc = *l_coeff;
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      if ( i==1 || i==2 || i==3 ) {
        int numLocalSideIds = bdryCellLocIds_[i].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          int numCellsSide = bdryCellLocIds_[i][j].size();
          int numBdryDofs = fidx_[j].size();
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
            for (int l = 0; l < numBdryDofs; ++l) {
              (*l_dbc)(cidx,fidx_[j][l]) = static_cast<Real>(0);
            }
          }
        }
      }
    }
    // COMPUTE NONLINEAR TERM
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valL_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_vol_->evaluateValue(valL_eval, l_dbc);
    Intrepid::FieldContainer<Real> d2phi_valU_eval(c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        d2phi_valU_eval(i,j) = (*valL_eval)(i,j)*delta2(i,j)
                              *(uScale_*std::exp((*valU_eval)(i,j))-vScale_*std::exp(-(*valU_eval)(i,j)));
      }
    }
    Intrepid::FieldContainer<Real> NLexpU(c,f,p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(NLexpU,
                                                                d2phi_valU_eval,
                                                                *(fe_vol_->N()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  NLexpU,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void Hessian_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Poisson_Boltzmann_ex02:Hessian_12: Hessian is zero.");
  }

  void Hessian_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Poisson_Boltzmann_ex02:Hessian_21: Hessian is zero.");
  }

  void Hessian_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Poisson_Boltzmann_ex02:Hessian_22: Hessian is zero.");
  }

  void RieszMap_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & riesz) {
    // GET DIMENSIONS
    const int c = fe_vol_->N()->dimension(0);
    const int f = fe_vol_->N()->dimension(1);
    // INITIALIZE RIESZ
    riesz = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    *riesz = *fe_vol_->stiffMat();
    Intrepid::RealSpaceTools<Real>::add(*riesz,*(fe_vol_->massMat()));
  }

  void RieszMap_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & riesz) {
    // GET DIMENSIONS
    const int c = fe_vol_->N()->dimension(0);
    const int f = fe_vol_->N()->dimension(1);
    // INITIALIZE RIESZ
    riesz = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    *riesz = *fe_vol_->stiffMat();
    Intrepid::RealSpaceTools<Real>::add(*riesz,*(fe_vol_->massMat()));
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

  const std::vector<std::vector<Teuchos::RCP<FE<Real> > > > getBdryFE(void) const {
    return fe_bdry_;
  }

  const Teuchos::RCP<Intrepid::FieldContainer<Real> > getCellNodes(void) const {
    return volCellNodes_;
  }

  const std::vector<std::vector<std::vector<int> > > getBdryCellLocIds(void) const {
    return bdryCellLocIds_;
  }

private:

  Real evaluateLambda2(const std::vector<Real> &x) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    return static_cast<Real>(2.5)*std::pow(static_cast<Real>(10),param[0]);
  }

  Real evaluateDelta2(const std::vector<Real> &x) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    return static_cast<Real>(1.45)*std::pow(static_cast<Real>(10),param[1]);
  }

  Real evaluateScale(const std::vector<Real> &x) const {
    return 1;
//    const std::vector<Real> param = PDE<Real>::getParameter();
//    const Real Y1 = static_cast<Real>(1) + static_cast<Real>(5)*(wy_-x[1])*param[2];
//    const Real Y2 = static_cast<Real>(1) + static_cast<Real>(5)*(wy_-x[1])*param[3];
//    const Real X1 = phi(x[0]);
//    const Real X2 = static_cast<Real>(1)-X1;
//    return Y1*X1 + Y2*X2;
  }

  Real phi(const Real x) const {
    const Real zero(0), one(1), two(2);
    const Real eps = std::sqrt(ROL::ROL_EPSILON<Real>());
    const Real m = one/(wx_-two*dx_);
    const Real v = (x < dx_+eps ? zero
                   : (x > wx_-dx_-eps ? one
                     : m*(x-dx_)));
    return v;
  }

  void computeCoefficients(Intrepid::FieldContainer<Real> &lambda2,
                           Intrepid::FieldContainer<Real> &delta2,
                           Intrepid::FieldContainer<Real> &scale) const {
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
        // Compute scaled Debye length
        lambda2(i,j) = evaluateLambda2(pt);
        delta2(i,j)  = evaluateDelta2(pt);
        scale(i,j)   = -evaluateScale(pt);
      }
    }
  }

  Real evaluateRobin(const Real u, const Real z, const std::vector<Real> &x,
                     const int sideset, const int locSideId,
                     const int deriv = 0, const int component = 1) const {
    const Real four(4), two(2), one(1);
    Real C(0);
    if (sideset==1 || sideset==3) {
      C = static_cast<Real>(1);
    }
    if (sideset==2) {
      C = static_cast<Real>(0.3);
    }
    Real f0 = std::log((C+std::sqrt(C*C+four))/two);
    Real f1 = one/std::sqrt(C*C+four);
    if ( deriv == 1 ) {
      return (component==1) ? robinCoeff_ : -robinCoeff_ * f1;
    }
    if ( deriv > 1 ) {
      return static_cast<Real>(0);
    }
    return robinCoeff_ * (u - (f0 + f1*(z-C)));
  }

  void computeRobin(Teuchos::RCP<Intrepid::FieldContainer<Real> > &robin,
                    const Teuchos::RCP<Intrepid::FieldContainer<Real> > &u,
                    const Teuchos::RCP<Intrepid::FieldContainer<Real> > &z,
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

  Real evaluateDirichlet(const std::vector<Real> & coords, int sideset, int locSideId) const {
    Real val(0);
    if (sideset==1) {
      val = static_cast<Real>(0.436);
    }
    if (sideset==2) {
      val = static_cast<Real>(0.407);
    }
    if (sideset==3) {
      val = static_cast<Real>(0.436);
    }
    return val;
  }

  void computeDirichlet(std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > & bdryCellDofValues) const {
    // Compute Dirichlet values at DOFs.
    int d = basisPtr_->getBaseCellTopology().getDimension();
    int numSidesets = bdryCellLocIds_.size();
    bdryCellDofValues.resize(numSidesets);
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellDofValues[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
        int f = basisPtr_->getCardinality();
        bdryCellDofValues[i][j] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
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
            (*bdryCellDofValues[i][j])(k, l) = evaluateDirichlet(dofpoint, i, j);
          }
        }
      }
    }
  }

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

}; // PDE_Poisson_Boltzmann

template <class Real>
class PDE_Doping : public PDE<Real> {
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
  Real a_, b_;

public:
  PDE_Doping(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("Problem").get("Order of FE discretization",1);
    if (basisOrder == 1) {
      basisPtr_ = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >);
    }
    else if (basisOrder == 2) {
      basisPtr_ = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >);
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
    a_  = parlist.sublist("Problem").get("Desired Lower Doping Value", 0.3);
    b_  = parlist.sublist("Problem").get("Desired Upper Doping Value", 1.0);
  }

  void residual(Teuchos::RCP<Intrepid::FieldContainer<Real> > & res,
                const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    const std::vector<Real> param = PDE<Real>::getParameter();
    // GET DIMENSIONS
    const int c = fe_vol_->gradN()->dimension(0);
    const int f = fe_vol_->gradN()->dimension(1);
    const int p = fe_vol_->gradN()->dimension(2);
    const int d = fe_vol_->gradN()->dimension(3);
    // INITIALIZE RESIDUAL
    res = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
    // EVALUATE STATE AND CONTROL OF FEM BASIS
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valU_eval, valZ_eval, gradU_eval;
    valU_eval  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    valZ_eval  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    gradU_eval = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    // ADD STIFFNESS TERM TO RESIDUAL
    Intrepid::FieldContainer<Real> lambda2_gradU_eval(c, p, d);
    Intrepid::RealSpaceTools<Real>::scale(lambda2_gradU_eval,*gradU_eval,param[2]);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  lambda2_gradU_eval,
                                                  *(fe_vol_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // ADD REACTION TERM TO RESIDUAL
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *valU_eval,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // ADD CONTROL TERM TO RESIDUAL
    Intrepid::FieldContainer<Real> valZ_scal(c, p);
    Intrepid::RealSpaceTools<Real>::scale(valZ_scal,*valZ_eval,static_cast<Real>(-1));
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  valZ_scal,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // APPLY DIRICHLET CONDITIONS
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > bdryCellDofValues;
    computeDirichlet(bdryCellDofValues);
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      if ( i==1 || i==2 || i==3 ) {
        int numLocalSideIds = bdryCellLocIds_[i].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          int numCellsSide = bdryCellLocIds_[i][j].size();
          int numBdryDofs = fidx_[j].size();
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
            for (int l = 0; l < numBdryDofs; ++l) {
              (*res)(cidx,fidx_[j][l]) = (*u_coeff)(cidx,fidx_[j][l]) - (*bdryCellDofValues[i][j])(k,l);
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
    const std::vector<Real> param = PDE<Real>::getParameter();
    // GET DIMENSIONS
    const int c = fe_vol_->gradN()->dimension(0);
    const int f = fe_vol_->gradN()->dimension(1);
    const int p = fe_vol_->gradN()->dimension(2);
    const int d = fe_vol_->gradN()->dimension(3);
    // INITIALIZE JACOBIAN
    jac = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    // ADD STIFFNESS TERM TO JACOBIAN
    Intrepid::FieldContainer<Real> lambda2_gradN_eval(c, f, p, d);
    Intrepid::RealSpaceTools<Real>::scale(lambda2_gradN_eval,*(fe_vol_->gradN()),param[2]);
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  lambda2_gradN_eval,
                                                  *(fe_vol_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // ADD REACTION TERM TO JACOBIAN
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *(fe_vol_->N()),
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    // APPLY DIRICHLET CONDITIONS
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      if ( i==1 || i==2 || i==3 ) {
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
    // GET DIMENSIONS
    const int c = fe_vol_->gradN()->dimension(0);
    const int f = fe_vol_->gradN()->dimension(1);
    const int p = fe_vol_->gradN()->dimension(2);
    // INITIALIZE JACOBIAN
    jac = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    // ADD CONTROL TERM
    Intrepid::FieldContainer<Real> valN_scal(c, f, p);
    Intrepid::RealSpaceTools<Real>::scale(valN_scal,*(fe_vol_->N()),static_cast<Real>(-1));
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  valN_scal,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // APPLY DIRICHLET CONDITIONS
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      if ( i==1 || i==2 || i==3 ) {
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

  void Hessian_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Poisson_Boltzmann_ex02:Hessian_11: Hessian is zero.");
  }

  void Hessian_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Poisson_Boltzmann_ex02:Hessian_12: Hessian is zero.");
  }

  void Hessian_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Poisson_Boltzmann_ex02:Hessian_21: Hessian is zero.");
  }

  void Hessian_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Poisson_Boltzmann_ex02:Hessian_22: Hessian is zero.");
  }

  void RieszMap_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & riesz) {
    // GET DIMENSIONS
    const int c = fe_vol_->N()->dimension(0);
    const int f = fe_vol_->N()->dimension(1);
    // INITIALIZE RIESZ
    riesz = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    *riesz = *fe_vol_->stiffMat();
    Intrepid::RealSpaceTools<Real>::add(*riesz,*(fe_vol_->massMat()));
  }

  void RieszMap_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & riesz) {
    // GET DIMENSIONS
    const int c = fe_vol_->N()->dimension(0);
    const int f = fe_vol_->N()->dimension(1);
    // INITIALIZE RIESZ
    riesz = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    *riesz = *fe_vol_->stiffMat();
    Intrepid::RealSpaceTools<Real>::add(*riesz,*(fe_vol_->massMat()));
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

  const std::vector<std::vector<Teuchos::RCP<FE<Real> > > > getBdryFE(void) const {
    return fe_bdry_;
  }

  const Teuchos::RCP<Intrepid::FieldContainer<Real> > getCellNodes(void) const {
    return volCellNodes_;
  }

  const std::vector<std::vector<std::vector<int> > > getBdryCellLocIds(void) const {
    return bdryCellLocIds_;
  }

private:

  Real evaluateDirichlet(const std::vector<Real> & coords, int sideset, int locSideId) const {
    Real val(0);
    if (sideset==1) {
      val = a_;
    }
    if (sideset==2) {
      val = b_;
    }
    if (sideset==3) {
      val = a_;
    }
    return val;
  }

  void computeDirichlet(std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > & bdryCellDofValues) const {
    // Compute Dirichlet values at DOFs.
    int d = basisPtr_->getBaseCellTopology().getDimension();
    int numSidesets = bdryCellLocIds_.size();
    bdryCellDofValues.resize(numSidesets);
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellDofValues[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
        int f = basisPtr_->getCardinality();
        bdryCellDofValues[i][j] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
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
            (*bdryCellDofValues[i][j])(k, l) = evaluateDirichlet(dofpoint, i, j);
          }
        }
      }
    }
  }

}; // PDE_Poisson_Boltzmann

#endif
