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

/*! \file  pde_ginzburg-landau.hpp
    \brief Implements the local PDE interface for the optimal control of
           simplified Ginzburg-Landau.
*/

#ifndef PDE_GINZBURGLANDAU_HPP
#define PDE_GINZBURGLANDAU_HPP

#include "../TOOLS/pde.hpp"
#include "../TOOLS/fe.hpp"
#include "../TOOLS/fieldhelper.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "Teuchos_RCP.hpp"


template <class Real>
class PDE_GinzburgLandau : public PDE<Real> {
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
  Teuchos::RCP<FE<Real> > fe_;
  std::vector<Teuchos::RCP<FE<Real> > > feBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int> > fidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > bdryCellDofValues_;
  // Field pattern, offsets, etc.
  std::vector<std::vector<int> > fieldPattern_;  // local Field/DOF pattern; set from DOF manager 
  int numFields_;                                // number of fields (equations in the PDE)
  int numDofs_;                                  // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                      // for each field, a counting offset
  std::vector<int> numFieldDofs_;                // for each field, number of degrees of freedom

  Real lambda_;

  Teuchos::RCP<FieldHelper<Real> > fieldHelper_;

  void computeMagneticPotential(const Teuchos::RCP<Intrepid::FieldContainer<Real> > &A) const {
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
   
    std::vector<Real> x(d), Ax(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (*fe_->cubPts())(i,j,k);
        }
        evaluateMagneticPotential(Ax,x);
        for (int k = 0; k < d; ++k) {
          (*A)(i,j,k) = Ax[k];
        }
      }
    } 
  }

  void computeForce(const Teuchos::RCP<Intrepid::FieldContainer<Real> > &F, const int component) const {
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
   
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (*fe_->cubPts())(i,j,k);
        }
        (*F)(i,j) = -evaluateForce(x,component);
      }
    } 
  }

  void computeNeumann(Teuchos::RCP<Intrepid::FieldContainer<Real> > &neumann,
                      const int locSideId,                              
                      const int component) const {
    const int c = feBdry_[locSideId]->gradN()->dimension(0);
    const int p = feBdry_[locSideId]->gradN()->dimension(2);
    const int d = fe_->gradN()->dimension(3);
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (*feBdry_[locSideId]->cubPts())(i,j,k);
        }
        (*neumann)(i,j) = evaluateNeumann(x,component)/lambda_;
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

public:
  PDE_GinzburgLandau(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    basisPtr_ = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();            // get the cell type from the basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                           // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 2);        // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                           // create default cubature
    int d = cellType.getDimension();

    basisPtrs_.clear();
    for (int i=0; i<2; ++i) {
      basisPtrs_.push_back(basisPtr_);  // Displacement component
    }

    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2); // set cubature degree, e.g., 2
    bdryCub_ = cubFactory.create(bdryCellType, bdryCubDegree);

    numDofs_ = 0;
    numFields_ = basisPtrs_.size();
    offset_.resize(numFields_);
    numFieldDofs_.resize(numFields_);
    for (int i=0; i<numFields_; ++i) {
      if (i==0) {
        offset_[i]  = 0;
      }
      else {
        offset_[i]  = offset_[i-1] + basisPtrs_[i-1]->getCardinality();
      }
      numFieldDofs_[i] = basisPtrs_[i]->getCardinality();
      numDofs_ += numFieldDofs_[i];
    }

    lambda_ = parlist.sublist("Problem").get("Current Loading",1.0);
  }

  virtual void evaluateMagneticPotential(std::vector<Real> &Ax, const std::vector<Real> &x) const = 0;

  virtual Real evaluateNeumann(const std::vector<Real> &x, const int component) const = 0;

  virtual Real evaluateForce(const std::vector<Real> &x, const int component) const = 0;

  void residual(Teuchos::RCP<Intrepid::FieldContainer<Real> > & res,
                const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
 
    // Initialize residuals.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > R(2);
    for (int i=0; i<2; ++i) {
      R[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f));
    }

    // Split u_coeff into components.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > U;
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > valU_eval(2);
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > gradU_eval(2);
    for (int i=0; i<2; ++i) {
      valU_eval[i]  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      fe_->evaluateValue(valU_eval[i], U[i]);
      gradU_eval[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
      fe_->evaluateGradient(gradU_eval[i], U[i]);
    }

    // Build force term
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > F(2);
    for (int i=0; i<2; ++i) {
      F[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      computeForce(F[i],i);
    }

    // Build magnetic potential
    Teuchos::RCP<Intrepid::FieldContainer<Real> > A
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    computeMagneticPotential(A);

    // Compute the magnitude of A
    Teuchos::RCP<Intrepid::FieldContainer<Real> > magA
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Intrepid::RealSpaceTools<Real>::dot(*magA,*A,*A);

    // Compute magnitude of U
    Teuchos::RCP<Intrepid::FieldContainer<Real> > magU
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > sqrU1
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*magU, *valU_eval[0],*valU_eval[0]);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*sqrU1,*valU_eval[1],*valU_eval[1]);
    Intrepid::RealSpaceTools<Real>::add(*magU,*sqrU1);
    Intrepid::RealSpaceTools<Real>::scale(*magU,lambda_);

    // Temporaries for residual evaluation
    Teuchos::RCP<Intrepid::FieldContainer<Real> > magAU
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > AgradU
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > AU
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > magUU
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    /*** Evaluate weak form of the residual. ***/
    int index = 0;
    for (int i = 0; i < 2; ++i) {
      index = (i+1)%2;
      // Symmetric term
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*magAU,*magA,*valU_eval[i]);
      Intrepid::RealSpaceTools<Real>::subtract(*magAU,*valU_eval[i]);
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *gradU_eval[i],     // grad U
                                                    *fe_->gradNdetJ(),  // grad N
                                                    Intrepid::COMP_CPP,
                                                    false);
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *magAU,             // (|A|^2 - 1)U
                                                    *fe_->NdetJ(),      // N
                                                    Intrepid::COMP_CPP,
                                                    true);
      // Nonlinear term
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*magUU,*magU,*valU_eval[i]);
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *magUU,             // |U|^2 U
                                                    *fe_->NdetJ(),      // N
                                                    Intrepid::COMP_CPP,
                                                    true);
      // Antisymmetric term
      Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(*AgradU,*A,*gradU_eval[index]);
      if (index) {
        Intrepid::RealSpaceTools<Real>::scale(*AgradU,static_cast<Real>(-1));
      }
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*AU,*valU_eval[index],*A);
      if (!index) {
        Intrepid::RealSpaceTools<Real>::scale(*AU,static_cast<Real>(-1));
      }
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *AgradU,            // -A . grad U
                                                    *fe_->NdetJ(),      // N
                                                    Intrepid::COMP_CPP,
                                                    true);
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *AU,                // A U
                                                    *fe_->gradNdetJ(),  // grad N
                                                    Intrepid::COMP_CPP,
                                                    true);
      // Force term
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *F[i],              // -F
                                                    *fe_->NdetJ(),      // N
                                                    Intrepid::COMP_CPP,
                                                    true);
    }

    // APPLY ROBIN CONTROLS: Sideset 0
    int sideset = 0;
    int numLocalSideIds = bdryCellLocIds_[sideset].size();
    const int numCubPerSide = bdryCub_->getNumPoints();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[sideset][j].size();
      if (numCellsSide) {
        for (int i = 0; i < 2; ++i) {
          // Get U coefficients on Robin boundary
          Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff_bdry
            = getBoundaryCoeff(*Z[i], sideset, j);
          // Evaluate U on FE basis
          Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval_bdry
            = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
          feBdry_[j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
          // Compute uncontrolled Neumann source
          Teuchos::RCP<Intrepid::FieldContainer<Real> > neumann
            = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
          computeNeumann(neumann,j,i);
          // Add uncontrolled Neumann source to control
          Intrepid::RealSpaceTools<Real>::add(*valZ_eval_bdry, *neumann);
          // Compute Neumann residual
          Intrepid::FieldContainer<Real> neumannRes(numCellsSide, f);
          Intrepid::FunctionSpaceTools::integrate<Real>(neumannRes,
                                                        *valZ_eval_bdry,
                                                        *(feBdry_[j]->NdetJ()),
                                                        Intrepid::COMP_CPP, false);
          // Add Neumann control residual to volume residual
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[sideset][j][k];
            for (int l = 0; l < f; ++l) { 
              (*R[i])(cidx,l) -= lambda_*neumannRes(k,l);
            }
          }
        }
      }
    }

    // Combine the residuals.
    fieldHelper_->combineFieldCoeff(res, R);
  }

  void Jacobian_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & jac,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
 
    // Initialize Jacobians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(2);
    for (int i=0; i<2; ++i) {
      for (int j=0; j<2; ++j) {
        J[i].push_back(Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,f)));
      }
    }

    // Split u_coeff into components.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > valU_eval(2);
    for (int i=0; i<2; ++i) {
      valU_eval[i]  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      fe_->evaluateValue(valU_eval[i], U[i]);
    }

    // Build magnetic potential
    Teuchos::RCP<Intrepid::FieldContainer<Real> > A
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    computeMagneticPotential(A);

    // Compute magnitude of magnetic potential
    Teuchos::RCP<Intrepid::FieldContainer<Real> > magA
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Intrepid::RealSpaceTools<Real>::dot(*magA,*A,*A);

    // Multiply magnitude of magnetic potential with basis function
    Teuchos::RCP<Intrepid::FieldContainer<Real> > magAN
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, p));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*magAN,*magA,*fe_->N());
    Intrepid::RealSpaceTools<Real>::subtract(*magAN,*fe_->N());

    // Dot magnetic potential with gradient of basis function
    Teuchos::RCP<Intrepid::FieldContainer<Real> > AgradN
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, p));
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*AgradN,*A,*fe_->gradN());

    // Compute jacobian of nonlinearity
    Teuchos::RCP<Intrepid::FieldContainer<Real> > sqrU
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > sqrUN(2);
    for (int i = 0; i < 2; ++i) {
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*sqrU,*valU_eval[i],*valU_eval[i]);
      Intrepid::RealSpaceTools<Real>::scale(*sqrU,lambda_);
      sqrUN[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, p));
      Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*sqrUN[i],*sqrU,*fe_->N());
    }
    Teuchos::RCP<Intrepid::FieldContainer<Real> > U0U1
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > U0U1N
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, p));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*U0U1,*valU_eval[0],*valU_eval[1]);
    Intrepid::RealSpaceTools<Real>::scale(*U0U1,static_cast<Real>(2)*lambda_);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*U0U1N,*U0U1,*fe_->N());

    // Temporaries for residual evaluation
    Teuchos::RCP<Intrepid::FieldContainer<Real> > sqrUN3
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, p));
    /*** Evaluate weak form of the Jacobian. ***/
    for (int i=0; i<2; ++i) {
      int index = (i+1)%2;
      // Diagonal linear terms
      *J[i][i] = *fe_->stiffMat();
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][i],
                                                    *magAN,             // (|A|^2-1) N
                                                    *fe_->NdetJ(),      // N
                                                    Intrepid::COMP_CPP,
                                                    true);
      // Diagonal nonlinear terms
      Intrepid::RealSpaceTools<Real>::scale(*sqrUN3,*sqrUN[i],static_cast<Real>(3));
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][i],
                                                    *sqrUN3,           // 3 U^2 N
                                                    *fe_->NdetJ(),     // N
                                                    Intrepid::COMP_CPP,
                                                    true);
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][i],
                                                    *sqrUN[index],     // U^2 N
                                                    *fe_->NdetJ(),     // N
                                                    Intrepid::COMP_CPP,
                                                    true);
      // Off-diagonal nonlinear terms
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][index],
                                                    *U0U1N,             // U U N
                                                    *fe_->NdetJ(),      // N
                                                    Intrepid::COMP_CPP,
                                                    false);
      
    }
    // Off-diagonal linear terms
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][1],
                                                  *AgradN,            // A . grad N
                                                  *fe_->NdetJ(),      // N
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][0],
                                                  *fe_->NdetJ(),      // N
                                                  *AgradN,            // A . grad N
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::RealSpaceTools<Real>::scale(*AgradN,static_cast<Real>(-1));
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][0],
                                                  *AgradN,            // -A . grad N
                                                  *fe_->NdetJ(),      // N
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][1],
                                                  *fe_->NdetJ(),      // N
                                                  *AgradN,            // -A . grad N
                                                  Intrepid::COMP_CPP,
                                                  true);

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);

  }


  void Jacobian_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & jac,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);

    // Initialize Jacobians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(2);
    for (int i=0; i<2; ++i) {
      for (int j=0; j<2; ++j) {
        J[i].push_back(Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,f)));
      }
    }

    // APPLY ROBIN CONTROL: Sideset 0
    int sideset = 0;
    int numLocalSideIds = bdryCellLocIds_[sideset].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[sideset][j].size();
      if (numCellsSide) {
        for (int i = 0; i < 2; ++i) {
          // Add Neumann control Jacobian to volume residual
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[sideset][j][k];
            for (int l = 0; l < f; ++l) { 
              for (int m = 0; m < f; ++m) { 
                (*J[i][i])(cidx,l,m) -= lambda_*(*feBdry_[j]->massMat())(k,l,m);
              }
            }
          }
        }
      }
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);

  }

  void Hessian_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
 
    // Initialize Jacobians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > H(2);
    for (int i=0; i<2; ++i) {
      for (int j=0; j<2; ++j) {
        H[i].push_back(Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,f)));
      }
    }

    // Split u_coeff into components.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);

    // Split l_coeff into components.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > L;
    fieldHelper_->splitFieldCoeff(L, l_coeff);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > valU_eval(2);
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > valL_eval(2);
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > UL(2);
    for (int i=0; i<2; ++i) {
      valU_eval[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      fe_->evaluateValue(valU_eval[i], U[i]);
      valL_eval[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      fe_->evaluateValue(valL_eval[i], L[i]);
      UL[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*UL[i],*valU_eval[i],*valL_eval[i]);
    }

    Teuchos::RCP<Intrepid::FieldContainer<Real> > U0L1_U1L0
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > U1L0
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > U0L1_U1L0_N
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,p));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*U0L1_U1L0,*valU_eval[0],*valL_eval[1]);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*U1L0,     *valU_eval[1],*valL_eval[0]);
    Intrepid::RealSpaceTools<Real>::add(*U0L1_U1L0,*U1L0);
    Intrepid::RealSpaceTools<Real>::scale(*U0L1_U1L0,static_cast<Real>(2)*lambda_);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*U0L1_U1L0_N,*U0L1_U1L0,*fe_->N());
    
    Teuchos::RCP<Intrepid::FieldContainer<Real> > diag
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > diag_N
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,p));

    for (int i = 0; i < 2; ++i) {
      int index = (i+1)%2;

      Intrepid::RealSpaceTools<Real>::scale(*diag,*UL[i],static_cast<Real>(3));
      Intrepid::RealSpaceTools<Real>::add(*diag,*UL[index]);
      Intrepid::RealSpaceTools<Real>::scale(*diag,static_cast<Real>(2)*lambda_);
      Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*diag_N,*diag,*fe_->N());
      Intrepid::FunctionSpaceTools::integrate<Real>(*H[i][i],
                                                    *diag_N,            // (6 U0 L0 + 2 U1 L1) N
                                                    *fe_->NdetJ(),      // N
                                                    Intrepid::COMP_CPP,
                                                    false);

      Intrepid::FunctionSpaceTools::integrate<Real>(*H[i][index],
                                                    *U0L1_U1L0_N,       // 2(U1 L0 + U0 L1) N
                                                    *fe_->NdetJ(),      // N
                                                    Intrepid::COMP_CPP,
                                                    false);
    }

    // Combine the hessians.
    fieldHelper_->combineFieldCoeff(hess, H);

  }

  void Hessian_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_GinzburgLandau::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_GinzburgLandau::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_GinzburgLandau::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & riesz) {
    //throw Exception::NotImplemented(">>> (PDE_TopoOpt::RieszMap_1): Not implemented.");

    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int d = fe_->gradN()->dimension(3);
 
    // Initialize Jacobians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(d);
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j) {
        J[i].push_back(Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,f)));
      }
    }

    for (int i=0; i<d; ++i) {
      *(J[i][i]) = *(fe_->stiffMat());
      Intrepid::RealSpaceTools<Real>::add(*(J[i][i]),*(fe_->massMat()));
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  void RieszMap_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & riesz) {
    //throw Exception::NotImplemented(">>> (PDE_TopoOpt::RieszMap_2): Not implemented.");

    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int d = fe_->gradN()->dimension(3);
 
    // Initialize Jacobians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(d);
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j) {
        J[i].push_back(Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,f)));
      }
    }

    for (int i=0; i<d; ++i) {
      *(J[i][i]) = *(fe_->massMat());
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
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
    fe_ = Teuchos::rcp(new FE<Real>(volCellNodes_,basisPtr_,cellCub_));
    fidx_ = fe_->getBoundaryDofs();
    // Construct boundary FE
    int sideset = 0;
    int numLocSides = bdryCellNodes[sideset].size();
    feBdry_.resize(numLocSides);
    for (int j = 0; j < numLocSides; ++j) {
      if (bdryCellNodes[sideset][j] != Teuchos::null) {
        feBdry_[j] = Teuchos::rcp(new FE<Real>(bdryCellNodes[sideset][j],basisPtr_,bdryCub_,j));
      }
    }
  }

  void setFieldPattern(const std::vector<std::vector<int> > & fieldPattern) {
    fieldPattern_ = fieldPattern;
    fieldHelper_ = Teuchos::rcp(new FieldHelper<Real>(numFields_, numDofs_, numFieldDofs_, fieldPattern_));
  }

  const Teuchos::RCP<FE<Real> > getFE(void) const {
    return fe_;
  }

  const std::vector<Teuchos::RCP<FE<Real> > > getBdryFE(void) const {
    return feBdry_;
  }

  const std::vector<std::vector<int> > getBdryCellLocIds(const int sideset = 0) const {
    return bdryCellLocIds_[sideset];
  }

  const Teuchos::RCP<FieldHelper<Real> > getFieldHelper(void) const {
    return fieldHelper_;
  }

}; // PDE_GinzburgLandau


#endif
