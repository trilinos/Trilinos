// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#include "ROL_Ptr.hpp"


template <class Real>
class PDE_GinzburgLandau : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>> basisPtr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real>> cellCub_;
  ROL::Ptr<Intrepid::Cubature<Real>> bdryCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real>> volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real>> fe_;
  std::vector<ROL::Ptr<FE<Real>>> feBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellDofValues_;
  // Field pattern, offsets, etc.
  std::vector<std::vector<int>> fieldPattern_;  // local Field/DOF pattern; set from DOF manager 
  int numFields_;                               // number of fields (equations in the PDE)
  int numDofs_;                                 // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                     // for each field, a counting offset
  std::vector<int> numFieldDofs_;               // for each field, number of degrees of freedom

  Real lambda_;

  ROL::Ptr<FieldHelper<Real>> fieldHelper_;

  void computeMagneticPotential(const ROL::Ptr<Intrepid::FieldContainer<Real>> &A) const {
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

  void computeForce(const ROL::Ptr<Intrepid::FieldContainer<Real>> &F, const int component) const {
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

  void computeNeumann(ROL::Ptr<Intrepid::FieldContainer<Real>> &neumann,
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

  ROL::Ptr<Intrepid::FieldContainer<Real>> getBoundaryCoeff(
      const Intrepid::FieldContainer<Real> &cell_coeff,
      int sideSet, int cell) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][cell];
    const int numCellsSide = bdryCellLocId.size();
    const int f = basisPtr_->getCardinality();
    
    ROL::Ptr<Intrepid::FieldContainer<Real>> bdry_coeff = 
      ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f);
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
    basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real,Intrepid::FieldContainer<Real>>>();
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
      if (i==0) offset_[i]  = 0;
      else      offset_[i]  = offset_[i-1] + basisPtrs_[i-1]->getCardinality();
      numFieldDofs_[i] = basisPtrs_[i]->getCardinality();
      numDofs_ += numFieldDofs_[i];
    }

    lambda_ = parlist.sublist("Problem").get("Current Loading",1.0);
  }

  virtual void evaluateMagneticPotential(std::vector<Real> &Ax, const std::vector<Real> &x) const = 0;

  virtual Real evaluateNeumann(const std::vector<Real> &x, const int component) const = 0;

  virtual Real evaluateForce(const std::vector<Real> &x, const int component) const = 0;

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real>>             &res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    const int d = fe_->gradN()->dimension(3);
 
    // Initialize residuals.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> R(2);
    R[0] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f);
    R[1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f);

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U, Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> valU_eval(2), gradU_eval(2), F(2);
    ROL::Ptr<Intrepid::FieldContainer<Real>> A, magA, magU, sqrU1, magAU, AgradU, AU, magUU;
    A      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    magA   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    magU   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    sqrU1  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    magAU  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    AgradU = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    AU     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    magUU  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i=0; i<2; ++i) {
      valU_eval[i]  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      gradU_eval[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      F[i]          = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
       // Evaluate/interpolate finite element fields on cells.
      fe_->evaluateValue(valU_eval[i], U[i]);
      fe_->evaluateGradient(gradU_eval[i], U[i]);
      // Build force term
      computeForce(F[i],i);
    }
    // Build magnetic potential
    computeMagneticPotential(A);
    // Compute the magnitude of A
    Intrepid::RealSpaceTools<Real>::dot(*magA,*A,*A);
    // Compute magnitude of U
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*magU, *valU_eval[0],*valU_eval[0]);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*sqrU1,*valU_eval[1],*valU_eval[1]);
    Intrepid::RealSpaceTools<Real>::add(*magU,*sqrU1);
    Intrepid::RealSpaceTools<Real>::scale(*magU,lambda_);

    /*** Evaluate weak form of the residual. ***/
    int index = 0;
    for (int i = 0; i < 2; ++i) {
      index = (i+1)%2;
      // Symmetric term
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *gradU_eval[i],     // grad U
                                                    *fe_->gradNdetJ(),  // grad N
                                                    Intrepid::COMP_CPP,
                                                    false);
      magAU->initialize();
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*magAU,*magA,*valU_eval[i]);
      Intrepid::RealSpaceTools<Real>::subtract(*magAU,*valU_eval[i]);
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *magAU,             // (|A|^2 - 1)U
                                                    *fe_->NdetJ(),      // N
                                                    Intrepid::COMP_CPP,
                                                    true);
      // Nonlinear term
      magUU->initialize();
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*magUU,*magU,*valU_eval[i]);
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *magUU,             // |U|^2 U
                                                    *fe_->NdetJ(),      // N
                                                    Intrepid::COMP_CPP,
                                                    true);
      // Antisymmetric term
      AgradU->initialize();
      Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(*AgradU,*A,*gradU_eval[index]);
      if (index==1) Intrepid::RealSpaceTools<Real>::scale(*AgradU,static_cast<Real>(-1));
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *AgradU,            // -A . grad U
                                                    *fe_->NdetJ(),      // N
                                                    Intrepid::COMP_CPP,
                                                    true);
      AU->initialize();
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*AU,*valU_eval[index],*A);
      if (index==0) Intrepid::RealSpaceTools<Real>::scale(*AU,static_cast<Real>(-1));
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
          ROL::Ptr<Intrepid::FieldContainer<Real>> z_coeff_bdry
            = getBoundaryCoeff(*Z[i], sideset, j);
          // Evaluate U on FE basis
          ROL::Ptr<Intrepid::FieldContainer<Real>> valZ_eval_bdry
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          feBdry_[j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
          // Compute uncontrolled Neumann source
          ROL::Ptr<Intrepid::FieldContainer<Real>> neumann
            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
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

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real>>             &jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    const int d = fe_->gradN()->dimension(3);
 
    // Initialize Jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J(2);
    for (int i=0; i<2; ++i) {
      for (int j=0; j<2; ++j) {
        J[i].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f));
      }
    }

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> valU_eval(2), sqrUN(2);
    ROL::Ptr<Intrepid::FieldContainer<Real>> A, magA, magAN, AgradN, sqrU, U0U1, U0U1N, sqrUN3;
    A      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    magA   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    magAN  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    AgradN = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    sqrU   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    U0U1   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    U0U1N  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    sqrUN3 = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    for (int i=0; i<2; ++i) {
      valU_eval[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      fe_->evaluateValue(valU_eval[i], U[i]);
      // Compute jacobian of nonlinearity
      sqrU->initialize();
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*sqrU,*valU_eval[i],*valU_eval[i]);
      Intrepid::RealSpaceTools<Real>::scale(*sqrU,lambda_);
      sqrUN[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
      Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*sqrUN[i],*sqrU,*fe_->N());
    }
    // Build magnetic potential
    computeMagneticPotential(A);
    // Compute magnitude of magnetic potential
    Intrepid::RealSpaceTools<Real>::dot(*magA,*A,*A);
    // Multiply magnitude of magnetic potential with basis function
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*magAN,*magA,*fe_->N());
    Intrepid::RealSpaceTools<Real>::subtract(*magAN,*fe_->N());
    // Dot magnetic potential with gradient of basis function
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*AgradN,*A,*fe_->gradN());
    // Compute jacobian of nonlinearity
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*U0U1,*valU_eval[0],*valU_eval[1]);
    Intrepid::RealSpaceTools<Real>::scale(*U0U1,static_cast<Real>(2)*lambda_);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*U0U1N,*U0U1,*fe_->N());

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
      sqrUN3->initialize();
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
      
      // Off-diagonal linear terms
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][index],
                                                    *AgradN,            // A . grad N
                                                    *fe_->NdetJ(),      // N
                                                    Intrepid::COMP_CPP,
                                                    true);
      Intrepid::RealSpaceTools<Real>::scale(*AgradN,static_cast<Real>(-1));
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][index],
                                                    *fe_->NdetJ(),      // N
                                                    *AgradN,            // -A . grad N
                                                    Intrepid::COMP_CPP,
                                                    true);
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }


  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real>>             &jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);

    // Initialize Jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J(2);
    for (int i=0; i<2; ++i) {
      for (int j=0; j<2; ++j) {
        J[i].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f));
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

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real>>             &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
 
    // Initialize Jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> H(2);
    for (int i=0; i<2; ++i) {
      for (int j=0; j<2; ++j) {
        H[i].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f));
      }
    }

    // Split u_coeff and l_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U, L;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(L, l_coeff);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> valU_eval(2), valL_eval(2), UL(2);
    ROL::Ptr<Intrepid::FieldContainer<Real>> U0L1_U1L0, U1L0, U0L1_U1L0_N, diag, diag_N;
    U0L1_U1L0   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    U1L0        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    U0L1_U1L0_N = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,p);
    diag        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    diag_N      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,p);
    for (int i=0; i<2; ++i) {
      valU_eval[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      fe_->evaluateValue(valU_eval[i], U[i]);
      valL_eval[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      fe_->evaluateValue(valL_eval[i], L[i]);
      UL[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*UL[i],*valU_eval[i],*valL_eval[i]);
    }
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*U0L1_U1L0,*valU_eval[0],*valL_eval[1]);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*U1L0,     *valU_eval[1],*valL_eval[0]);
    Intrepid::RealSpaceTools<Real>::add(*U0L1_U1L0,*U1L0);
    Intrepid::RealSpaceTools<Real>::scale(*U0L1_U1L0,static_cast<Real>(2)*lambda_);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*U0L1_U1L0_N,*U0L1_U1L0,*fe_->N());

    for (int i = 0; i < 2; ++i) {
      int index = (i+1)%2;
      diag->initialize(); diag_N->initialize();
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

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real>>             &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_GinzburgLandau::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real>>             &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_GinzburgLandau::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real>>             &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_GinzburgLandau::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    //throw Exception::NotImplemented(">>> (PDE_TopoOpt::RieszMap_1): Not implemented.");

    // Retrieve dimensions.
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
 
    // Initialize Jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J(2);
    for (int i=0; i<2; ++i) {
      for (int j=0; j<2; ++j) {
        J[i].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f));
      }
    }

    for (int i=0; i<2; ++i) {
      *(J[i][i]) = *(fe_->stiffMat());
      Intrepid::RealSpaceTools<Real>::add(*(J[i][i]),*(fe_->massMat()));
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    //throw Exception::NotImplemented(">>> (PDE_TopoOpt::RieszMap_2): Not implemented.");

    // Retrieve dimensions.
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
 
    // Initialize Jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J(2);
    for (int i=0; i<2; ++i) {
      for (int j=0; j<2; ++j) {
        J[i].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f));
      }
    }

    for (int i=0; i<2; ++i) {
      *(J[i][i]) = *(fe_->massMat());
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
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
    fidx_ = fe_->getBoundaryDofs();
    // Construct boundary FE
    int sideset = 0;
    int numLocSides = bdryCellNodes[sideset].size();
    feBdry_.resize(numLocSides);
    for (int j = 0; j < numLocSides; ++j) {
      if (bdryCellNodes[sideset][j] != ROL::nullPtr) {
        feBdry_[j] = ROL::makePtr<FE<Real>>(bdryCellNodes[sideset][j],basisPtr_,bdryCub_,j);
      }
    }
  }

  void setFieldPattern(const std::vector<std::vector<int>> &fieldPattern) {
    fieldPattern_ = fieldPattern;
    fieldHelper_ = ROL::makePtr<FieldHelper<Real>>(numFields_, numDofs_, numFieldDofs_, fieldPattern_);
  }

  const ROL::Ptr<FE<Real>> getFE(void) const {
    return fe_;
  }

  const std::vector<ROL::Ptr<FE<Real>>> getBdryFE(void) const {
    return feBdry_;
  }

  const std::vector<std::vector<int>> getBdryCellLocIds(const int sideset = 0) const {
    return bdryCellLocIds_[sideset];
  }

  const ROL::Ptr<FieldHelper<Real>> getFieldHelper(void) const {
    return fieldHelper_;
  }

}; // PDE_GinzburgLandau


#endif
