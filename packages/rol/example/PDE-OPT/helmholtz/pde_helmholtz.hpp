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

#ifndef PDE_HELMHOLTZ_HPP
#define PDE_HELMHOLTZ_HPP

#include "../TOOLS/pde.hpp"
#include "../TOOLS/fe.hpp"
#include "../TOOLS/fieldhelper.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "Teuchos_RCP.hpp"


template <class Real>
class PDE_Helmholtz : public PDE<Real> {
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

  Real waveNumber_;
  Real RoiWaveNumber_;
  Real innerAnnulusRadius_;
  Real outerAnnulusRadius_;
  Real RoiRadius_;

  Teuchos::RCP<Intrepid::FieldContainer<Real> > ctrlWeight_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > ctrlJac_;

  Teuchos::RCP<FieldHelper<Real> > fieldHelper_;

  void computeRefractiveIndex(const Teuchos::RCP<Intrepid::FieldContainer<Real> > &kappa) const {
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
   
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (*fe_->cubPts())(i,j,k);
        }
        (*kappa)(i,j) = std::pow(evaluateRefractiveIndex(x), 2);
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
        (*F)(i,j) = evaluateForce(x,component);
      }
    }
  }

  void computeControlWeight(void) {
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
   
    ctrlWeight_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));

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

  void buildControlJacobian(void) {
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);

    // Build force/control term
    Teuchos::RCP<Intrepid::FieldContainer<Real> > F
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, p));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*F, *ctrlWeight_, *fe_->N());
    ctrlJac_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    Intrepid::FunctionSpaceTools::integrate<Real>(*ctrlJac_,
                                                  *F,
                                                  *fe_->NdetJ(),
                                                  Intrepid::COMP_CPP,
                                                  false);
  }

public:
  PDE_Helmholtz(Teuchos::ParameterList &parlist) {
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

    waveNumber_         = parlist.sublist("Problem").get("Wave Number",10.0);
    RoiWaveNumber_      = parlist.sublist("Problem").get("ROI Wave Number",10.0);

    Real dist2annulus   = parlist.sublist("Problem").get("Distance to Control Annulus",0.5);
    Real annulusWidth   = parlist.sublist("Problem").get("Control Annulus Width",0.1);
    RoiRadius_          = parlist.sublist("Problem").get("ROI Radius",2.0);
    innerAnnulusRadius_ = RoiRadius_ + dist2annulus;
    outerAnnulusRadius_ = innerAnnulusRadius_ + annulusWidth;
  }

  virtual Real evaluateRefractiveIndex(const std::vector<Real> &x) const {
    Real xnorm(0), val(0);
    const int d = x.size();
    for (int i = 0; i < d; ++i) {
      xnorm += x[i]*x[i];
    }
    xnorm = std::sqrt(xnorm);
    val   = (xnorm <= RoiRadius_) ? RoiWaveNumber_ : waveNumber_;
    return val;
  }

  virtual bool insideControlDomain(const std::vector<Real> &x) const {
    Real xnorm(0);
    const int d = x.size();
    for (int i = 0; i < d; ++i) {
      xnorm += x[i]*x[i];
    }
    xnorm = std::sqrt(xnorm);
    return (xnorm <= outerAnnulusRadius_ && xnorm >= innerAnnulusRadius_);
  }

  virtual Real evaluateForce(const std::vector<Real> &x, const int component) const {
    Real val(0);
    return val;
  }

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
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > valZ_eval(2);
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > valU_eval(2);
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > gradU_eval(2);
    for (int i=0; i<2; ++i) {
      valZ_eval[i]  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      fe_->evaluateValue(valZ_eval[i], Z[i]);
      valU_eval[i]  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      fe_->evaluateValue(valU_eval[i], U[i]);
      gradU_eval[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
      fe_->evaluateGradient(gradU_eval[i], U[i]);
    }

    // Build force/control term
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > F(2);
    for (int i=0; i<2; ++i) {
      F[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      computeForce(F[i],i);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > wZ
        = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*wZ, *ctrlWeight_, *valZ_eval[i]);
      Intrepid::RealSpaceTools<Real>::add(*F[i],*wZ);
    }

    // Build wave number
    Teuchos::RCP<Intrepid::FieldContainer<Real> > kappa
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    computeRefractiveIndex(kappa);
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > kappaU(2);
    for (int i=0; i<2; ++i) {
      kappaU[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*kappaU[i],*kappa,*valU_eval[i]);
    }

    /*******************************************************************/
    /*** Evaluate weak form of the residual.****************************/
    /*** a(ur,vr) - a(ui,vi) + b(ur,vr) + b(ui,vi) = zr(vr) - zi(vi) ***/
    /*** a(ur,vi) + a(ui,vr) - b(ur,vi) + b(ui,vr) = zr(vi) + zi(vr) ***/
    /*******************************************************************/
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > stiff(2);
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > mass(2);
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > load(2);
    for (int i=0; i<2; ++i) {
      stiff[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
      mass[i]  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
      load[i]  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));

      Intrepid::FunctionSpaceTools::integrate<Real>(*stiff[i],
                                                    *gradU_eval[i],     // grad U
                                                    *fe_->gradNdetJ(),  // grad N
                                                    Intrepid::COMP_CPP,
                                                    false);
      Intrepid::FunctionSpaceTools::integrate<Real>(*mass[i],
                                                    *kappaU[i],         // -kappa2 U
                                                    *fe_->NdetJ(),  // N
                                                    Intrepid::COMP_CPP,
                                                    false);
      Intrepid::FunctionSpaceTools::integrate<Real>(*load[i],
                                                    *F[i],              // F
                                                    *fe_->NdetJ(),      // N
                                                    Intrepid::COMP_CPP,
                                                    false);
    }

    *R[0] = *stiff[0];
    Intrepid::RealSpaceTools<Real>::subtract(*R[0],*stiff[1]);
    Intrepid::RealSpaceTools<Real>::subtract(*R[0],*mass[0]);
    Intrepid::RealSpaceTools<Real>::add(*R[0],*mass[1]);
    Intrepid::RealSpaceTools<Real>::subtract(*R[0],*load[0]);
    Intrepid::RealSpaceTools<Real>::add(*R[0],*load[1]);

    *R[1] = *stiff[0];
    Intrepid::RealSpaceTools<Real>::add(*R[1],*stiff[1]);
    Intrepid::RealSpaceTools<Real>::subtract(*R[1],*mass[0]);
    Intrepid::RealSpaceTools<Real>::subtract(*R[1],*mass[1]);
    Intrepid::RealSpaceTools<Real>::subtract(*R[1],*load[0]);
    Intrepid::RealSpaceTools<Real>::subtract(*R[1],*load[1]);

    // APPLY ROBIN CONTROLS: Sideset 0
    int sideset = 0;
    int numLocalSideIds = bdryCellLocIds_[sideset].size();
    const int numCubPerSide = bdryCub_->getNumPoints();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[sideset][j].size();
      if (numCellsSide) {
        std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > robinRes(2);
        for (int i = 0; i < 2; ++i) {
          robinRes[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, f));
          // Get U coefficients on Robin boundary
          Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff_bdry
            = getBoundaryCoeff(*U[i], sideset, j);
          // Evaluate U on FE basis
          Teuchos::RCP<Intrepid::FieldContainer<Real> > valU_eval_bdry
            = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
          feBdry_[j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
          // Compute Neumann residual
          Intrepid::FunctionSpaceTools::integrate<Real>(*robinRes[i],
                                                        *valU_eval_bdry,
                                                        *(feBdry_[j]->NdetJ()),
                                                        Intrepid::COMP_CPP, false);
        }
        // Add Neumann control residual to volume residual
        for (int k = 0; k < numCellsSide; ++k) {
          int cidx = bdryCellLocIds_[sideset][j][k];
          for (int l = 0; l < f; ++l) { 
            (*R[0])(cidx,l) += waveNumber_ * ((*robinRes[0])(k,l) + (*robinRes[1])(k,l));
            (*R[1])(cidx,l) += waveNumber_ * ((*robinRes[1])(k,l) - (*robinRes[0])(k,l));
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
 
    // Initialize Jacobians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(2);
    for (int i=0; i<2; ++i) {
      for (int j=0; j<2; ++j) {
        J[i].push_back(Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,f)));
      }
    }

    // Build wave number
    Teuchos::RCP<Intrepid::FieldContainer<Real> > kappa
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    computeRefractiveIndex(kappa);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > kappaN
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, p));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*kappaN,*kappa,*fe_->N());

    /*** Evaluate weak form of the Jacobian. ***/
    Teuchos::RCP<Intrepid::FieldContainer<Real> > kappaM
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, f));
    Intrepid::FunctionSpaceTools::integrate<Real>(*kappaM,
                                                  *kappaN,
                                                  *fe_->NdetJ(),
                                                  Intrepid::COMP_CPP,
                                                  false);

    *J[0][0] = *fe_->stiffMat();
    Intrepid::RealSpaceTools<Real>::subtract(*J[0][0],*kappaM);

    *J[0][1] = *fe_->stiffMat();
    Intrepid::RealSpaceTools<Real>::subtract(*J[0][1],*kappaM);
    Intrepid::RealSpaceTools<Real>::scale(*J[0][1],static_cast<Real>(-1));
    
    *J[1][0] = *fe_->stiffMat();
    Intrepid::RealSpaceTools<Real>::subtract(*J[1][0],*kappaM);

    *J[1][1] = *fe_->stiffMat();
    Intrepid::RealSpaceTools<Real>::subtract(*J[1][1],*kappaM);

    // APPLY ROBIN CONTROL: Sideset 0
    int sideset = 0;
    int numLocalSideIds = bdryCellLocIds_[sideset].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[sideset][j].size();
      if (numCellsSide) {
        // Add Neumann control Jacobian to volume residual
        for (int k = 0; k < numCellsSide; ++k) {
          int cidx = bdryCellLocIds_[sideset][j][k];
          for (int l = 0; l < f; ++l) { 
            for (int m = 0; m < f; ++m) { 
              (*J[0][0])(cidx,l,m) += waveNumber_ * (*feBdry_[j]->massMat())(k,l,m);
              (*J[0][1])(cidx,l,m) += waveNumber_ * (*feBdry_[j]->massMat())(k,l,m);
              (*J[1][0])(cidx,l,m) -= waveNumber_ * (*feBdry_[j]->massMat())(k,l,m);
              (*J[1][1])(cidx,l,m) += waveNumber_ * (*feBdry_[j]->massMat())(k,l,m);
            }
          }
        }
      }
    }

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

    /*** Evaluate weak form of the Jacobian. ***/
    *J[0][0] = *ctrlJac_;
    Intrepid::RealSpaceTools<Real>::scale(*J[0][0],static_cast<Real>(-1));

    *J[0][1] = *ctrlJac_;
    
    *J[1][0] = *ctrlJac_;
    Intrepid::RealSpaceTools<Real>::scale(*J[1][0],static_cast<Real>(-1));

    *J[1][1] = *ctrlJac_;
    Intrepid::RealSpaceTools<Real>::scale(*J[1][1],static_cast<Real>(-1));

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }

  void Hessian_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Helmholtz::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Helmholtz::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Helmholtz::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Helmholtz::Hessian_22): Hessian is zero.");
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
    // Compute control weight
    computeControlWeight();
    buildControlJacobian();
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

}; // PDE_Helmholtz


#endif
