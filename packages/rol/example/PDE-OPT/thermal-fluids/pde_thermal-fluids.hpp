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

/*! \file  pde_thermal-fluids.hpp
    \brief Implements the local PDE interface for the Navier-Stokes control problem.
*/

#ifndef PDE_THERMALFLUIDS_HPP
#define PDE_THERMALFLUIDS_HPP

#include "../TOOLS/pde.hpp"
#include "../TOOLS/fe.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "Teuchos_RCP.hpp"

template<class Real>
class FieldHelper {
  private:
  const int numFields_, numDofs_;
  const std::vector<int> numFieldDofs_;
  const std::vector<std::vector<int> > fieldPattern_;

  public:
  FieldHelper(const int numFields, const int numDofs,
              const std::vector<int> &numFieldDofs,
              const std::vector<std::vector<int> > &fieldPattern)
    : numFields_(numFields), numDofs_(numDofs),
      numFieldDofs_(numFieldDofs), fieldPattern_(fieldPattern) {}

  void splitFieldCoeff(std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > & U,
                       const Teuchos::RCP<const Intrepid::FieldContainer<Real> >   & u_coeff) const {
    U.resize(numFields_);
    int c = u_coeff->dimension(0);
    for (int i=0; i<numFields_; ++i) {
      U[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,numFieldDofs_[i]));
      for (int j=0; j<c; ++j) {
        for (int k=0; k<numFieldDofs_[i]; ++k) {
          //U[i](j,k) = u_coeff(j,offset[i]+k);
          (*U[i])(j,k) = (*u_coeff)(j,fieldPattern_[i][k]);
        }
      }
    }
  }

  void combineFieldCoeff(Teuchos::RCP<Intrepid::FieldContainer<Real> >   & res,
                         const std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > & R) const {
    int c = R[0]->dimension(0);  // number of cells
    res = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, numDofs_));        
    for (int i=0; i<numFields_; ++i) {
      for (int j=0; j<c; ++j) {
        for (int k=0; k<numFieldDofs_[i]; ++k) {
          (*res)(j,fieldPattern_[i][k]) = (*R[i])(j,k);
        }
      }
    }
  }

  void combineFieldCoeff(Teuchos::RCP<Intrepid::FieldContainer<Real> >   & jac,
                         const std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > & J) const {
    int c = J[0][0]->dimension(0);  // number of cells
    jac = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, numDofs_, numDofs_));        
    for (int i=0; i<numFields_; ++i) {
      for (int j=0; j<numFields_; ++j) {
        for (int k=0; k<c; ++k) {
          for (int l=0; l<numFieldDofs_[i]; ++l) {
            for (int m=0; m<numFieldDofs_[j]; ++m) {
              (*jac)(k,fieldPattern_[i][l],fieldPattern_[j][m]) = (*J[i][j])(k,l,m);
            }
          }
        }
      }
    }
  }

  int numFields(void) const {
    return numFields_;
  }

};

template <class Real>
class PDE_ThermalFluids : public PDE<Real> {
private:
  // Finite element basis information
  Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtrVel_;
  Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtrPrs_;
  Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtrThr_;
  std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;
  // Cell cubature information
  Teuchos::RCP<Intrepid::Cubature<Real> > cellCub_;
  Teuchos::RCP<Intrepid::Cubature<Real> > bdryCub_;
  // Cell node information
  Teuchos::RCP<Intrepid::FieldContainer<Real> > volCellNodes_;
  std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;
  std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
  // Finite element definition
  Teuchos::RCP<FE<Real> > feVel_;
  Teuchos::RCP<FE<Real> > fePrs_;
  Teuchos::RCP<FE<Real> > feThr_;
  std::vector<std::vector<Teuchos::RCP<FE<Real> > > > feVelBdry_;
  std::vector<std::vector<Teuchos::RCP<FE<Real> > > > fePrsBdry_;
  std::vector<std::vector<Teuchos::RCP<FE<Real> > > > feThrBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int> > fvidx_;
  std::vector<std::vector<int> > fpidx_;
  std::vector<std::vector<int> > fhidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > bdryCellDofValues_;
  // Field pattern, offsets, etc.
  std::vector<std::vector<int> > fieldPattern_;  // local Field/DOF pattern; set from DOF manager 
  int numFields_;                                // number of fields (equations in the PDE)
  int numDofs_;                                  // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                      // for each field, a counting offset
  std::vector<int> numFieldDofs_;                // for each field, number of degrees of freedom
  
  // Problem parameters.
  Real Re_, Pr_, Gr_, h_;
  const Real grav_;

  Teuchos::RCP<FieldHelper<Real> > fieldHelper_;

  Real velocityDirichletFunc(const std::vector<Real> & coords, int sideset, int locSideId, int dir) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    Real val(0);
    if ((sideset==4) && (dir==0)) {
      if ( param.size() ) {
        Real zero(0), half(0.5), one(1), two(2), pi(M_PI), four(4), c1(1), c2(-0.5);
        if ( param[0] == zero ) {
          val = half*std::sin(two*pi * (c1*coords[1] + c2));
        }
        else {
          Real num = (std::exp(four*param[0]*(c1*coords[1] + c2)) - one);
          Real den = (std::exp(two*param[0])-one);
          val = half*std::sin(pi * num/den);
        }
      }
      else {
        val = static_cast<Real>(8)
             *(coords[1]-static_cast<Real>(0.5))
             *(static_cast<Real>(1)-coords[1]);
      }
    }
    else if (((sideset==1) || (sideset==2)) && (dir==0)) {
      val = (static_cast<Real>(1)-coords[1]) * coords[1];
    }
    return val;
  }

  void computeVelocityDirichlet(void) {
    // Compute Dirichlet values at DOFs.
    int d = basisPtrVel_->getBaseCellTopology().getDimension();
    int numSidesets = bdryCellLocIds_.size();
    bdryCellDofValues_.resize(numSidesets);
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellDofValues_[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
        int f = basisPtrVel_->getCardinality();
        bdryCellDofValues_[i][j] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, d));
        Teuchos::RCP<Intrepid::FieldContainer<Real> > coords =
          Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, d));
        if (c > 0) {
          feVel_->computeDofCoords(coords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<f; ++l) {
            std::vector<Real> dofpoint(d);
            //std::cout << "Sideset " << i << " LocalSide " << j << "  Cell " << k << "  Field " << l << "  Coord ";
            for (int m=0; m<d; ++m) {
              dofpoint[m] = (*coords)(k, l, m);
              //std::cout << dofpoint[m] << "  ";
            }

            for (int m=0; m<d; ++m) {
              (*bdryCellDofValues_[i][j])(k, l, m) = velocityDirichletFunc(dofpoint, i, j, m);
              //std::cout << "  " << m << "-Value " << DirichletFunc(dofpoint, i, j, m);
            }
            //std::cout << std::endl;
          }
        }
      }
    }
  }

  Real viscosityFunc(const std::vector<Real> &x) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    if ( param.size() ) {
      return static_cast<Real>(1)/Re_;
    }
    return static_cast<Real>(1)/Re_;
  }

  Real conductivityFunc(const std::vector<Real> &x) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    if ( param.size() ) {
      return static_cast<Real>(1)/(Re_*Pr_);
    }
    return static_cast<Real>(1)/(Re_*Pr_);
  }

  Real gravityFunc(const std::vector<Real> &x) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    if ( param.size() ) {
      return grav_*Gr_/(Re_*Re_);
    }
    return grav_*Gr_/(Re_*Re_);
  }

  Real thermalRobinFunc(const std::vector<Real> &x) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    if ( param.size() ) {
      return h_ * conductivityFunc(x);
    }
    return h_ * conductivityFunc(x);
  }

  void computeCoefficients(const Teuchos::RCP<Intrepid::FieldContainer<Real> > &nu,
                           const Teuchos::RCP<Intrepid::FieldContainer<Real> > &grav,
                           const Teuchos::RCP<Intrepid::FieldContainer<Real> > &kappa) const {
    int c = feVel_->gradN()->dimension(0);
    int p = feVel_->gradN()->dimension(2);
    int d = feVel_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*feVel_->cubPts())(i,j,k);
        }
        // Compute spatially dependent coefficients
        (*nu)(i,j)    = viscosityFunc(pt);
        (*grav)(i,j)  = gravityFunc(pt);
        (*kappa)(i,j) = conductivityFunc(pt);
      }
    }
  }

  void computeThermalRobin(Teuchos::RCP<Intrepid::FieldContainer<Real> > &robin,
                           const Teuchos::RCP<Intrepid::FieldContainer<Real> > &u,
                           const Teuchos::RCP<Intrepid::FieldContainer<Real> > &z,
                           const int sideset,
                           const int locSideId,
                           const int deriv = 0,
                           const int component = 1) const {
    const int c = u->dimension(0);
    const int p = u->dimension(1);
    const int d = feThr_->gradN()->dimension(3);
    std::vector<Real> x(d);
    Real h(0);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (*feThrBdry_[sideset][locSideId]->cubPts())(i,j,k);
        }
        h = thermalRobinFunc(x);
        if ( deriv == 0 ) {
          (*robin)(i,j) = h * ((*u)(i,j) - (*z)(i,j));
        }
        else if ( deriv == 1 ) {
          (*robin)(i,j) = ((component==0) ? h : -h);
        }
        else {
          (*robin)(i,j) = static_cast<Real>(0);
        }
      }
    }
  }

  Teuchos::RCP<Intrepid::FieldContainer<Real> > getThermalBoundaryCoeff(
      const Intrepid::FieldContainer<Real> & cell_coeff,
      int sideSet, int cell) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][cell];
    const int numCellsSide = bdryCellLocId.size();
    const int f = basisPtrThr_->getCardinality();
    
    Teuchos::RCP<Intrepid::FieldContainer<Real > > bdry_coeff
      = Teuchos::rcp(new Intrepid::FieldContainer<Real > (numCellsSide, f));
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        (*bdry_coeff)(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

public:
  PDE_ThermalFluids(Teuchos::ParameterList &parlist) : grav_(-1) {
    // Finite element fields.
    basisPtrVel_ = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >);
    basisPtrPrs_ = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >);
    basisPtrThr_ = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >);
    basisPtrs_.clear();
    basisPtrs_.push_back(basisPtrVel_); // Velocity X
    basisPtrs_.push_back(basisPtrVel_); // Velocity Y
    basisPtrs_.push_back(basisPtrPrs_); // Pressure
    basisPtrs_.push_back(basisPtrThr_); // Heat
    // Volume quadrature rules.
    shards::CellTopology cellType = basisPtrs_[0]->getBaseCellTopology();        // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                           // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 4);        // set cubature degree, e.g., 4
    cellCub_ = cubFactory.create(cellType, cubDegree);                           // create default cubature
    // Boundary quadrature rules.
    int d = cellType.getDimension();
    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",4); // set cubature degree, e.g., 4
    bdryCub_ = cubFactory.create(bdryCellType, bdryCubDegree);

    // Other problem parameters.
    Re_ = parlist.sublist("Problem").get("Reynolds Number",   200.0);
    Pr_ = parlist.sublist("Problem").get("Prandtl Number",    0.72);
    Gr_ = parlist.sublist("Problem").get("Grashof Number",    40000.0);
    h_  = parlist.sublist("Problem").get("Robin Coefficient", 1.0);

    numDofs_ = 0;
    numFields_ = basisPtrs_.size();
    offset_.resize(numFields_);
    numFieldDofs_.resize(numFields_);
    for (int i=0; i<numFields_; ++i) {
      if (i==0) {
        offset_[i] = 0;
      }
      else {
        offset_[i] = offset_[i-1] + basisPtrs_[i-1]->getCardinality();
      }
      numFieldDofs_[i] = basisPtrs_[i]->getCardinality();
      numDofs_ += numFieldDofs_[i];
    }
  }

  void residual(Teuchos::RCP<Intrepid::FieldContainer<Real> > & res,
                const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Retrieve dimensions.
    int c  = u_coeff->dimension(0);
    int p  = cellCub_->getNumPoints();
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize residuals.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > R(d+2);
    for (int i = 0; i < d; ++i) {
      R[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv));
    }
    R[d]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp));
    R[d+1] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh));

    // Split u_coeff into components.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > U;
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // Evaluate problem coefficients at quadrature points.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > nu
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > grav
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > kappa
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    computeCoefficients(nu,grav,kappa);

    // Evaluate/interpolate finite element fields on cells.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valVelX_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valVelY_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valPres_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valHeat_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradVelX_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradVelY_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradHeat_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    feVel_->evaluateValue(valVelX_eval, U[0]);
    feVel_->evaluateValue(valVelY_eval, U[1]);
    fePrs_->evaluateValue(valPres_eval, U[2]);
    feThr_->evaluateValue(valHeat_eval, U[3]);
    feVel_->evaluateGradient(gradVelX_eval, U[0]);
    feVel_->evaluateGradient(gradVelY_eval, U[1]);
    feThr_->evaluateGradient(gradHeat_eval, U[3]);

    // Assemble the velocity vector and its divergence.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valVel_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > divVel_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valVel_eval)(i,j,0) = (*valVelX_eval)(i,j);
        (*valVel_eval)(i,j,1) = (*valVelY_eval)(i,j);
        (*divVel_eval)(i,j)   = (*gradVelX_eval)(i,j,0) + (*gradVelY_eval)(i,j,1);
      }
    }

    /**************************************************************************/
    /*** NAVIER STOKES ********************************************************/
    /**************************************************************************/
    Teuchos::RCP<Intrepid::FieldContainer<Real> > nuGradVelX_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > nuGradVelY_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valVelDotgradVelX_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valVelDotgradVelY_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gravForce
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    // Multiply velocity gradients with viscosity.
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(*nuGradVelX_eval, *nu, *gradVelX_eval);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(*nuGradVelY_eval, *nu, *gradVelY_eval);
    // Compute nonlinear terms in the Navier-Stokes equations.
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(*valVelDotgradVelX_eval, *valVel_eval, *gradVelX_eval);
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(*valVelDotgradVelY_eval, *valVel_eval, *gradVelY_eval);
    // Negative pressure
    Intrepid::RealSpaceTools<Real>::scale(*valPres_eval,static_cast<Real>(-1));
    // Compute gravitational force.
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*gravForce, *grav, *valHeat_eval);

    /**************************************************************************/
    /*** THERMAL **************************************************************/
    /**************************************************************************/
    Teuchos::RCP<Intrepid::FieldContainer<Real> > kappaGradHeat_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > VelDotGradHeat_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    // Multiply temperature gradient with conductivity
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(*kappaGradHeat_eval, *kappa, *gradHeat_eval);
    // Dot multiply scaled velocity with temperature gradient
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(*VelDotGradHeat_eval, *valVel_eval, *gradHeat_eval);

    /*** Evaluate weak form of the residual. ***/
    // X component of velocity equation.
    Intrepid::FunctionSpaceTools::integrate<Real>((*R[0]),
                                                  *nuGradVelX_eval,         // nu gradUX
                                                  *(feVel_->gradNdetJ()),   // gradPhi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>((*R[0]),
                                                  *valVelDotgradVelX_eval,  // (U . gradUX)
                                                  *(feVel_->NdetJ()),       // Phi
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>((*R[0]),
                                                  *valPres_eval,            // p
                                                  *(feVel_->DNDdetJ(0)),    // dPhi/dx
                                                  Intrepid::COMP_CPP,
                                                  true);
    // Y component of velocity equation.
    Intrepid::FunctionSpaceTools::integrate<Real>((*R[1]),
                                                  *nuGradVelY_eval,         // nu gradUY
                                                  *(feVel_->gradNdetJ()),   // gradPhi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>((*R[1]),
                                                  *valVelDotgradVelY_eval,  // (U . gradUY)
                                                  *(feVel_->NdetJ()),       // Phi
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>((*R[1]),
                                                  *valPres_eval,            // p
                                                  *(feVel_->DNDdetJ(1)),    // dPhi/dy
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>((*R[1]),
                                                  *gravForce,               // grav T
                                                  *(feVel_->NdetJ()),       // Phi
                                                  Intrepid::COMP_CPP,
                                                  true);
    // Pressure equation.
    Intrepid::FunctionSpaceTools::integrate<Real>((*R[2]),
                                                  *divVel_eval,             // divU
                                                  *(fePrs_->NdetJ()),       // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::RealSpaceTools<Real>::scale((*R[2]),static_cast<Real>(-1));
    // Heat equation.
    Intrepid::FunctionSpaceTools::integrate<Real>((*R[3]),
                                                  *kappaGradHeat_eval,      // kappa gradT
                                                  *(feThr_->gradNdetJ()),   // gradPhi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>((*R[3]),
                                                  *VelDotGradHeat_eval,     // U . gradT
                                                  *(feThr_->NdetJ()),       // Phi
                                                  Intrepid::COMP_CPP,
                                                  true);

    // APPLY BOUNDARY CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    const int numCubPerSide = bdryCub_->getNumPoints();
    if (numSideSets > 0) {
      // ROBIN CONDITIONS
      for (int i = 0; i < numSideSets; ++i) {
        // Control boundaries
        if (i==0 || i==3) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            Teuchos::RCP<Intrepid::FieldContainer<Real> > robinRes;
            if (numCellsSide) {
              // Get U and Z coefficients on Robin boundary
              Teuchos::RCP<Intrepid::FieldContainer<Real > > u_coeff_bdry
                = getThermalBoundaryCoeff(*(U[d+1]), i, j);
              Teuchos::RCP<Intrepid::FieldContainer<Real > > z_coeff_bdry
                = getThermalBoundaryCoeff(*(Z[d+1]), i, j);
              // Evaluate U and Z on FE basis
              Teuchos::RCP<Intrepid::FieldContainer<Real > > valU_eval_bdry
                = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
              Teuchos::RCP<Intrepid::FieldContainer<Real > > valZ_eval_bdry
                = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
              feThrBdry_[i][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
              feThrBdry_[i][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
              // Compute Stefan-Boltzmann residual
              Teuchos::RCP< Intrepid::FieldContainer<Real> > robinVal
                = Teuchos::rcp( new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
              computeThermalRobin(robinVal,valU_eval_bdry,valZ_eval_bdry,i,j,0);
              // Evaluate boundary integral
              robinRes = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, fh));
              Intrepid::FunctionSpaceTools::integrate<Real>(*robinRes,
                                                            *robinVal,
                                                            *(feThrBdry_[i][j]->NdetJ()),
                                                            Intrepid::COMP_CPP, false);
            }
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              // Thermal boundary conditions.
              if (numCellsSide) {
                for (int l = 0; l < fh; ++l) {
                  (*R[d+1])(cidx,l) += (*robinRes)(k,l);
                }
              }
            }
          }
        }
      }
      // DIRICHLET CONDITIONS
      computeVelocityDirichlet();
      for (int i = 0; i < numSideSets; ++i) {
        // Control boundaries
        if (i==0 || i==3) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            Teuchos::RCP<Intrepid::FieldContainer<Real> > robinRes;
            int numVBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numVBdryDofs; ++l) {
                // Velocity boundary conditions.
                for (int m = 0; m < d; ++m) {
                  (*R[m])(cidx,fvidx_[j][l]) = (*U[m])(cidx,fvidx_[j][l]) - (*bdryCellDofValues_[i][j])(k,fvidx_[j][l],m);
                }
              }
            }
          }
        }
        // Outflow boundaries
        if (i==1 || i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numVBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numVBdryDofs; ++l) {
                // Velocity boundary conditions.
                for (int m = 0; m < d; ++m) {
                  (*R[m])(cidx,fvidx_[j][l]) = (*U[m])(cidx,fvidx_[j][l]) - (*bdryCellDofValues_[i][j])(k,fvidx_[j][l],m);
                }
              }
              // Thermal boundary conditions.
            }
          }
        }
        // Inflow boundaries
        if (i==4) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numVBdryDofs = fvidx_[j].size();
            int numHBdryDofs = fhidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numVBdryDofs; ++l) {
                // Velocity boundary conditions.
                for (int m = 0; m < d; ++m) {
                  (*R[m])(cidx,fvidx_[j][l]) = (*U[m])(cidx,fvidx_[j][l]) - (*bdryCellDofValues_[i][j])(k,fvidx_[j][l],m);
                }
              }
              // Thermal boundary conditions.
              for (int l = 0; l < numHBdryDofs; ++l) {
                (*R[d+1])(cidx,fhidx_[j][l]) = (*U[d+1])(cidx,fhidx_[j][l]);
              }
            }
          }
        }
        // Step boundaries
        if (i==5 || i==6 || i==7 || i==8) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numVBdryDofs = (i==7) ? 1 : fvidx_[j].size();
            int numHBdryDofs = (i==7) ? 1 : fhidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numVBdryDofs; ++l) {
                // Velocity boundary conditions.
                for (int m = 0; m < d; ++m) {
                  (*R[m])(cidx,fvidx_[j][l]) = (*U[m])(cidx,fvidx_[j][l]) - (*bdryCellDofValues_[i][j])(k,fvidx_[j][l],m);
                }
              }
              // Thermal boundary conditions.
              for (int l = 0; l < numHBdryDofs; ++l) {
                (*R[d+1])(cidx,fhidx_[j][l]) = (*U[d+1])(cidx,fhidx_[j][l]) - static_cast<Real>(1);
              }
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
    int c  = u_coeff->dimension(0);
    int p  = cellCub_->getNumPoints();
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(d+2);
    for (int i = 0; i < d+2; ++i) {
      J[i].resize(d+2);
    }
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        J[i][j] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fv));
      }
    }
    for (int i = 0; i < d; ++i) {
      J[d][i]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fv));
      J[i][d]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fp));
      J[d+1][i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, fv));
      J[i][d+1] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fh));
    }
    J[d][d]     = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fp));
    J[d+1][d]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, fp));
    J[d][d+1]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fh));
    J[d+1][d+1] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, fh));

    // Split u_coeff into components.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > U;
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // Evaluate problem coefficients at quadrature points.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > nu
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > grav
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > kappa
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    computeCoefficients(nu,grav,kappa);

    // Evaluate/interpolate finite element fields on cells.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valVelX_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valVelY_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradVelX_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradVelY_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradHeat_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    feVel_->evaluateValue(valVelX_eval, U[0]);
    feVel_->evaluateValue(valVelY_eval, U[1]);
    feVel_->evaluateGradient(gradVelX_eval, U[0]);
    feVel_->evaluateGradient(gradVelY_eval, U[1]);
    feThr_->evaluateGradient(gradHeat_eval, U[3]);

    // Assemble the velocity vector and its divergence.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valVel_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > ddxVelX_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > ddyVelX_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > ddxVelY_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > ddyVelY_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradHeatX_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradHeatY_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valVel_eval)(i,j,0)  = (*valVelX_eval)(i,j);
        (*valVel_eval)(i,j,1)  = (*valVelY_eval)(i,j);
        (*ddxVelX_eval)(i,j)   = (*gradVelX_eval)(i,j,0);
        (*ddyVelX_eval)(i,j)   = (*gradVelX_eval)(i,j,1);
        (*ddxVelY_eval)(i,j)   = (*gradVelY_eval)(i,j,0);
        (*ddyVelY_eval)(i,j)   = (*gradVelY_eval)(i,j,1);
        (*gradHeatX_eval)(i,j) = (*gradHeat_eval)(i,j,0);
        (*gradHeatY_eval)(i,j) = (*gradHeat_eval)(i,j,1);
      }
    }

    /**************************************************************************/
    /*** NAVIER STOKES ********************************************************/
    /**************************************************************************/
    Teuchos::RCP<Intrepid::FieldContainer<Real> > nuGradPhiX_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, p, d));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > nuGradPhiY_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, p, d));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valVelDotgradPhiX_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valVelDotgradPhiY_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > ddyVelYPhiY_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > ddxVelYPhiX_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > ddxVelXPhiX_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > ddyVelXPhiY_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gravPhi
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, p));
    Intrepid::FieldContainer<Real> negPrsPhi(*(fePrs_->N()));
    // Multiply velocity gradients with viscosity.
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(*nuGradPhiX_eval, *nu, *(feVel_->gradN()));
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(*nuGradPhiY_eval, *nu, *(feVel_->gradN()));
    // Compute nonlinear terms in the Navier-Stokes equations.
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*valVelDotgradPhiX_eval, *valVel_eval, *(feVel_->gradN()));
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*valVelDotgradPhiY_eval, *valVel_eval, *(feVel_->gradN()));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*ddxVelXPhiX_eval, *ddxVelX_eval, *(feVel_->N()));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*ddyVelXPhiY_eval, *ddyVelX_eval, *(feVel_->N()));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*ddxVelYPhiX_eval, *ddxVelY_eval, *(feVel_->N()));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*ddyVelYPhiY_eval, *ddyVelY_eval, *(feVel_->N()));
    // Negative pressure basis.
    Intrepid::RealSpaceTools<Real>::scale(negPrsPhi,static_cast<Real>(-1));
    // Compute gravity Jacobian
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*gravPhi, *grav, *(feThr_->N()));

    /**************************************************************************/
    /*** THERMAL **************************************************************/
    /**************************************************************************/
    Teuchos::RCP<Intrepid::FieldContainer<Real> > kappaGradPhi
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, p, d));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > ddxHeatPhi
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > ddyHeatPhi
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > VelPhi
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, p));
    // Compute kappa times gradient of basis.
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(*kappaGradPhi, *kappa, *(feThr_->gradN()));
    // Thermal-velocity Jacobians.
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*ddxHeatPhi, *gradHeatX_eval, *feThr_->N());
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*ddyHeatPhi, *gradHeatY_eval, *feThr_->N());
    // Compute scaled velocity.
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*VelPhi, *valVel_eval, *feThr_->gradN());

    /*** Evaluate weak form of the Jacobian. ***/
    // X component of velocity equation.
    Intrepid::FunctionSpaceTools::integrate<Real>((*J[0][0]),
                                                  *nuGradPhiX_eval,         // nu gradPhi
                                                  *(feVel_->gradNdetJ()),   // gradPhi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>((*J[0][0]),
                                                  *(feVel_->NdetJ()),       // Phi
                                                  *valVelDotgradPhiX_eval,  // (U . gradPhiX)
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>((*J[0][0]),
                                                  *(feVel_->NdetJ()),       // Phi
                                                  *ddxVelXPhiX_eval,        // (Phi . gradU)
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>((*J[0][1]),
                                                  *(feVel_->NdetJ()),       // Phi
                                                  *ddyVelXPhiY_eval,        // (Phi . gradU)
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>((*J[0][2]),
                                                  *(feVel_->DNDdetJ(0)),    // dPhi/dx
                                                  negPrsPhi,                // -Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    // Y component of velocity equation.
    Intrepid::FunctionSpaceTools::integrate<Real>((*J[1][0]),
                                                  *(feVel_->NdetJ()),       // Phi
                                                  *ddxVelYPhiX_eval,        // (Phi . gradU)
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>((*J[1][1]),
                                                  *nuGradPhiY_eval,         // nu gradPhi
                                                  *(feVel_->gradNdetJ()),   // gradPhi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>((*J[1][1]),
                                                  *(feVel_->NdetJ()),       // Phi
                                                  *valVelDotgradPhiY_eval,  // (U . gradPhiX)
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>((*J[1][1]),
                                                  *(feVel_->NdetJ()),       // Phi
                                                  *ddyVelYPhiY_eval,        // (Phi . gradU)
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>((*J[1][2]),
                                                  *(feVel_->DNDdetJ(1)),    // dPhi/dx
                                                  negPrsPhi,                // -Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>((*J[1][3]),
                                                  *(feVel_->NdetJ()),       // Phi
                                                  *gravPhi,                 // grav Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    // Pressure equation.
    Intrepid::FunctionSpaceTools::integrate<Real>((*J[2][0]),
                                                  *(fePrs_->NdetJ()),       // Phi
                                                  *(feVel_->DND(0)),        // dPhi/dx
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>((*J[2][1]),
                                                  *(fePrs_->NdetJ()),       // Phi
                                                  *(feVel_->DND(1)),        // dPhi/dx
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::RealSpaceTools<Real>::scale((*J[2][0]),static_cast<Real>(-1));
    Intrepid::RealSpaceTools<Real>::scale((*J[2][1]),static_cast<Real>(-1));
    // Heat equation.
    Intrepid::FunctionSpaceTools::integrate<Real>((*J[3][0]),
                                                  *ddxHeatPhi,              // gradT Phi
                                                  *(feVel_->NdetJ()),       // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>((*J[3][1]),
                                                  *ddyHeatPhi,              // gradT Phi
                                                  *(feVel_->NdetJ()),       // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>((*J[3][3]),
                                                  *kappaGradPhi,            // kappa gradPhi
                                                  *(feThr_->gradNdetJ()),   // gradPhi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>((*J[3][3]),
                                                  *(feThr_->NdetJ()),       // Phi
                                                  *VelPhi,                  // U . gradPhi
                                                  Intrepid::COMP_CPP,
                                                  true);

    // APPLY DIRICHLET CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    const int numCubPerSide = bdryCub_->getNumPoints();
    if (numSideSets > 0) {
      // ROBIN CONDITIONS
      for (int i = 0; i < numSideSets; ++i) {
        // Control boundaries
        if (i==0 || i==3) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            Teuchos::RCP<Intrepid::FieldContainer<Real> > robinJac;
            if (numCellsSide) {
              // Get U and Z coefficients on Robin boundary
              Teuchos::RCP<Intrepid::FieldContainer<Real > > u_coeff_bdry
                = getThermalBoundaryCoeff(*(U[d+1]), i, j);
              Teuchos::RCP<Intrepid::FieldContainer<Real > > z_coeff_bdry
                = getThermalBoundaryCoeff(*(Z[d+1]), i, j);
              // Evaluate U and Z on FE basis
              Teuchos::RCP<Intrepid::FieldContainer<Real > > valU_eval_bdry
                = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
              Teuchos::RCP<Intrepid::FieldContainer<Real > > valZ_eval_bdry
                = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
              feThrBdry_[i][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
              feThrBdry_[i][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
              // Compute Stefan-Boltzmann residual
              Teuchos::RCP< Intrepid::FieldContainer<Real> > robinVal1
                = Teuchos::rcp( new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
              Intrepid::FieldContainer<Real> robinVal(numCellsSide, fh, numCubPerSide);
              computeThermalRobin(robinVal1,valU_eval_bdry,valZ_eval_bdry,i,j,1,0);
              Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(robinVal,
                                                                          *robinVal1,
                                                                          *(feThrBdry_[i][j]->N()));
              // Evaluate boundary integral
              robinJac = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, fh, fh));
              Intrepid::FunctionSpaceTools::integrate<Real>(*robinJac,
                                                            robinVal,
                                                            *(feThrBdry_[i][j]->NdetJ()),
                                                            Intrepid::COMP_CPP, false);
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < fh; ++l) {
                  for (int m = 0; m < fh; ++m) {
                    (*J[d+1][d+1])(cidx,l,m) += (*robinJac)(k,l,m);
                  }
                }
              }
            }
          }
        }
      }
      // DIRICHLET CONDITIONS
      for (int i = 0; i < numSideSets; ++i) {
        // Control boundaries
        if (i==0 || i==3) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            Teuchos::RCP<Intrepid::FieldContainer<Real> > robinJac;
            int numVBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numVBdryDofs; ++l) {
                for (int m = 0; m < fv; ++m) {
                  for (int n = 0; n < d; ++n) {
                    for (int p = 0; p < d; ++p) {
                      (*J[n][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                    }
                    (*J[n][n])(cidx,fvidx_[j][l],fvidx_[j][l]) = static_cast<Real>(1);
                  }
                }
                for (int m = 0; m < fp; ++m) {
                  for (int n = 0; n < d; ++n) {
                    (*J[n][d])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
                for (int m = 0; m < fh; ++m) {
                  for (int n = 0; n < d; ++n) {
                    (*J[n][d+1])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
              }
            }
          }
        }
        // Outflow boundaries
        if (i==1 || i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m = 0; m < fv; ++m) {
                  for (int n = 0; n < d; ++n) {
                    for (int p = 0; p < d; ++p) {
                      (*J[n][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                    }
                    (*J[n][n])(cidx,fvidx_[j][l],fvidx_[j][l]) = static_cast<Real>(1);
                  }
                }
                for (int m = 0; m < fp; ++m) {
                  for (int n=0; n < d; ++n) {
                    (*J[n][d])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
                for (int m = 0; m < fh; ++m) {
                  for (int n=0; n < d; ++n) {
                    (*J[n][d+1])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
              }
            }
          }
        }
        // Inflow boundaries
        if (i==4) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numVBdryDofs = fvidx_[j].size();
            int numHBdryDofs = fhidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numVBdryDofs; ++l) {
                for (int m = 0; m < fv; ++m) {
                  for (int n = 0; n < d; ++n) {
                    for (int p = 0; p < d; ++p) {
                      (*J[n][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                    }
                    (*J[n][n])(cidx,fvidx_[j][l],fvidx_[j][l]) = static_cast<Real>(1);
                  }
                }
                for (int m = 0; m < fp; ++m) {
                  for (int n = 0; n < d; ++n) {
                    (*J[n][d])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
                for (int m = 0; m < fh; ++m) {
                  for (int n = 0; n < d; ++n) {
                    (*J[n][d+1])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
              }
              for (int l = 0; l < numHBdryDofs; ++l) {
                for (int m = 0; m < fv; ++m) {
                  for (int n = 0; n < d; ++n) {
                    (*J[d+1][n])(cidx,fhidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
                for (int m = 0; m < fp; ++m) {
                  (*J[d+1][d])(cidx,fhidx_[j][l],m) = static_cast<Real>(0);
                }
                for (int m = 0; m < fh; ++m) {
                  (*J[d+1][d+1])(cidx,fhidx_[j][l],m) = static_cast<Real>(0);
                }
                (*J[d+1][d+1])(cidx,fhidx_[j][l],fhidx_[j][l]) = static_cast<Real>(1);
              }
            }
          }
        }
        // Step boundaries
        if (i==5 || i==6 || i==7 || i==8) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numVBdryDofs = (i==7) ? 1 : fvidx_[j].size();
            int numHBdryDofs = (i==7) ? 1 : fhidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numVBdryDofs; ++l) {
                for (int m = 0; m < fv; ++m) {
                  for (int n = 0; n < d; ++n) {
                    for (int p = 0; p < d; ++p) {
                      (*J[n][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                    }
                    (*J[n][n])(cidx,fvidx_[j][l],fvidx_[j][l]) = static_cast<Real>(1);
                  }
                }
                for (int m = 0; m < fp; ++m) {
                  for (int n = 0; n < d; ++n) {
                    (*J[n][d])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
                for (int m = 0; m < fh; ++m) {
                  for (int n = 0; n < d; ++n) {
                    (*J[n][d+1])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
              }
              for (int l = 0; l < numHBdryDofs; ++l) {
                for (int m = 0; m < fv; ++m) {
                  for (int n = 0; n < d; ++n) {
                    (*J[d+1][n])(cidx,fhidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
                for (int m = 0; m < fp; ++m) {
                  (*J[d+1][d])(cidx,fhidx_[j][l],m) = static_cast<Real>(0);
                }
                for (int m = 0; m < fh; ++m) {
                  (*J[d+1][d+1])(cidx,fhidx_[j][l],m) = static_cast<Real>(0);
                }
                (*J[d+1][d+1])(cidx,fhidx_[j][l],fhidx_[j][l]) = static_cast<Real>(1);
              }
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
    int c  = u_coeff->dimension(0);
    //int p  = cellCub_->getNumPoints();
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(d+2);
    for (int i = 0; i < d+2; ++i) {
      J[i].resize(d+2);
    }
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        J[i][j] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fv));
      }
    }
    for (int i = 0; i < d; ++i) {
      J[d][i]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fv));
      J[i][d]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fp));
      J[d+1][i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, fv));
      J[i][d+1] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fh));
    }
    J[d][d]     = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fp));
    J[d+1][d]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, fp));
    J[d][d+1]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fh));
    J[d+1][d+1] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, fh));

    // Split u_coeff into components.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > U;
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // APPLY DIRICHLET CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    const int numCubPerSide = bdryCub_->getNumPoints();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        // Control boundaries
        if (i==0 || i==3) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            if (numCellsSide) {
              // Get U and Z coefficients on Robin boundary
              Teuchos::RCP<Intrepid::FieldContainer<Real > > u_coeff_bdry
                = getThermalBoundaryCoeff(*(U[d+1]), i, j);
              Teuchos::RCP<Intrepid::FieldContainer<Real > > z_coeff_bdry
                = getThermalBoundaryCoeff(*(Z[d+1]), i, j);
              // Evaluate U and Z on FE basis
              Teuchos::RCP<Intrepid::FieldContainer<Real > > valU_eval_bdry
                = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
              Teuchos::RCP<Intrepid::FieldContainer<Real > > valZ_eval_bdry
                = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
              feThrBdry_[i][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
              feThrBdry_[i][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
              // Compute Stefan-Boltzmann residual
              Teuchos::RCP< Intrepid::FieldContainer<Real> > robinVal1
                = Teuchos::rcp( new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
              Intrepid::FieldContainer<Real> robinVal(numCellsSide, fh, numCubPerSide);
              computeThermalRobin(robinVal1,valU_eval_bdry,valZ_eval_bdry,i,j,1,1);
              Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(robinVal,
                                                                          *robinVal1,
                                                                          *(feThrBdry_[i][j]->N()));
              // Evaluate boundary integral
              Intrepid::FieldContainer<Real> robinJac(numCellsSide, fh, fh);
              Intrepid::FunctionSpaceTools::integrate<Real>(robinJac,
                                                            robinVal,
                                                            *(feThrBdry_[i][j]->NdetJ()),
                                                            Intrepid::COMP_CPP, false);
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < fh; ++l) {
                  for (int m = 0; m < fh; ++m) {
                    (*J[d+1][d+1])(cidx,l,m) += robinJac(k,l,m);
                  }
                }
              }
            }
          }
        }
      }
      // DIRICHLET CONDITIONS
      for (int i = 0; i < numSideSets; ++i) {
        // Control boundaries
        if (i==0 || i==3) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            Teuchos::RCP<Intrepid::FieldContainer<Real> > robinJac;
            int numVBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numVBdryDofs; ++l) {
                for (int m = 0; m < fv; ++m) {
                  for (int n = 0; n < d; ++n) {
                    for (int p = 0; p < d; ++p) {
                      (*J[n][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                    }
                  }
                }
                for (int m = 0; m < fp; ++m) {
                  for (int n = 0; n < d; ++n) {
                    (*J[n][d])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
                for (int m = 0; m < fh; ++m) {
                  for (int n = 0; n < d; ++n) {
                    (*J[n][d+1])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
              }
            }
          }
        }
        // Outflow boundaries
        if (i==1 || i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m = 0; m < fv; ++m) {
                  for (int n = 0; n < d; ++n) {
                    for (int p = 0; p < d; ++p) {
                      (*J[n][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                    }
                  }
                }
                for (int m = 0; m < fp; ++m) {
                  for (int n=0; n < d; ++n) {
                    (*J[n][d])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
                for (int m = 0; m < fh; ++m) {
                  for (int n=0; n < d; ++n) {
                    (*J[n][d+1])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
              }
            }
          }
        }
        // Inflow boundaries
        if (i==4) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numVBdryDofs = fvidx_[j].size();
            int numHBdryDofs = fhidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numVBdryDofs; ++l) {
                for (int m = 0; m < fv; ++m) {
                  for (int n = 0; n < d; ++n) {
                    for (int p = 0; p < d; ++p) {
                      (*J[n][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                    }
                  }
                }
                for (int m = 0; m < fp; ++m) {
                  for (int n = 0; n < d; ++n) {
                    (*J[n][d])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
                for (int m = 0; m < fh; ++m) {
                  for (int n = 0; n < d; ++n) {
                    (*J[n][d+1])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
              }
              for (int l = 0; l < numHBdryDofs; ++l) {
                for (int m = 0; m < fv; ++m) {
                  for (int n = 0; n < d; ++n) {
                    (*J[d+1][n])(cidx,fhidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
                for (int m = 0; m < fp; ++m) {
                  (*J[d+1][d])(cidx,fhidx_[j][l],m) = static_cast<Real>(0);
                }
                for (int m = 0; m < fh; ++m) {
                  (*J[d+1][d+1])(cidx,fhidx_[j][l],m) = static_cast<Real>(0);
                }
              }
            }
          }
        }
        // Step boundaries
        if (i==5 || i==6 || i==7 || i==8) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numVBdryDofs = (i==7) ? 1 : fvidx_[j].size();
            int numHBdryDofs = (i==7) ? 1 : fhidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numVBdryDofs; ++l) {
                for (int m = 0; m < fv; ++m) {
                  for (int n = 0; n < d; ++n) {
                    for (int p = 0; p < d; ++p) {
                      (*J[n][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                    }
                  }
                }
                for (int m = 0; m < fp; ++m) {
                  for (int n = 0; n < d; ++n) {
                    (*J[n][d])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
                for (int m = 0; m < fh; ++m) {
                  for (int n = 0; n < d; ++n) {
                    (*J[n][d+1])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
              }
              for (int l = 0; l < numHBdryDofs; ++l) {
                for (int m = 0; m < fv; ++m) {
                  for (int n = 0; n < d; ++n) {
                    (*J[d+1][n])(cidx,fhidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
                for (int m = 0; m < fp; ++m) {
                  (*J[d+1][d])(cidx,fhidx_[j][l],m) = static_cast<Real>(0);
                }
                for (int m = 0; m < fh; ++m) {
                  (*J[d+1][d+1])(cidx,fhidx_[j][l],m) = static_cast<Real>(0);
                }
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
//    throw Exception::NotImplemented("");
    // Retrieve dimensions.
    int c  = u_coeff->dimension(0);
    int p  = cellCub_->getNumPoints();
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > H(d+2);
    for (int i = 0; i < d+2; ++i) {
      H[i].resize(d+2);
    }
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        H[i][j] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fv));
      }
    }
    for (int i = 0; i < d; ++i) {
      H[d][i]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fv));
      H[i][d]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fp));
      H[d+1][i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, fv));
      H[i][d+1] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fh));
    }
    H[d][d]     = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fp));
    H[d+1][d]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, fp));
    H[d][d+1]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fh));
    H[d+1][d+1] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, fh));

    // Split l_coeff into components.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > L;
    fieldHelper_->splitFieldCoeff(L, l_coeff);

    // APPLY DIRICHLET CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      // DIRICHLET CONDITIONS
      for (int i = 0; i < numSideSets; ++i) {
        // Control boundaries
        if (i==0 || i==3) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numVBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numVBdryDofs; ++l) {
                for (int m = 0; m < d; ++m) {
                  (*L[m])(cidx,fvidx_[j][l]) = static_cast<Real>(0);
                }
              }
            }
          }
        }
        // Outflow boundaries
        if (i==1 || i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numVBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numVBdryDofs; ++l) {
                for (int m = 0; m < d; ++m) {
                  (*L[m])(cidx,fvidx_[j][l]) = static_cast<Real>(0);
                }
              }
            }
          }
        }
        // Inflow boundaries
        if (i==4) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numVBdryDofs = fvidx_[j].size();
            int numHBdryDofs = fhidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numVBdryDofs; ++l) {
                for (int m = 0; m < d; ++m) {
                  (*L[m])(cidx,fvidx_[j][l]) = static_cast<Real>(0);
                }
              }
              for (int l = 0; l < numHBdryDofs; ++l) {
                (*L[d+1])(cidx,fhidx_[j][l]) = static_cast<Real>(0);
              }
            }
          }
        }
        // Step boundaries
        if (i==5 || i==6 || i==7 || i==8) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numVBdryDofs = (i==7) ? 1 : fvidx_[j].size();
            int numHBdryDofs = (i==7) ? 1 : fhidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numVBdryDofs; ++l) {
                for (int m = 0; m < d; ++m) {
                  (*L[m])(cidx,fvidx_[j][l]) = static_cast<Real>(0);
                }
              }
              for (int l = 0; l < numHBdryDofs; ++l) {
                (*L[d+1])(cidx,fhidx_[j][l]) = static_cast<Real>(0);
              }
            }
          }
        }
      }
    }

    // Evaluate/interpolate finite element fields on cells.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valVelX_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valVelY_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valHeat_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    feVel_->evaluateValue(valVelX_eval, L[0]);
    feVel_->evaluateValue(valVelY_eval, L[1]);
    feThr_->evaluateValue(valHeat_eval, L[3]);

    // Compute nonlinear terms in the Navier-Stokes equations.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valVelXPhi_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valVelYPhi_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, p));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*valVelXPhi_eval, *valVelX_eval, *(feVel_->N()));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*valVelYPhi_eval, *valVelY_eval, *(feVel_->N()));

    // Compute nonlinear terms in the thermal equation.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valHeatVelPhi =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, p));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*valHeatVelPhi, *valHeat_eval, *(feVel_->N()));

    /*** Evaluate weak form of the Hessian. ***/
    // X component of velocity equation.
    Intrepid::FunctionSpaceTools::integrate<Real>((*H[0][0]),
                                                  *valVelXPhi_eval,        // L Phi
                                                  *(feVel_->DNDdetJ(0)),   // dPhi/dx
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>((*H[0][0]),
                                                  *(feVel_->DNDdetJ(0)),   // dPhi/dx
                                                  *valVelXPhi_eval,        // L Phi
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>((*H[0][1]),
                                                  *(feVel_->DNDdetJ(1)),   // dPhi/dy
                                                  *valVelXPhi_eval,        // L Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>((*H[0][1]),
                                                  *valVelYPhi_eval,        // L Phi
                                                  *(feVel_->DNDdetJ(0)),   // dPhi/dx
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>((*H[0][3]),
                                                  *valHeatVelPhi,          // L Phi
                                                  *(feThr_->DNDdetJ(0)),   // dPhi/dx
                                                  Intrepid::COMP_CPP,
                                                  false);
    // Y component of velocity equation.
    Intrepid::FunctionSpaceTools::integrate<Real>((*H[1][0]),
                                                  *(feVel_->DNDdetJ(0)),   // dPhi/dx
                                                  *valVelYPhi_eval,        // L Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>((*H[1][0]),
                                                  *valVelXPhi_eval,        // L Phi
                                                  *(feVel_->DNDdetJ(1)),   // dPhi/dy
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>((*H[1][1]),
                                                  *valVelYPhi_eval,        // L Phi
                                                  *(feVel_->DNDdetJ(1)),   // dPhi/dy
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>((*H[1][1]),
                                                  *(feVel_->DNDdetJ(1)),   // dPhi/dy
                                                  *valVelYPhi_eval,        // L Phi
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>((*H[1][3]),
                                                  *valHeatVelPhi,          // L Phi
                                                  *(feThr_->DNDdetJ(1)),   // dPhi/dy
                                                  Intrepid::COMP_CPP,
                                                  false);
    // Thermal equation.
    Intrepid::FunctionSpaceTools::integrate<Real>((*H[3][0]),
                                                  *(feThr_->DNDdetJ(0)),   // dPhi/dx
                                                  *valHeatVelPhi,          // L Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>((*H[3][1]),
                                                  *(feThr_->DNDdetJ(1)),   // dPhi/dy
                                                  *valHeatVelPhi,          // L Phi
                                                  Intrepid::COMP_CPP,
                                                  false);

    // Combine the Hessians.
    fieldHelper_->combineFieldCoeff(hess, H);

  }

  void Hessian_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & riesz) {
    // Retrieve dimensions.
    int c  = feVel_->N()->dimension(0);
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(d+2);
    for (int i = 0; i < d+2; ++i) {
      J[i].resize(d+2);
    }
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        J[i][j] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fv));
      }
    }
    for (int i = 0; i < d; ++i) {
      J[d][i]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fv));
      J[i][d]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fp));
      J[d+1][i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, fv));
      J[i][d+1] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fh));
    }
    J[d][d]     = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fp));
    J[d+1][d]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, fp));
    J[d][d+1]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fh));
    J[d+1][d+1] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, fh));

    for (int i = 0; i < d; ++i) {
      *(J[i][i]) = *(feVel_->stiffMat());
      Intrepid::RealSpaceTools<Real>::add(*(J[i][i]),*(feVel_->massMat()));
    }
    *(J[d][d]) = *(fePrs_->massMat());
    *(J[d+1][d+1]) = *(feThr_->stiffMat());
    Intrepid::RealSpaceTools<Real>::add(*(J[d+1][d+1]),*(feThr_->massMat()));

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  void RieszMap_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & riesz) {
    // Retrieve dimensions.
    int c  = feVel_->N()->dimension(0);
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(d+2);
    for (int i = 0; i < d+2; ++i) {
      J[i].resize(d+2);
    }
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        J[i][j] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fv));
      }
    }
    for (int i = 0; i < d; ++i) {
      J[d][i]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fv));
      J[i][d]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fp));
      J[d+1][i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, fv));
      J[i][d+1] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fh));
    }
    J[d][d]     = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fp));
    J[d+1][d]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, fp));
    J[d][d+1]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fh));
    J[d+1][d+1] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, fh));


    for (int i = 0; i < d; ++i) {
      *(J[i][i]) = *(feVel_->massMat());
    }
    *(J[d][d]) = *(fePrs_->massMat());
    *(J[d+1][d+1]) = *(feThr_->massMat());

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
    feVel_ = Teuchos::rcp(new FE<Real>(volCellNodes_,basisPtrVel_,cellCub_));
    fePrs_ = Teuchos::rcp(new FE<Real>(volCellNodes_,basisPtrPrs_,cellCub_));
    feThr_ = Teuchos::rcp(new FE<Real>(volCellNodes_,basisPtrThr_,cellCub_));
    // Get boundary degrees of freedom.
    fvidx_ = feVel_->getBoundaryDofs();
    fpidx_ = fePrs_->getBoundaryDofs();
    fhidx_ = feThr_->getBoundaryDofs();
    // Construct boundary FEs
    const int numSideSets = bdryCellNodes.size();
    feVelBdry_.resize(numSideSets);
    fePrsBdry_.resize(numSideSets);
    feThrBdry_.resize(numSideSets);
    for (int i = 0; i < numSideSets; ++i) {
      int numLocSides = bdryCellNodes[i].size();
      feVelBdry_[i].resize(numLocSides);
      fePrsBdry_[i].resize(numLocSides);
      feThrBdry_[i].resize(numLocSides);
      for (int j = 0; j < numLocSides; ++j) {
        if (bdryCellNodes[i][j] != Teuchos::null) {
          feVelBdry_[i][j] = Teuchos::rcp(new FE<Real>(bdryCellNodes[i][j],basisPtrVel_,bdryCub_,j));
          fePrsBdry_[i][j] = Teuchos::rcp(new FE<Real>(bdryCellNodes[i][j],basisPtrPrs_,bdryCub_,j));
          feThrBdry_[i][j] = Teuchos::rcp(new FE<Real>(bdryCellNodes[i][j],basisPtrThr_,bdryCub_,j));
        }
      }
    }
  }

  void setFieldPattern(const std::vector<std::vector<int> > & fieldPattern) {
    fieldPattern_ = fieldPattern;
    fieldHelper_ = Teuchos::rcp(new FieldHelper<Real>(numFields_, numDofs_, numFieldDofs_, fieldPattern_));
  }

  const Teuchos::RCP<FE<Real> > getVelocityFE(void) const {
    return feVel_;
  }

  const Teuchos::RCP<FE<Real> > getPressureFE(void) const {
    return fePrs_;
  }

  const Teuchos::RCP<FE<Real> > getThermalFE(void) const {
    return feThr_;
  }

  const std::vector<std::vector<Teuchos::RCP<FE<Real> > > > getVelocityBdryFE(void) const {
    return feVelBdry_;
  }

  const std::vector<std::vector<Teuchos::RCP<FE<Real> > > > getPressureBdryFE(void) const {
    return fePrsBdry_;
  }

  const std::vector<std::vector<Teuchos::RCP<FE<Real> > > > getThermalBdryFE(void) const {
    return feThrBdry_;
  }

  const std::vector<std::vector<std::vector<int> > > getBdryCellLocIds(void) const {
    return bdryCellLocIds_;
  }

  const Teuchos::RCP<FieldHelper<Real> > getFieldHelper(void) const {
    return fieldHelper_;
  }

}; // PDE_ThermalFluids

#endif
