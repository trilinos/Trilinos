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
#include "../TOOLS/fieldhelper.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "Teuchos_RCP.hpp"

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
  std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > bdryCellVDofValues_;
  std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > bdryCellTDofValues_;
  // Field pattern, offsets, etc.
  std::vector<std::vector<int> > fieldPattern_;  // local Field/DOF pattern; set from DOF manager 
  int numFields_;                                // number of fields (equations in the PDE)
  int numDofs_;                                  // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                      // for each field, a counting offset
  std::vector<int> numFieldDofs_;                // for each field, number of degrees of freedom
  
  // Problem parameters.
  Real Re_, Pr_, Gr_, h_, Tinflow_, Tstep_;
  const Real grav_;

  Teuchos::RCP<FieldHelper<Real> > fieldHelper_;

  Real velocityDirichletFunc(const std::vector<Real> & coords, int sideset, int locSideId, int dir) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    Real val(0);
    if ((sideset==4) && (dir==0)) {
      if ( param.size() >= 85 ) {
        Real zero(0), one(1), two(2), pi(M_PI), four(4), c1(1), c2(-0.5);
        if ( param[4] == zero ) {
          val = std::sin(two*pi * (c1*coords[1] + c2));
        }
        else {
          Real num = (std::exp(four*param[4]*(c1*coords[1] + c2)) - one);
          Real den = (std::exp(two*param[4])-one);
          val = std::sin(pi * num/den);
        }
      }
      else {
        val = static_cast<Real>(8)
             *(coords[1]-static_cast<Real>(0.5))
             *(static_cast<Real>(1)-coords[1]);
      }
    }
    else if (((sideset==1) || (sideset==2)) && (dir==0)) {
      if ( param.size() >= 86 ) {
        Real zero(0), half(0.5), one(1), two(2), pi(M_PI), four(4), c1(0.5), c2(0);
        if ( param[5] == zero ) {
          val = half*half*std::sin(two*pi * (c1*coords[1] + c2));
        }
        else {
          Real num = (std::exp(four*param[5]*(c1*coords[1] + c2)) - one);
          Real den = (std::exp(two*param[5])-one);
          val = half*half*std::sin(pi * num/den);
        }
      }
      else {
        val = (static_cast<Real>(1)-coords[1]) * coords[1];
      }
    }
    return val;
  }

  Real thermalDirichletFunc(const std::vector<Real> & coords, int sideset, int locSideId) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    Real val(0);
    if (sideset==4) {
      if ( param.size() >= 87 ) {
        Real l(-0.01), u(0.01);
        val = (u-l)*param[6] + l;
      }
      else {
        val = Tinflow_;
      }
    }
    else if ((sideset==5) || (sideset==6) || (sideset==7) || (sideset==8)) {
      if ( param.size() >= 88 ) {
        Real l(0.99), u(1.01);
        val = (u-l)*param[7] + l;
      }
      else {
        val = Tstep_;
      }
    }
    return val;
  }

  void computeDirichlet(void) {
    // Compute Dirichlet values at DOFs.
    int d  = basisPtrVel_->getBaseCellTopology().getDimension();
    int fv = basisPtrVel_->getCardinality();
    int ft = basisPtrThr_->getCardinality();
    int numSidesets = bdryCellLocIds_.size();
    bdryCellVDofValues_.resize(numSidesets);
    bdryCellTDofValues_.resize(numSidesets);
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellVDofValues_[i].resize(numLocSides);
      bdryCellTDofValues_[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
        bdryCellVDofValues_[i][j] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, d));
        bdryCellTDofValues_[i][j] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, ft));
        Teuchos::RCP<Intrepid::FieldContainer<Real> > Vcoords, Tcoords;
        Vcoords = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, d));
        Tcoords = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, ft, d));
        if (c > 0) {
          feVel_->computeDofCoords(Vcoords, bdryCellNodes_[i][j]);
          feThr_->computeDofCoords(Tcoords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<fv; ++l) {
            std::vector<Real> dofpoint(d);
            //std::cout << "Sideset " << i << " LocalSide " << j << "  Cell " << k << "  Field " << l << "  Coord ";
            for (int m=0; m<d; ++m) {
              dofpoint[m] = (*Vcoords)(k, l, m);
              //std::cout << dofpoint[m] << "  ";
            }

            for (int m=0; m<d; ++m) {
              (*bdryCellVDofValues_[i][j])(k, l, m) = velocityDirichletFunc(dofpoint, i, j, m);
              //std::cout << "  " << m << "-Value " << DirichletFunc(dofpoint, i, j, m);
            }
            //std::cout << std::endl;
          }
          for (int l=0; l<ft; ++l) {
            std::vector<Real> dofpoint(d);
            //std::cout << "Sideset " << i << " LocalSide " << j << "  Cell " << k << "  Field " << l << "  Coord ";
            for (int m=0; m<d; ++m) {
              dofpoint[m] = (*Tcoords)(k, l, m);
              //std::cout << dofpoint[m] << "  ";
            }

            (*bdryCellTDofValues_[i][j])(k, l) = thermalDirichletFunc(dofpoint, i, j);
            //std::cout << "    Value " << DirichletFunc(dofpoint, i, j);
            //std::cout << std::endl;
          }
        }
      }
    }
  }

  Real ReynoldsNumber(void) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    if ( param.size() >= 1 ) {
      Real l(190), u(210);
      return (u-l)*param[0] + l;
    }
    return Re_;
  }

  Real PrandtlNumber(void) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    if ( param.size() >= 2 ) {
      Real l(0.62), u(0.82);
      return (u-l)*param[1] + l;
    }
    return Pr_;
  }

  Real GrashofNumber(void) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    if ( param.size() >= 3 ) {
      Real l(38000), u(50000);
      return (u-l)*param[2] + l;
    }
    return Gr_;
  }
  
  Real ThicknessNumber(void) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    if ( param.size() >= 4 ) {
      Real l(1), u(15);
      return (u-l)*param[3] + l;
    }
    return h_;
  }

  Real viscosityFunc(const std::vector<Real> &x) const {
    Real ReynoldsNum = ReynoldsNumber();
    return static_cast<Real>(1)/ReynoldsNum;
  }

  Real conductivityFunc(const std::vector<Real> &x) const {
    Real ReynoldsNum = ReynoldsNumber();
    Real PrandtlNum = PrandtlNumber();
    return static_cast<Real>(1)/(ReynoldsNum*PrandtlNum);
  }

  Real gravityFunc(const std::vector<Real> &x) const {
    Real ReynoldsNum = ReynoldsNumber();
    Real GrashofNum = GrashofNumber();
    return grav_*GrashofNum/(ReynoldsNum*ReynoldsNum);
  }

  Real thermalRobinFunc(const std::vector<Real> &x) const {
    Real H = ThicknessNumber();
    return H * conductivityFunc(x);
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
    Real h(0), kappa(0);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (*feThrBdry_[sideset][locSideId]->cubPts())(i,j,k);
        }
        h = thermalRobinFunc(x);
        kappa = conductivityFunc(x);
        if ( deriv == 0 ) {
          (*robin)(i,j) = kappa * h * ((*u)(i,j) - (*z)(i,j));
        }
        else if ( deriv == 1 ) {
          (*robin)(i,j) = ((component==0) ? kappa * h : -kappa * h);
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
    // Finite element fields -- NOT DIMENSION INDEPENDENT!
    basisPtrVel_ = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >);
    basisPtrPrs_ = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >);
    basisPtrThr_ = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >);
    // Volume quadrature rules.
    shards::CellTopology cellType = basisPtrVel_->getBaseCellTopology();         // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                           // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 4);        // set cubature degree, e.g., 4
    cellCub_ = cubFactory.create(cellType, cubDegree);                           // create default cubature
    // Boundary quadrature rules.
    int d = cellType.getDimension();
    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",4); // set cubature degree, e.g., 4
    bdryCub_ = cubFactory.create(bdryCellType, bdryCubDegree);
    // Fill finite element basis container
    basisPtrs_.clear();
    for (int i = 0; i < d; ++i) {
      basisPtrs_.push_back(basisPtrVel_); // Velocity
    }
    basisPtrs_.push_back(basisPtrPrs_); // Pressure
    basisPtrs_.push_back(basisPtrThr_); // Heat

    // Other problem parameters.
    Re_      = parlist.sublist("Problem").get("Reynolds Number",    200.0);
    Pr_      = parlist.sublist("Problem").get("Prandtl Number",     0.72);
    Gr_      = parlist.sublist("Problem").get("Grashof Number",     40000.0);
    h_       = parlist.sublist("Problem").get("Robin Coefficient",  1.0);
    Tinflow_ = parlist.sublist("Problem").get("Inflow Temperature", 0.0);
    Tstep_   = parlist.sublist("Problem").get("Step Temperature",   1.0);

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
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > U, Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // Evaluate problem coefficients at quadrature points.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > nu, grav, kappa;
    nu    = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    grav  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    kappa = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    computeCoefficients(nu,grav,kappa);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > valVel_vec(d);
    for (int i = 0; i < d; ++i) {
      valVel_vec[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      feVel_->evaluateValue(valVel_vec[i], U[i]);
    }
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valPres_eval, valHeat_eval;
    valPres_eval = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    valHeat_eval = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fePrs_->evaluateValue(valPres_eval, U[d]);
    feThr_->evaluateValue(valHeat_eval, U[d+1]);
    // Evaluate/interpolate gradient of finite element fields on cells.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > gradVel_vec(d);
    for (int i = 0; i < d; ++i) {
      gradVel_vec[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
      feVel_->evaluateGradient(gradVel_vec[i], U[i]);
    }
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradHeat_eval;
    gradHeat_eval = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    feThr_->evaluateGradient(gradHeat_eval, U[d+1]);

    // Assemble the velocity vector and its divergence.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valVel_eval, divVel_eval;
    valVel_eval = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    divVel_eval = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*divVel_eval)(i,j) = static_cast<Real>(0);
        for (int k = 0; k < d; ++k) {
          (*valVel_eval)(i,j,k) = (*valVel_vec[k])(i,j);
          (*divVel_eval)(i,j)  += (*gradVel_vec[k])(i,j,k);
        }
      }
    }

    /**************************************************************************/
    /*** NAVIER STOKES ********************************************************/
    /**************************************************************************/
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > nuGradVel_vec(d), valVelDotgradVel_vec(d);
    for (int i = 0; i < d; ++i) {
      // Multiply velocity gradients with viscosity.
      nuGradVel_vec[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
      Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(*nuGradVel_vec[i], *nu, *gradVel_vec[i]);
      // Compute nonlinear terms in the Navier-Stokes equations.
      valVelDotgradVel_vec[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(*valVelDotgradVel_vec[i], *valVel_eval, *gradVel_vec[i]);
    }
    // Negative pressure
    Intrepid::RealSpaceTools<Real>::scale(*valPres_eval,static_cast<Real>(-1));
    // Compute gravitational force.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gravForce;
    gravForce = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*gravForce, *grav, *valHeat_eval);

    /**************************************************************************/
    /*** THERMAL **************************************************************/
    /**************************************************************************/
    Teuchos::RCP<Intrepid::FieldContainer<Real> > kappaGradHeat_eval, velDotGradHeat_eval;
    kappaGradHeat_eval  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    velDotGradHeat_eval = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    // Multiply temperature gradient with conductivity
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(*kappaGradHeat_eval, *kappa, *gradHeat_eval);
    // Dot multiply scaled velocity with temperature gradient
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(*velDotGradHeat_eval, *valVel_eval, *gradHeat_eval);

    /**************************************************************************/
    /*** EVALUATE WEAK FORM OF RESIDUAL ***************************************/
    /**************************************************************************/
    // Velocity equation.
    for (int i = 0; i < d; ++i) {
      Intrepid::FunctionSpaceTools::integrate<Real>((*R[i]),
                                                    *nuGradVel_vec[i],        // nu gradUX
                                                    *(feVel_->gradNdetJ()),   // gradPhi
                                                    Intrepid::COMP_CPP,
                                                    false);
      Intrepid::FunctionSpaceTools::integrate<Real>((*R[i]),
                                                    *valVelDotgradVel_vec[i], // (U . gradUX)
                                                    *(feVel_->NdetJ()),       // Phi
                                                    Intrepid::COMP_CPP,
                                                    true);
      Intrepid::FunctionSpaceTools::integrate<Real>((*R[i]),
                                                    *valPres_eval,            // p
                                                    *(feVel_->DNDdetJ(i)),    // dPhi/dx
                                                    Intrepid::COMP_CPP,
                                                    true);
    }
    // Apply gravitational force.
    Intrepid::FunctionSpaceTools::integrate<Real>((*R[d-1]),
                                                  *gravForce,               // grav T
                                                  *(feVel_->NdetJ()),       // Phi
                                                  Intrepid::COMP_CPP,
                                                  true);
    // Pressure equation.
    Intrepid::FunctionSpaceTools::integrate<Real>((*R[d]),
                                                  *divVel_eval,             // divU
                                                  *(fePrs_->NdetJ()),       // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::RealSpaceTools<Real>::scale((*R[d]),static_cast<Real>(-1));
    // Heat equation.
    Intrepid::FunctionSpaceTools::integrate<Real>((*R[d+1]),
                                                  *kappaGradHeat_eval,      // kappa gradT
                                                  *(feThr_->gradNdetJ()),   // gradPhi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>((*R[d+1]),
                                                  *velDotGradHeat_eval,     // U . gradT
                                                  *(feThr_->NdetJ()),       // Phi
                                                  Intrepid::COMP_CPP,
                                                  true);

    /**************************************************************************/
    /*** APPLY BOUNDARY CONDITIONS ********************************************/
    /**************************************************************************/
    // --> Control boundaries: i=0,3
    // --> Outflow boundaries: i=1,2
    // -->  Inflow boundaries: i=4
    // -->    Step boundaries: i=5,6,7,8
    // -->   Pressure Pinning: i=9
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
            if (numCellsSide) {
              // Get U and Z coefficients on Robin boundary
              Teuchos::RCP<Intrepid::FieldContainer<Real > > u_coeff_bdry, z_coeff_bdry;
              u_coeff_bdry = getThermalBoundaryCoeff(*(U[d+1]), i, j);
              z_coeff_bdry = getThermalBoundaryCoeff(*(Z[d+1]), i, j);
              // Evaluate U and Z on FE basis
              Teuchos::RCP<Intrepid::FieldContainer<Real > > valU_eval_bdry, valZ_eval_bdry;
              valU_eval_bdry = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
              valZ_eval_bdry = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
              feThrBdry_[i][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
              feThrBdry_[i][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
              // Compute Stefan-Boltzmann residual
              Teuchos::RCP<Intrepid::FieldContainer<Real> > robinVal;
              robinVal = Teuchos::rcp( new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
              computeThermalRobin(robinVal,valU_eval_bdry,valZ_eval_bdry,i,j,0);
              // Evaluate boundary integral
              Teuchos::RCP<Intrepid::FieldContainer<Real> > robinRes;
              robinRes = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, fh));
              Intrepid::FunctionSpaceTools::integrate<Real>(*robinRes,
                                                            *robinVal,
                                                            *(feThrBdry_[i][j]->NdetJ()),
                                                            Intrepid::COMP_CPP, false);
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                // Thermal boundary conditions.
                for (int l = 0; l < fh; ++l) {
                  (*R[d+1])(cidx,l) += (*robinRes)(k,l);
                }
              }
            }
          }
        }
      }
      // DIRICHLET CONDITIONS
      computeDirichlet();
      // Velocity Boundary Conditions
      for (int i = 0; i < numSideSets; ++i) {
        if ( (i==0 || i==3) || (i==1 || i==2) || (i==4) || (i==5 || i==6 || i==8) ) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numVBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numVBdryDofs; ++l) {
                for (int m = 0; m < d; ++m) {
                  (*R[m])(cidx,fvidx_[j][l]) = (*U[m])(cidx,fvidx_[j][l]) - (*bdryCellVDofValues_[i][j])(k,fvidx_[j][l],m);
                }
              }
            }
          }
        }
        if (i==7) {
          int j = 0, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) { 
            int cidx = bdryCellLocIds_[i][0][0];
            for (int m=0; m < d; ++m) {
              (*R[m])(cidx,fvidx_[j][l]) = (*U[m])(cidx,fvidx_[j][l]) - (*bdryCellVDofValues_[i][j])(0,fvidx_[j][l],m);
            }
          }
        }
        // Thermal Boundary Conditions
        if ( (i==4) || (i==5 || i==6 || i==8) ) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numHBdryDofs = fhidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numHBdryDofs; ++l) {
                (*R[d+1])(cidx,fhidx_[j][l]) = (*U[d+1])(cidx,fhidx_[j][l]) - (*bdryCellTDofValues_[i][j])(k,fhidx_[j][l]);
              }
            }
          }
        }
        if (i==7) {
          int j = 0, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) { 
            int cidx = bdryCellLocIds_[i][0][0];
            (*R[d+1])(cidx,fhidx_[j][l]) = (*U[d+1])(cidx,fhidx_[j][l]) - (*bdryCellTDofValues_[i][j])(0,fhidx_[j][l]);
          }
        }
        // Pressure pinning
        if (i==9) {
          int j = 2, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) { 
            int cidx = bdryCellLocIds_[i][0][0];
            (*R[d])(cidx,fpidx_[j][l]) = (*U[d])(cidx,fpidx_[j][l]);
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
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > U, Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // Evaluate problem coefficients at quadrature points.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > nu, grav, kappa;
    nu    = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    grav  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    kappa = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    computeCoefficients(nu,grav,kappa);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > valVel_vec(d);
    for (int i = 0; i < d; ++i) {
      valVel_vec[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      feVel_->evaluateValue(valVel_vec[i], U[i]);
    }
    // Evaluate/interpolate gradient of finite element fields on cells.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > gradVel_vec(d);
    for (int i = 0; i < d; ++i) {
      gradVel_vec[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
      feVel_->evaluateGradient(gradVel_vec[i], U[i]);
    }
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradHeat_eval;
    gradHeat_eval = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    feThr_->evaluateGradient(gradHeat_eval, U[d+1]);

    // Assemble the velocity vector and its divergence.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valVel_eval;
    valVel_eval = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > dVel_vec(d);
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > gradHeat_vec(d);
    for (int k = 0; k < d; ++k) {
      dVel_vec[k].resize(d);
      for (int l = 0; l < d; ++l) {
        dVel_vec[k][l] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
        for (int i = 0; i < c; ++i) {
          for (int j = 0; j < p; ++j) {
            (*dVel_vec[k][l])(i,j) = (*gradVel_vec[k])(i,j,l);
          }
        }
      }
      gradHeat_vec[k] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          (*valVel_eval)(i,j,k)  = (*valVel_vec[k])(i,j);
          (*gradHeat_vec[k])(i,j) = (*gradHeat_eval)(i,j,k);
        }
      }
    }

    /**************************************************************************/
    /*** NAVIER STOKES ********************************************************/
    /**************************************************************************/
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > dVelPhi_vec(d);
    for (int i = 0; i < d; ++i) {
      // Compute nonlinear terms in the Navier-Stokes equations.
      dVelPhi_vec[i].resize(d);
      for (int j = 0; j < d; ++j) {
        dVelPhi_vec[i][j] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, p));
        Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*dVelPhi_vec[i][j], *dVel_vec[i][j], *(feVel_->N()));
      }
    }
    // Multiply velocity gradients with viscosity.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > nuGradPhi_eval;
    nuGradPhi_eval = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, p, d));
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(*nuGradPhi_eval, *nu, *(feVel_->gradN()));
    // Compute nonlinear terms in the Navier-Stokes equations.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valVelDotgradPhi_eval;
    valVelDotgradPhi_eval = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, p));
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*valVelDotgradPhi_eval, *valVel_eval, *(feVel_->gradN()));
    // Negative pressure basis.
    Intrepid::FieldContainer<Real> negPrsPhi(*(fePrs_->N()));
    Intrepid::RealSpaceTools<Real>::scale(negPrsPhi,static_cast<Real>(-1));
    // Compute gravity Jacobian
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gravPhi;
    gravPhi = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, p));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*gravPhi, *grav, *(feThr_->N()));

    /**************************************************************************/
    /*** THERMAL **************************************************************/
    /**************************************************************************/
    Teuchos::RCP<Intrepid::FieldContainer<Real> > kappaGradPhi, VelPhi;
    kappaGradPhi = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, p, d));
    VelPhi = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, p));
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > dHeatPhi(d);
    // Compute kappa times gradient of basis.
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(*kappaGradPhi, *kappa, *(feThr_->gradN()));
    // Compute scaled velocity.
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*VelPhi, *valVel_eval, *feThr_->gradN());
    // Thermal-velocity Jacobians.
    for (int i = 0; i < d; ++i) {
      dHeatPhi[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fh, p));
      Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*dHeatPhi[i], *gradHeat_vec[i], *feThr_->N());
    }

    /*** Evaluate weak form of the Jacobian. ***/
    // X component of velocity equation.
    for (int i = 0; i < d; ++i) {
      // Velocity components
      for (int j = 0; j < d; ++j) {
        Intrepid::FunctionSpaceTools::integrate<Real>((*J[i][j]),
                                                      *(feVel_->NdetJ()),     // Phi
                                                      *dVelPhi_vec[i][j],     // (Phi . gradU)
                                                      Intrepid::COMP_CPP,
                                                      false);
      }
      Intrepid::FunctionSpaceTools::integrate<Real>((*J[i][i]),
                                                    *nuGradPhi_eval,          // nu gradPhi
                                                    *(feVel_->gradNdetJ()),   // gradPhi
                                                    Intrepid::COMP_CPP,
                                                    true);
      Intrepid::FunctionSpaceTools::integrate<Real>((*J[i][i]),
                                                    *(feVel_->NdetJ()),       // Phi
                                                    *valVelDotgradPhi_eval,   // (U . gradPhiX)
                                                    Intrepid::COMP_CPP,
                                                    true);
      // Pressure components
      Intrepid::FunctionSpaceTools::integrate<Real>((*J[i][d]),
                                                    *(feVel_->DNDdetJ(i)),    // dPhi/dx
                                                    negPrsPhi,                // -Phi
                                                    Intrepid::COMP_CPP,
                                                    false);
      Intrepid::FunctionSpaceTools::integrate<Real>((*J[d][i]),
                                                    *(fePrs_->NdetJ()),       // Phi
                                                    *(feVel_->DND(i)),        // dPhi/dx
                                                    Intrepid::COMP_CPP,
                                                    false);
      Intrepid::RealSpaceTools<Real>::scale((*J[d][i]),static_cast<Real>(-1));
      // Thermal components
      Intrepid::FunctionSpaceTools::integrate<Real>((*J[d+1][i]),
                                                    *dHeatPhi[i],             // gradT Phi
                                                    *(feVel_->NdetJ()),       // Phi
                                                    Intrepid::COMP_CPP,
                                                    false);
    }
    // Gravitational component
    Intrepid::FunctionSpaceTools::integrate<Real>((*J[d-1][d+1]),
                                                  *(feVel_->NdetJ()),       // Phi
                                                  *gravPhi,                 // grav Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    // Thermal components
    Intrepid::FunctionSpaceTools::integrate<Real>((*J[d+1][d+1]),
                                                  *kappaGradPhi,            // kappa gradPhi
                                                  *(feThr_->gradNdetJ()),   // gradPhi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>((*J[d+1][d+1]),
                                                  *(feThr_->NdetJ()),       // Phi
                                                  *VelPhi,                  // U . gradPhi
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
            Teuchos::RCP<Intrepid::FieldContainer<Real> > robinJac;
            if (numCellsSide) {
              // Get U and Z coefficients on Robin boundary
              Teuchos::RCP<Intrepid::FieldContainer<Real > > u_coeff_bdry, z_coeff_bdry;
              u_coeff_bdry = getThermalBoundaryCoeff(*(U[d+1]), i, j);
              z_coeff_bdry = getThermalBoundaryCoeff(*(Z[d+1]), i, j);
              // Evaluate U and Z on FE basis
              Teuchos::RCP<Intrepid::FieldContainer<Real > > valU_eval_bdry, valZ_eval_bdry;
              valU_eval_bdry = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
              valZ_eval_bdry = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
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
        // Velocity Boundary Conditions
        if ( (i==0 || i==3) || (i==1 || i==2) || (i==4) || (i==5 || i==6 || i==8) ) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
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
                for (int m = 0; m < d; ++m) {
                  for (int n = 0; n < fp; ++n) {
                    (*J[m][d])(cidx,fvidx_[j][l],n) = static_cast<Real>(0);
                  }
                  for (int n = 0; n < fh; ++n) {
                    (*J[m][d+1])(cidx,fvidx_[j][l],n) = static_cast<Real>(0);
                  }
                }
              }
            }
          }
        }
        if (i==7) {
          int j = 0, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) { 
            int cidx = bdryCellLocIds_[i][0][0];
            for (int m=0; m < fv; ++m) {
              for (int n=0; n < d; ++n) {
                for (int p=0; p < d; ++p) {
                  (*J[n][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                }
                (*J[n][n])(cidx,fvidx_[j][l],fvidx_[j][l]) = static_cast<Real>(1);
              }
            }
            for (int m=0; m < d; ++m) {
              for (int n=0; n < fp; ++n) {
                (*J[m][d])(cidx,fvidx_[j][l],n) = static_cast<Real>(0);
              }
              for (int n = 0; n < fh; ++n) {
                (*J[m][d+1])(cidx,fvidx_[j][l],n) = static_cast<Real>(0);
              }
            }
          }
        }
        // Thermal boundary conditions
        if ( (i==4) || (i==5 || i==6 || i==8) ) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numHBdryDofs = (i==7) ? 1 : fhidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
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
        if (i==7) {
          int j = 0, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) { 
            int cidx = bdryCellLocIds_[i][0][0];
            for (int m=0; m < fv; ++m) {
              for (int n=0; n < d; ++n) {
                for (int p=0; p < d; ++p) {
                  (*J[d+1][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                }
              }
            }
            for (int m=0; m < fp; ++m) {
              (*J[d+1][d])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
            }
            for (int m = 0; m < fh; ++m) {
              (*J[d+1][d+1])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
            }
            (*J[d+1][d+1])(cidx,fhidx_[j][l],fhidx_[j][l]) = static_cast<Real>(1);
          }
        }
        // Pressure pinning
        if (i==9) {
          int j = 2, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) {
            int cidx = bdryCellLocIds_[i][0][0];
            for (int m = 0; m < fv; ++m) {
              for (int n = 0; n < d; ++n) {
                (*J[d][n])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
              }
            }
            for (int m = 0; m < fp; ++m) {
              (*J[d][d])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
            }
            (*J[d][d])(cidx,fpidx_[j][l],fpidx_[j][l]) = static_cast<Real>(1);
            for (int m = 0; m < fh; ++m) {
              (*J[d][d+1])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
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
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > U, Z;
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
        // Velocity Boundary Conditions
        if ( (i==0 || i==3) || (i==1 || i==2) || (i==4) || (i==5 || i==6 || i==8) ) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
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
                for (int m = 0; m < d; ++m) {
                  for (int n = 0; n < fp; ++n) {
                    (*J[m][d])(cidx,fvidx_[j][l],n) = static_cast<Real>(0);
                  }
                  for (int n = 0; n < fh; ++n) {
                    (*J[m][d+1])(cidx,fvidx_[j][l],n) = static_cast<Real>(0);
                  }
                }
              }
            }
          }
        }
        if (i==7) {
          int j = 0, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) { 
            int cidx = bdryCellLocIds_[i][0][0];
            for (int m=0; m < fv; ++m) {
              for (int n=0; n < d; ++n) {
                for (int p=0; p < d; ++p) {
                  (*J[n][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                }
              }
            }
            for (int m=0; m < d; ++m) {
              for (int n=0; n < fp; ++n) {
                (*J[m][d])(cidx,fvidx_[j][l],n) = static_cast<Real>(0);
              }
              for (int n = 0; n < fh; ++n) {
                (*J[m][d+1])(cidx,fvidx_[j][l],n) = static_cast<Real>(0);
              }
            }
          }
        }
        // Thermal boundary conditions
        if ( (i==4) || (i==5 || i==6 || i==8) ) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numHBdryDofs = (i==7) ? 1 : fhidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
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
        if (i==7) {
          int j = 0, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) { 
            int cidx = bdryCellLocIds_[i][0][0];
            for (int m=0; m < fv; ++m) {
              for (int n=0; n < d; ++n) {
                for (int p=0; p < d; ++p) {
                  (*J[d+1][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                }
              }
            }
            for (int m=0; m < fp; ++m) {
              (*J[d+1][d])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
            }
            for (int m = 0; m < fh; ++m) {
              (*J[d+1][d+1])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
            }
          }
        }
        // Pressure pinning
        if (i==9) {
          int j = 2, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) {
            int cidx = bdryCellLocIds_[i][0][0];
            for (int m = 0; m < fv; ++m) {
              for (int n = 0; n < d; ++n) {
                (*J[d][n])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
              }
            }
            for (int m = 0; m < fp; ++m) {
              (*J[d][d])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
            }
            for (int m = 0; m < fh; ++m) {
              (*J[d][d+1])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
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
        // Velocity boundary conditions
        if ( (i==0 || i==3) || (i==1 || i==2) || (i==4) || (i==5 || i==6 || i==8) ) {
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
        if (i==7) {
          int j = 0, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) { 
            int cidx = bdryCellLocIds_[i][0][0];
            for (int m=0; m < d; ++m) {
              (*L[m])(cidx,fvidx_[j][l]) = static_cast<Real>(0);
            }
          }
        }
        // Inflow boundaries
        if ( (i==4) || (i==5 || i==6 || i==8) ) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numHBdryDofs = fhidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numHBdryDofs; ++l) {
                (*L[d+1])(cidx,fhidx_[j][l]) = static_cast<Real>(0);
              }
            }
          }
        }
        if (i==7) {
          int j = 0, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) { 
            int cidx = bdryCellLocIds_[i][0][0];
            (*L[d+1])(cidx,fhidx_[j][l]) = static_cast<Real>(0);
          }
        }
        // Pressure pinning
        if (i==9) {
          int j = 2, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) { 
            int cidx = bdryCellLocIds_[i][0][0];
            (*L[d])(cidx,fpidx_[j][l]) = static_cast<Real>(0);
          }
        }
      }
    }

    // Evaluate/interpolate finite element fields on cells.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > valVel_vec(d);
    for (int i = 0; i < d; ++i) {
      valVel_vec[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      feVel_->evaluateValue(valVel_vec[i], L[i]);
    }
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valHeat_eval;
    valHeat_eval = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    feThr_->evaluateValue(valHeat_eval, L[3]);

    // Compute nonlinear terms in the Navier-Stokes equations.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > valVelPhi_vec(d);
    for (int i = 0; i < d; ++i) {
      valVelPhi_vec[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, p));
      Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*valVelPhi_vec[i], *valVel_vec[i], *(feVel_->N()));
    }

    // Compute nonlinear terms in the thermal equation.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valHeatVelPhi;
    valHeatVelPhi = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, p));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*valHeatVelPhi, *valHeat_eval, *(feVel_->N()));

    /*** Evaluate weak form of the Hessian. ***/
    for (int i = 0; i < d; ++i) {
      // velocity equation.
      for (int j = 0; j < d; ++j) {
        Intrepid::FunctionSpaceTools::integrate<Real>((*H[i][j]),
                                                      *valVelPhi_vec[j],       // L Phi
                                                      *(feVel_->DNDdetJ(i)),   // dPhi/dx
                                                      Intrepid::COMP_CPP,
                                                      false);
        Intrepid::FunctionSpaceTools::integrate<Real>((*H[i][j]),
                                                      *(feVel_->DNDdetJ(j)),   // dPhi/dx
                                                      *valVelPhi_vec[i],       // L Phi
                                                      Intrepid::COMP_CPP,
                                                      true);
      }
      Intrepid::FunctionSpaceTools::integrate<Real>((*H[i][d+1]),
                                                    *valHeatVelPhi,          // L Phi
                                                    *(feThr_->DNDdetJ(i)),   // dPhi/dx
                                                    Intrepid::COMP_CPP,
                                                    false);
      // Thermal equation.
      Intrepid::FunctionSpaceTools::integrate<Real>((*H[d+1][i]),
                                                    *(feThr_->DNDdetJ(i)),   // dPhi/dx
                                                    *valHeatVelPhi,          // L Phi
                                                    Intrepid::COMP_CPP,
                                                    false);
    }

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
