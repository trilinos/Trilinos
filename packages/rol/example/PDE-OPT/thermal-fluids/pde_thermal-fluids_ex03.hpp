// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde_thermal-fluids.hpp
    \brief Implements the local PDE interface for the Navier-Stokes control problem.
*/

#ifndef PDE_THERMALFLUIDS_EX03_HPP
#define PDE_THERMALFLUIDS_EX03_HPP

#include "../TOOLS/pde.hpp"
#include "../TOOLS/fe.hpp"
#include "../TOOLS/fieldhelper.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real>
class PDE_ThermalFluids_ex03 : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtrVel_;
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtrPrs_;
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtrThr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real> > cellCub_;
  ROL::Ptr<Intrepid::Cubature<Real> > bdryCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real> > volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;
  std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real> > feVel_;
  ROL::Ptr<FE<Real> > fePrs_;
  ROL::Ptr<FE<Real> > feThr_;
  std::vector<std::vector<ROL::Ptr<FE<Real> > > > feVelBdry_;
  std::vector<std::vector<ROL::Ptr<FE<Real> > > > fePrsBdry_;
  std::vector<std::vector<ROL::Ptr<FE<Real> > > > feThrBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int> > fvidx_;
  std::vector<std::vector<int> > fpidx_;
  std::vector<std::vector<int> > fhidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellVDofValues_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellTDofValues_;
  // Field pattern, offsets, etc.
  std::vector<std::vector<int> > fieldPattern_;  // local Field/DOF pattern; set from DOF manager 
  int numFields_;                                // number of fields (equations in the PDE)
  int numDofs_;                                  // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                      // for each field, a counting offset
  std::vector<int> numFieldDofs_;                // for each field, number of degrees of freedom
  
  // Problem parameters.
  Real Re_, Pr_, Gr_, h_;
  const Real grav_;
  int Nbottom_, Nleft_, Nright_;
  Real ReScale_, PrScale_, GrScale_, hScale_, TScale_;
  bool pinPressure_;

  ROL::Ptr<FieldHelper<Real> > fieldHelper_;

  Real velocityDirichletFunc(const std::vector<Real> & coords, int sideset, int locSideId, int dir) const {
    Real val(0);
    if ((sideset==3) && (dir==1)) {
      const Real one(1), two(2), three(3), four(4), x = coords[0];
      if (x <= one/three) {
        val = two*(one/three - x)*x;
      }
      else if (x > one/three && x < two/three) {
        val = -four*(x-one/three)*(two/three-x);
      }
      else if (x >= two/three) {
        val = two*(x-two/three)*(one-x);
      }
    }
    return val;
  }

  Real thermalDirichletFunc(const std::vector<Real> & coords, int sideset, int locSideId) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    Real val(0);
    if (sideset==0) {
      val = static_cast<Real>(1);
      if (param.size()) {
        Real root2(std::sqrt(2.0)), pi(M_PI), ln2(std::log(2.0));
        for (int i = 0; i < Nbottom_; ++i) {
          Real di(i+1);
          val += TScale_ * param[i]/ln2 * (root2 * std::sin(di * pi * coords[0]))/(di * pi);
        }
      }
    }
    else if (sideset==5) {
      val = static_cast<Real>(0);
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
        bdryCellVDofValues_[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, d);
        bdryCellTDofValues_[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, ft);
        ROL::Ptr<Intrepid::FieldContainer<Real> > Vcoords, Tcoords;
        Vcoords = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, d);
        Tcoords = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, ft, d);
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
    Real val = Re_;
    if (param.size()) {
      int offset = Nbottom_ + Nright_ + Nleft_;
      Real one(1);
      val /= (one + ReScale_*param[offset]);
    }
    return val;
  }

  Real PrandtlNumber(void) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    Real val = Pr_;
    if (param.size()) {
      int offset = Nbottom_ + Nright_ + Nleft_;
      Real one(1);
      val *= (one + ReScale_*param[offset])/(one + PrScale_*param[offset+1]);
    }
    return val;
  }

  Real GrashofNumber(void) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    Real val = Gr_;
    if (param.size()) {
      int offset = Nbottom_ + Nright_ + Nleft_;
      Real one(1), muscale = one + ReScale_*param[offset];
      val *= (one + GrScale_*param[offset+2])/(muscale*muscale);
    }
    return val;
  }
  
  Real ThicknessNumber(const std::vector<Real> &x, const int sideset) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
    Real val = h_;
    if ( param.size()) {
      if ( sideset == 2 ) {
        int offset = Nbottom_;
        Real root2(std::sqrt(2.0)), pi(M_PI), ln2(std::log(2.0));
        for (int i = 0; i < Nleft_; ++i) {
          Real di(i+1);
          val += hScale_ * param[offset + i]/ln2 * (root2 * std::sin(di * pi * x[1]))/(di * pi);
        }
      }
      else if ( sideset == 1 ) {
        int offset = Nbottom_ + Nleft_;
        Real root2(std::sqrt(2.0)), pi(M_PI), ln2(std::log(2.0));
        for (int i = 0; i < Nright_; ++i) {
          Real di(i+1);
          val += hScale_ * param[offset + i]/ln2 * (root2 * std::sin(di * pi * x[1]))/(di * pi);
        }
      }
    }
    return val;
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

  Real thermalRobinFunc(const std::vector<Real> &x, const int sideset) const {
    return ThicknessNumber(x, sideset) * conductivityFunc(x);
  }

  void computeCoefficients(const ROL::Ptr<Intrepid::FieldContainer<Real> > &nu,
                           const ROL::Ptr<Intrepid::FieldContainer<Real> > &grav,
                           const ROL::Ptr<Intrepid::FieldContainer<Real> > &kappa) const {
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

  void computeThermalRobin(ROL::Ptr<Intrepid::FieldContainer<Real> > &robin,
                           const ROL::Ptr<Intrepid::FieldContainer<Real> > &u,
                           const ROL::Ptr<Intrepid::FieldContainer<Real> > &z,
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
        h = thermalRobinFunc(x, sideset);
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

  ROL::Ptr<Intrepid::FieldContainer<Real> > getThermalBoundaryCoeff(
      const Intrepid::FieldContainer<Real> & cell_coeff,
      int sideSet, int cell) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][cell];
    const int numCellsSide = bdryCellLocId.size();
    const int f = basisPtrThr_->getCardinality();
    
    ROL::Ptr<Intrepid::FieldContainer<Real > > bdry_coeff
      = ROL::makePtr<Intrepid::FieldContainer<Real >>(numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        (*bdry_coeff)(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

public:
  PDE_ThermalFluids_ex03(Teuchos::ParameterList &parlist) : grav_(-1) {
    // Finite element fields -- NOT DIMENSION INDEPENDENT!
    basisPtrVel_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >>();
    basisPtrPrs_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >>();
    basisPtrThr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >>();
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
    // Stochastic Expansion Information.
    Nbottom_ = parlist.sublist("Problem").get("Bottom KL Truncation Order",5);
    Nleft_   = parlist.sublist("Problem").get("Left KL Truncation Order",5);
    Nright_  = parlist.sublist("Problem").get("Right KL Truncation Order",5);
    // Stochastic scales
    ReScale_ = parlist.sublist("Problem").get("Reynolds Number Noise Scale",0.05);
    PrScale_ = parlist.sublist("Problem").get("Prandtl Number Noise Scale",0.05);
    GrScale_ = parlist.sublist("Problem").get("Grashof Number Noise Scale",0.05);
    hScale_  = parlist.sublist("Problem").get("Robin Noise Scale",0.2);
    TScale_  = parlist.sublist("Problem").get("Bottom Temperature Noise Scale",0.2);
    // Pin pressure
    pinPressure_ = parlist.sublist("Problem").get("Pin Pressure",true);

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

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real> > & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c  = u_coeff->dimension(0);
    int p  = cellCub_->getNumPoints();
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize residuals.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > R(d+2);
    for (int i = 0; i < d; ++i) {
      R[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    R[d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    R[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U, Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // Evaluate problem coefficients at quadrature points.
    ROL::Ptr<Intrepid::FieldContainer<Real> > nu, grav, kappa;
    nu    = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    grav  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    kappa = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    computeCoefficients(nu,grav,kappa);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > valVel_vec(d);
    for (int i = 0; i < d; ++i) {
      valVel_vec[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      feVel_->evaluateValue(valVel_vec[i], U[i]);
    }
    ROL::Ptr<Intrepid::FieldContainer<Real> > valPres_eval, valHeat_eval;
    valPres_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    valHeat_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fePrs_->evaluateValue(valPres_eval, U[d]);
    feThr_->evaluateValue(valHeat_eval, U[d+1]);
    // Evaluate/interpolate gradient of finite element fields on cells.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > gradVel_vec(d);
    for (int i = 0; i < d; ++i) {
      gradVel_vec[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      feVel_->evaluateGradient(gradVel_vec[i], U[i]);
    }
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradHeat_eval;
    gradHeat_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    feThr_->evaluateGradient(gradHeat_eval, U[d+1]);

    // Assemble the velocity vector and its divergence.
    ROL::Ptr<Intrepid::FieldContainer<Real> > valVel_eval, divVel_eval;
    valVel_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    divVel_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
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
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > nuGradVel_vec(d), valVelDotgradVel_vec(d);
    for (int i = 0; i < d; ++i) {
      // Multiply velocity gradients with viscosity.
      nuGradVel_vec[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(*nuGradVel_vec[i], *nu, *gradVel_vec[i]);
      // Compute nonlinear terms in the Navier-Stokes equations.
      valVelDotgradVel_vec[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(*valVelDotgradVel_vec[i], *valVel_eval, *gradVel_vec[i]);
    }
    // Negative pressure
    Intrepid::RealSpaceTools<Real>::scale(*valPres_eval,static_cast<Real>(-1));
    // Compute gravitational force.
    ROL::Ptr<Intrepid::FieldContainer<Real> > gravForce;
    gravForce = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*gravForce, *grav, *valHeat_eval);

    /**************************************************************************/
    /*** THERMAL **************************************************************/
    /**************************************************************************/
    ROL::Ptr<Intrepid::FieldContainer<Real> > kappaGradHeat_eval, velDotGradHeat_eval;
    kappaGradHeat_eval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    velDotGradHeat_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
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
    // -->        Control boundaries: i=1,2
    // -->        No slip boundaries: i=0,1,2
    // --> Outflow/Inflow boundaries: i=3
    // -->              Pressure Pin: i=4
    int numSideSets = bdryCellLocIds_.size();
    const int numCubPerSide = bdryCub_->getNumPoints();
    if (numSideSets > 0) {
      // ROBIN CONDITIONS
      for (int i = 0; i < numSideSets; ++i) {
        // Control boundaries
        if (i==1 || i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            if (numCellsSide) {
              // Get U and Z coefficients on Robin boundary
              ROL::Ptr<Intrepid::FieldContainer<Real > > u_coeff_bdry, z_coeff_bdry;
              u_coeff_bdry = getThermalBoundaryCoeff(*(U[d+1]), i, j);
              z_coeff_bdry = getThermalBoundaryCoeff(*(Z[d+1]), i, j);
              // Evaluate U and Z on FE basis
              ROL::Ptr<Intrepid::FieldContainer<Real > > valU_eval_bdry, valZ_eval_bdry;
              valU_eval_bdry = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
              valZ_eval_bdry = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
              feThrBdry_[i][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
              feThrBdry_[i][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
              // Compute Stefan-Boltzmann residual
              ROL::Ptr<Intrepid::FieldContainer<Real> > robinVal;
              robinVal = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
              computeThermalRobin(robinVal,valU_eval_bdry,valZ_eval_bdry,i,j,0);
              // Evaluate boundary integral
              ROL::Ptr<Intrepid::FieldContainer<Real> > robinRes;
              robinRes = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, fh);
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
        if (i<4) {
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
        // Thermal Boundary Conditions
        if ( (i==0) || (i==5) ) {
          int numLocalSideIds = bdryCellLocIds_[i].size();   // Number of sides per cell
          for (int j = 0; j < numLocalSideIds; ++j) {        // Loop over sides of cell: Quad = {0, 1, 2, 3}
            int numCellsSide = bdryCellLocIds_[i][j].size(); // Number of cells with side j
            int numHBdryDofs = fhidx_[j].size();             // Number of thermal boundary DOFs
            for (int k = 0; k < numCellsSide; ++k) {         // Loop over cells with side j
              int cidx = bdryCellLocIds_[i][j][k];           // Cell index
              for (int l = 0; l < numHBdryDofs; ++l) {       // Loop over all fields of cell k on side j
                (*R[d+1])(cidx,fhidx_[j][l]) = (*U[d+1])(cidx,fhidx_[j][l]) - (*bdryCellTDofValues_[i][j])(k,fhidx_[j][l]);
              }
            }
          }
        }
        // Pressure pinning
        if (i==7 && pinPressure_) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int l = 0, cidx = bdryCellLocIds_[i][j][k];
              (*R[d])(cidx,fpidx_[j][l]) = (*U[d])(cidx,fpidx_[j][l]) - static_cast<Real>(0);
            }
          }
        }
      }
    }

    // Combine the residuals.
    fieldHelper_->combineFieldCoeff(res, R);
  }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c  = u_coeff->dimension(0);
    int p  = cellCub_->getNumPoints();
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(d+2);
    for (int i = 0; i < d+2; ++i) {
      J[i].resize(d+2);
    }
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fv);
      }
      J[d][i]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fv);
      J[i][d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fp);
      J[d+1][i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fv);
      J[i][d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fh);
    }
    J[d][d]     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fp);
    J[d+1][d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fp);
    J[d][d+1]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fh);
    J[d+1][d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fh);

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U, Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // Evaluate problem coefficients at quadrature points.
    ROL::Ptr<Intrepid::FieldContainer<Real> > nu, grav, kappa;
    nu    = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    grav  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    kappa = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    computeCoefficients(nu,grav,kappa);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > valVel_vec(d);
    for (int i = 0; i < d; ++i) {
      valVel_vec[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      feVel_->evaluateValue(valVel_vec[i], U[i]);
    }
    // Evaluate/interpolate gradient of finite element fields on cells.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > gradVel_vec(d);
    for (int i = 0; i < d; ++i) {
      gradVel_vec[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      feVel_->evaluateGradient(gradVel_vec[i], U[i]);
    }
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradHeat_eval;
    gradHeat_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    feThr_->evaluateGradient(gradHeat_eval, U[d+1]);

    // Assemble the velocity vector and its divergence.
    ROL::Ptr<Intrepid::FieldContainer<Real> > valVel_eval;
    valVel_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > dVel_vec(d);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > gradHeat_vec(d);
    for (int k = 0; k < d; ++k) {
      dVel_vec[k].resize(d);
      for (int l = 0; l < d; ++l) {
        dVel_vec[k][l] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
        for (int i = 0; i < c; ++i) {
          for (int j = 0; j < p; ++j) {
            (*dVel_vec[k][l])(i,j) = (*gradVel_vec[k])(i,j,l);
          }
        }
      }
      gradHeat_vec[k] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
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
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > dVelPhi_vec(d);
    for (int i = 0; i < d; ++i) {
      // Compute nonlinear terms in the Navier-Stokes equations.
      dVelPhi_vec[i].resize(d);
      for (int j = 0; j < d; ++j) {
        dVelPhi_vec[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p);
        Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*dVelPhi_vec[i][j], *dVel_vec[i][j], *(feVel_->N()));
      }
    }
    // Multiply velocity gradients with viscosity.
    ROL::Ptr<Intrepid::FieldContainer<Real> > nuGradPhi_eval;
    nuGradPhi_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(*nuGradPhi_eval, *nu, *(feVel_->gradN()));
    // Compute nonlinear terms in the Navier-Stokes equations.
    ROL::Ptr<Intrepid::FieldContainer<Real> > valVelDotgradPhi_eval;
    valVelDotgradPhi_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*valVelDotgradPhi_eval, *valVel_eval, *(feVel_->gradN()));
    // Negative pressure basis.
    Intrepid::FieldContainer<Real> negPrsPhi(*(fePrs_->N()));
    Intrepid::RealSpaceTools<Real>::scale(negPrsPhi,static_cast<Real>(-1));
    // Compute gravity Jacobian
    ROL::Ptr<Intrepid::FieldContainer<Real> > gravPhi;
    gravPhi = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*gravPhi, *grav, *(feThr_->N()));

    /**************************************************************************/
    /*** THERMAL **************************************************************/
    /**************************************************************************/
    ROL::Ptr<Intrepid::FieldContainer<Real> > kappaGradPhi, VelPhi;
    kappaGradPhi = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, p, d);
    VelPhi = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, p);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > dHeatPhi(d);
    // Compute kappa times gradient of basis.
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(*kappaGradPhi, *kappa, *(feThr_->gradN()));
    // Compute scaled velocity.
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*VelPhi, *valVel_eval, *feThr_->gradN());
    // Thermal-velocity Jacobians.
    for (int i = 0; i < d; ++i) {
      dHeatPhi[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, p);
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
        if (i==1 || i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            ROL::Ptr<Intrepid::FieldContainer<Real> > robinJac;
            if (numCellsSide) {
              // Get U and Z coefficients on Robin boundary
              ROL::Ptr<Intrepid::FieldContainer<Real > > u_coeff_bdry, z_coeff_bdry;
              u_coeff_bdry = getThermalBoundaryCoeff(*(U[d+1]), i, j);
              z_coeff_bdry = getThermalBoundaryCoeff(*(Z[d+1]), i, j);
              // Evaluate U and Z on FE basis
              ROL::Ptr<Intrepid::FieldContainer<Real > > valU_eval_bdry, valZ_eval_bdry;
              valU_eval_bdry = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
              valZ_eval_bdry = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
              feThrBdry_[i][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
              feThrBdry_[i][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
              // Compute Stefan-Boltzmann residual
              ROL::Ptr< Intrepid::FieldContainer<Real> > robinVal1
                = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
              Intrepid::FieldContainer<Real> robinVal(numCellsSide, fh, numCubPerSide);
              computeThermalRobin(robinVal1,valU_eval_bdry,valZ_eval_bdry,i,j,1,0);
              Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(robinVal,
                                                                          *robinVal1,
                                                                          *(feThrBdry_[i][j]->N()));
              // Evaluate boundary integral
              robinJac = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, fh, fh);
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
        if (i<4) {
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
        // Thermal boundary conditions
        if ( (i==0) || (i==5) ) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numHBdryDofs = fhidx_[j].size();
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
        // Pressure pinning
        if (i==7 && pinPressure_) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int l = 0, cidx = bdryCellLocIds_[i][j][k];
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
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);

  }


  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c  = u_coeff->dimension(0);
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(d+2);
    for (int i = 0; i < d+2; ++i) {
      J[i].resize(d+2);
    }
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fv);
      }
      J[d][i]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fv);
      J[i][d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fp);
      J[d+1][i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fv);
      J[i][d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fh);
    }
    J[d][d]     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fp);
    J[d+1][d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fp);
    J[d][d+1]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fh);
    J[d+1][d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fh);

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U, Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // APPLY DIRICHLET CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    const int numCubPerSide = bdryCub_->getNumPoints();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        // Control boundaries
        if (i==1 || i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            if (numCellsSide) {
              // Get U and Z coefficients on Robin boundary
              ROL::Ptr<Intrepid::FieldContainer<Real > > u_coeff_bdry
                = getThermalBoundaryCoeff(*(U[d+1]), i, j);
              ROL::Ptr<Intrepid::FieldContainer<Real > > z_coeff_bdry
                = getThermalBoundaryCoeff(*(Z[d+1]), i, j);
              // Evaluate U and Z on FE basis
              ROL::Ptr<Intrepid::FieldContainer<Real > > valU_eval_bdry
                = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
              ROL::Ptr<Intrepid::FieldContainer<Real > > valZ_eval_bdry
                = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
              feThrBdry_[i][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
              feThrBdry_[i][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
              // Compute Stefan-Boltzmann residual
              ROL::Ptr< Intrepid::FieldContainer<Real> > robinVal1
                = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
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
        if (i<4) {
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
        // Thermal boundary conditions
        if ( (i==0) || (i==5) ) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numHBdryDofs = fhidx_[j].size();
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
        // Pressure pinning
        if (i==7 && pinPressure_) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int l = 0, cidx = bdryCellLocIds_[i][j][k];
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
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);

  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
//    throw Exception::NotImplemented("");
    // Retrieve dimensions.
    int c  = u_coeff->dimension(0);
    int p  = cellCub_->getNumPoints();
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > H(d+2);
    for (int i = 0; i < d+2; ++i) {
      H[i].resize(d+2);
    }
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        H[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fv);
      }
    }
    for (int i = 0; i < d; ++i) {
      H[d][i]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fv);
      H[i][d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fp);
      H[d+1][i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fv);
      H[i][d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fh);
    }
    H[d][d]     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fp);
    H[d+1][d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fp);
    H[d][d+1]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fh);
    H[d+1][d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fh);

    // Split l_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > L;
    fieldHelper_->splitFieldCoeff(L, l_coeff);

    // APPLY DIRICHLET CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      // DIRICHLET CONDITIONS
      for (int i = 0; i < numSideSets; ++i) {
        // Velocity boundary conditions
        if (i<4) {
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
        // Thermal boundaries
        if ( (i==0) || (i==5) ) {
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
        // Pressure pinning
        if (i==7 && pinPressure_) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int l = 0, cidx = bdryCellLocIds_[i][j][k];
              (*L[d])(cidx,fpidx_[j][l]) = static_cast<Real>(0);
            }
          }
        }
      }
    }

    // Evaluate/interpolate finite element fields on cells.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > valVel_vec(d);
    for (int i = 0; i < d; ++i) {
      valVel_vec[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      feVel_->evaluateValue(valVel_vec[i], L[i]);
    }
    ROL::Ptr<Intrepid::FieldContainer<Real> > valHeat_eval;
    valHeat_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    feThr_->evaluateValue(valHeat_eval, L[3]);

    // Compute nonlinear terms in the Navier-Stokes equations.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > valVelPhi_vec(d);
    for (int i = 0; i < d; ++i) {
      valVelPhi_vec[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p);
      Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*valVelPhi_vec[i], *valVel_vec[i], *(feVel_->N()));
    }

    // Compute nonlinear terms in the thermal equation.
    ROL::Ptr<Intrepid::FieldContainer<Real> > valHeatVelPhi;
    valHeatVelPhi = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p);
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

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    // Retrieve dimensions.
    int c  = feVel_->N()->dimension(0);
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(d+2);
    for (int i = 0; i < d+2; ++i) {
      J[i].resize(d+2);
    }
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fv);
      }
      J[d][i]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fv);
      J[i][d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fp);
      J[d+1][i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fv);
      J[i][d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fh);
    }
    J[d][d]     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fp);
    J[d+1][d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fp);
    J[d][d+1]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fh);
    J[d+1][d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fh);

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

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    // Retrieve dimensions.
    int c  = feVel_->N()->dimension(0);
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(d+2);
    for (int i = 0; i < d+2; ++i) {
      J[i].resize(d+2);
    }
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fv);
      }
      J[d][i]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fv);
      J[i][d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fp);
      J[d+1][i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fv);
      J[i][d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fh);
    }
    J[d][d]     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fp);
    J[d+1][d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fp);
    J[d][d+1]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fh);
    J[d+1][d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fh);


    for (int i = 0; i < d; ++i) {
      *(J[i][i]) = *(feVel_->massMat());
    }
    *(J[d][d]) = *(fePrs_->massMat());
    *(J[d+1][d+1]) = *(feThr_->massMat());

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > getFields() {
    return basisPtrs_;
  }

  void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real> > &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds) {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    feVel_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrVel_,cellCub_);
    fePrs_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrPrs_,cellCub_);
    feThr_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrThr_,cellCub_);
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
        if (bdryCellNodes[i][j] != ROL::nullPtr) {
          feVelBdry_[i][j] = ROL::makePtr<FE<Real>>(bdryCellNodes[i][j],basisPtrVel_,bdryCub_,j);
          fePrsBdry_[i][j] = ROL::makePtr<FE<Real>>(bdryCellNodes[i][j],basisPtrPrs_,bdryCub_,j);
          feThrBdry_[i][j] = ROL::makePtr<FE<Real>>(bdryCellNodes[i][j],basisPtrThr_,bdryCub_,j);
        }
      }
    }
  }

  void setFieldPattern(const std::vector<std::vector<int> > & fieldPattern) {
    fieldPattern_ = fieldPattern;
    fieldHelper_ = ROL::makePtr<FieldHelper<Real>>(numFields_, numDofs_, numFieldDofs_, fieldPattern_);
  }

  const ROL::Ptr<FE<Real> > getVelocityFE(void) const {
    return feVel_;
  }

  const ROL::Ptr<FE<Real> > getPressureFE(void) const {
    return fePrs_;
  }

  const ROL::Ptr<FE<Real> > getThermalFE(void) const {
    return feThr_;
  }

  const std::vector<std::vector<ROL::Ptr<FE<Real> > > > getVelocityBdryFE(void) const {
    return feVelBdry_;
  }

  const std::vector<std::vector<ROL::Ptr<FE<Real> > > > getPressureBdryFE(void) const {
    return fePrsBdry_;
  }

  const std::vector<std::vector<ROL::Ptr<FE<Real> > > > getThermalBdryFE(void) const {
    return feThrBdry_;
  }

  const std::vector<std::vector<std::vector<int> > > getBdryCellLocIds(void) const {
    return bdryCellLocIds_;
  }

  const ROL::Ptr<FieldHelper<Real> > getFieldHelper(void) const {
    return fieldHelper_;
  }

}; // PDE_ThermalFluids

#endif
