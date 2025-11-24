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

#ifndef PDE_THERMALFLUIDSK_HPP
#define PDE_THERMALFLUIDSK_HPP

#include "../TOOLS/pdeK.hpp"
#include "../TOOLS/feK.hpp"
#include "../TOOLS/fieldhelperK.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real, class DeviceType>
class PDE_ThermalFluids : public PDE<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  // Finite element basis information
  basis_ptr basisPtrVel_, basisPtrPrs_, basisPtrThr_;
  std::vector<basis_ptr> basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid2::Cubature<DeviceType,Real,Real>> cellCub_, bdryCub_;
  // Cell node information
  scalar_view volCellNodes_;
  std::vector<std::vector<scalar_view>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<fe_type> feVel_, fePrs_, feThr_;
  std::vector<std::vector<ROL::Ptr<fe_type>>> feVelBdry_, fePrsBdry_, feThrBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fvidx_, fpidx_, fhidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<scalar_view>> bdryCellVDofValues_, bdryCellTDofValues_;
  // Field pattern, offsets, etc.
  std::vector<std::vector<int>> fieldPattern_;  // local Field/DOF pattern; set from DOF manager 
  int numFields_;                                // number of fields (equations in the PDE)
  int numDofs_;                                  // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                      // for each field, a counting offset
  std::vector<int> numFieldDofs_;                // for each field, number of degrees of freedom
  
  // Problem parameters.
  Real Re_, Pr_, Gr_, h_, Tinflow_, Tstep_;
  const Real grav_;

  ROL::Ptr<FieldHelper<Real,DeviceType>> fieldHelper_;

  Real velocityDirichletFunc(const std::vector<Real> & coords, int sideset, int locSideId, int dir) const {
    const std::vector<Real> param = PDE<Real,DeviceType>::getParameter();
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
    const std::vector<Real> param = PDE<Real,DeviceType>::getParameter();
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
        bdryCellVDofValues_[i][j] = scalar_view("bdryCellVDofValues", c, fv, d);
        bdryCellTDofValues_[i][j] = scalar_view("bdryCellTDofValues", c, ft);
        scalar_view Vcoords("Vcoords", c, fv, d);
        scalar_view Tcoords("Tcoords", c, ft, d);
        if (c > 0) {
          feVel_->computeDofCoords(Vcoords, bdryCellNodes_[i][j]);
          feThr_->computeDofCoords(Tcoords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<fv; ++l) {
            std::vector<Real> dofpoint(d);
            //std::cout << "Sideset " << i << " LocalSide " << j << "  Cell " << k << "  Field " << l << "  Coord ";
            for (int m=0; m<d; ++m) {
              dofpoint[m] = Vcoords(k, l, m);
              //std::cout << dofpoint[m] << "  ";
            }

            for (int m=0; m<d; ++m) {
              (bdryCellVDofValues_[i][j])(k, l, m) = velocityDirichletFunc(dofpoint, i, j, m);
              //std::cout << "  " << m << "-Value " << DirichletFunc(dofpoint, i, j, m);
            }
            //std::cout << std::endl;
          }
          for (int l=0; l<ft; ++l) {
            std::vector<Real> dofpoint(d);
            //std::cout << "Sideset " << i << " LocalSide " << j << "  Cell " << k << "  Field " << l << "  Coord ";
            for (int m=0; m<d; ++m) {
              dofpoint[m] = Tcoords(k, l, m);
              //std::cout << dofpoint[m] << "  ";
            }

            (bdryCellTDofValues_[i][j])(k, l) = thermalDirichletFunc(dofpoint, i, j);
            //std::cout << "    Value " << DirichletFunc(dofpoint, i, j);
            //std::cout << std::endl;
          }
        }
      }
    }
  }

  Real ReynoldsNumber(void) const {
    const std::vector<Real> param = PDE<Real,DeviceType>::getParameter();
    if ( param.size() >= 1 ) {
      Real l(190), u(210);
      return (u-l)*param[0] + l;
    }
    return Re_;
  }

  Real PrandtlNumber(void) const {
    const std::vector<Real> param = PDE<Real,DeviceType>::getParameter();
    if ( param.size() >= 2 ) {
      Real l(0.62), u(0.82);
      return (u-l)*param[1] + l;
    }
    return Pr_;
  }

  Real GrashofNumber(void) const {
    const std::vector<Real> param = PDE<Real,DeviceType>::getParameter();
    if ( param.size() >= 3 ) {
      Real l(38000), u(50000);
      return (u-l)*param[2] + l;
    }
    return Gr_;
  }
  
  Real ThicknessNumber(void) const {
    const std::vector<Real> param = PDE<Real,DeviceType>::getParameter();
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

  void computeCoefficients(scalar_view &nu, scalar_view &grav, scalar_view &kappa) const {
    int c = feVel_->gradN().extent_int(0);
    int p = feVel_->gradN().extent_int(2);
    int d = feVel_->gradN().extent_int(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k)
          pt[k] = (feVel_->cubPts())(i,j,k);
        // Compute spatially dependent coefficients
        nu(i,j)    = viscosityFunc(pt);
        grav(i,j)  = gravityFunc(pt);
        kappa(i,j) = conductivityFunc(pt);
      }
    }
  }

  void computeThermalRobin(scalar_view &robin,
                           const scalar_view u,
                           const scalar_view z,
                           int sideset,
                           int locSideId,
                           int deriv = 0,
                           int component = 1) const {
    const int c = u.extent_int(0);
    const int p = u.extent_int(1);
    const int d = feThr_->gradN().extent_int(3);
    std::vector<Real> x(d);
    Real h(0), kappa(0);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k)
          x[k] = (feThrBdry_[sideset][locSideId]->cubPts())(i,j,k);
        h = thermalRobinFunc(x);
        kappa = conductivityFunc(x);
        if ( deriv == 0 )
          robin(i,j) = kappa * h * (u(i,j) - z(i,j));
        else if ( deriv == 1 )
          robin(i,j) = ((component==0) ? kappa * h : -kappa * h);
        else
          robin(i,j) = static_cast<Real>(0);
      }
    }
  }

  scalar_view getThermalBoundaryCoeff(const scalar_view cell_coeff,
                                      int sideSet, int cell) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][cell];
    const int numCellsSide = bdryCellLocId.size();
    const int f = basisPtrThr_->getCardinality();
    
    scalar_view bdry_coeff("bdry_coeff", numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j)
        bdry_coeff(i, j) = cell_coeff(bdryCellLocId[i], j);
    }
    return bdry_coeff;
  }

public:
  PDE_ThermalFluids(ROL::ParameterList &parlist) : grav_(-1) {
    // Finite element fields -- NOT DIMENSION INDEPENDENT!
    basisPtrVel_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<DeviceType,Real,Real>>();
    basisPtrPrs_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType,Real,Real>>();
    basisPtrThr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<DeviceType,Real,Real>>();
    // Volume quadrature rules.
    shards::CellTopology cellType = basisPtrVel_->getBaseCellTopology();     // get the cell type from any basis
    Intrepid2::DefaultCubatureFactory cubFactory;                            // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 4);    // set cubature degree, e.g., 4
    cellCub_ = cubFactory.create<DeviceType,Real,Real>(cellType, cubDegree); // create default cubature
    // Boundary quadrature rules.
    int d = cellType.getDimension();
    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",4); // set cubature degree, e.g., 4
    bdryCub_ = cubFactory.create<DeviceType,Real,Real>(bdryCellType, bdryCubDegree);
    // Fill finite element basis container
    basisPtrs_.clear();
    basisPtrs_.resize(d,basisPtrVel_);  // Velocity
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
      if (i==0)
        offset_[i] = 0;
      else
        offset_[i] = offset_[i-1] + basisPtrs_[i-1]->getCardinality();
      numFieldDofs_[i] = basisPtrs_[i]->getCardinality();
      numDofs_ += numFieldDofs_[i];
    }
  }

  void residual(scalar_view & res,
                const scalar_view u_coeff,
                const scalar_view z_coeff = scalar_view(),
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    int c  = u_coeff.extent_int(0);
    int p  = cellCub_->getNumPoints();
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize residuals.
    std::vector<scalar_view> R(d+2);
    for (int i = 0; i < d; ++i)
      R[i] = scalar_view("res", c, fv);
    R[d]   = scalar_view("res", c, fp);
    R[d+1] = scalar_view("res", c, fh);

    // Split u_coeff into components.
    std::vector<scalar_view> U, Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // Evaluate problem coefficients at quadrature points.
    scalar_view nu("nu", c, p);
    scalar_view grav("grav", c, p);
    scalar_view kappa("kappa", c, p);
    computeCoefficients(nu,grav,kappa);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<scalar_view> valVel_vec(d);
    for (int i = 0; i < d; ++i) {
      valVel_vec[i] = scalar_view("valVel_vec", c, p);
      feVel_->evaluateValue(valVel_vec[i], U[i]);
    }
    scalar_view valPres_eval("valPres_eval", c, p);
    scalar_view valHeat_eval("valHeat_eval", c, p);
    fePrs_->evaluateValue(valPres_eval, U[d]);
    feThr_->evaluateValue(valHeat_eval, U[d+1]);
    // Evaluate/interpolate gradient of finite element fields on cells.
    std::vector<scalar_view> gradVel_vec(d);
    for (int i = 0; i < d; ++i) {
      gradVel_vec[i] = scalar_view("gradVel_vec", c, p, d);
      feVel_->evaluateGradient(gradVel_vec[i], U[i]);
    }
    scalar_view gradHeat_eval("gradHeat_eval", c, p, d);
    feThr_->evaluateGradient(gradHeat_eval, U[d+1]);

    // Assemble the velocity vector and its divergence.
    scalar_view valVel_eval("valVel_eval", c, p, d);
    scalar_view divVel_eval("divVel_eval", c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        divVel_eval(i,j) = static_cast<Real>(0);
        for (int k = 0; k < d; ++k) {
          valVel_eval(i,j,k) = (valVel_vec[k])(i,j);
          divVel_eval(i,j)  += (gradVel_vec[k])(i,j,k);
        }
      }
    }

    /**************************************************************************/
    /*** NAVIER STOKES ********************************************************/
    /**************************************************************************/
    std::vector<scalar_view> nuGradVel_vec(d), valVelDotgradVel_vec(d);
    for (int i = 0; i < d; ++i) {
      // Multiply velocity gradients with viscosity.
      nuGradVel_vec[i] = scalar_view("nuGradVel_vec", c, p, d);
      fst::tensorMultiplyDataData(nuGradVel_vec[i], nu, gradVel_vec[i]);
      // Compute nonlinear terms in the Navier-Stokes equations.
      valVelDotgradVel_vec[i] = scalar_view("valVelDotgradVel_vec", c, p);
      fst::dotMultiplyDataData(valVelDotgradVel_vec[i], valVel_eval, gradVel_vec[i]);
    }
    // Negative pressure
    rst::scale(valPres_eval,static_cast<Real>(-1));
    // Compute gravitational force.
    scalar_view gravForce("gravForce", c, p);
    fst::scalarMultiplyDataData(gravForce, grav, valHeat_eval);

    /**************************************************************************/
    /*** THERMAL **************************************************************/
    /**************************************************************************/
    scalar_view kappaGradHeat_eval("kappaGradHeat_eval", c, p, d);
    scalar_view velDotGradHeat_eval("velDotGradHeat_eval", c, p);
    // Multiply temperature gradient with conductivity
    fst::tensorMultiplyDataData(kappaGradHeat_eval, kappa, gradHeat_eval);
    // Dot multiply scaled velocity with temperature gradient
    fst::dotMultiplyDataData(velDotGradHeat_eval, valVel_eval, gradHeat_eval);

    /**************************************************************************/
    /*** EVALUATE WEAK FORM OF RESIDUAL ***************************************/
    /**************************************************************************/
    // Velocity equation.
    for (int i = 0; i < d; ++i) {
      fst::integrate(R[i],nuGradVel_vec[i],feVel_->gradNdetJ(),false);
      fst::integrate(R[i],valVelDotgradVel_vec[i],feVel_->NdetJ(),true);
      fst::integrate(R[i],valPres_eval,feVel_->DNDdetJ(i),true);
    }
    // Apply gravitational force.
    fst::integrate(R[d-1],gravForce,feVel_->NdetJ(),true);
    // Pressure equation.
    fst::integrate(R[d],divVel_eval,fePrs_->NdetJ(),false);
    rst::scale(R[d],static_cast<Real>(-1));
    // Heat equation.
    fst::integrate(R[d+1],kappaGradHeat_eval,feThr_->gradNdetJ(),false);
    fst::integrate(R[d+1],velDotGradHeat_eval,feThr_->NdetJ(),true);

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
              scalar_view u_coeff_bdry = getThermalBoundaryCoeff(U[d+1], i, j);
              scalar_view z_coeff_bdry = getThermalBoundaryCoeff(Z[d+1], i, j);
              // Evaluate U and Z on FE basis
              scalar_view valU_eval_bdry("valU_eval_bdry", numCellsSide, numCubPerSide);
              scalar_view valZ_eval_bdry("valZ_eval_bdry", numCellsSide, numCubPerSide);
              feThrBdry_[i][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
              feThrBdry_[i][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
              // Compute Stefan-Boltzmann residual
              scalar_view robinVal("robinVal", numCellsSide, numCubPerSide);
              computeThermalRobin(robinVal,valU_eval_bdry,valZ_eval_bdry,i,j,0);
              // Evaluate boundary integral
              scalar_view robinRes("robinRes", numCellsSide, fh);
              fst::integrate(robinRes,robinVal,feThrBdry_[i][j]->NdetJ(),false);
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                // Thermal boundary conditions.
                for (int l = 0; l < fh; ++l)
                  (R[d+1])(cidx,l) += robinRes(k,l);
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
                  (R[m])(cidx,fvidx_[j][l]) = (U[m])(cidx,fvidx_[j][l]) - (bdryCellVDofValues_[i][j])(k,fvidx_[j][l],m);
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
              (R[m])(cidx,fvidx_[j][l]) = (U[m])(cidx,fvidx_[j][l]) - (bdryCellVDofValues_[i][j])(0,fvidx_[j][l],m);
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
                (R[d+1])(cidx,fhidx_[j][l]) = (U[d+1])(cidx,fhidx_[j][l]) - (bdryCellTDofValues_[i][j])(k,fhidx_[j][l]);
              }
            }
          }
        }
        if (i==7) {
          int j = 0, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) { 
            int cidx = bdryCellLocIds_[i][0][0];
            (R[d+1])(cidx,fhidx_[j][l]) = (U[d+1])(cidx,fhidx_[j][l]) - (bdryCellTDofValues_[i][j])(0,fhidx_[j][l]);
          }
        }
        // Pressure pinning
        if (i==9) {
          int j = 2, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) { 
            int cidx = bdryCellLocIds_[i][0][0];
            (R[d])(cidx,fpidx_[j][l]) = (U[d])(cidx,fpidx_[j][l]);
          }
        }
      }
    }

    // Combine the residuals.
    fieldHelper_->combineFieldCoeff(res, R);
  }

  void Jacobian_1(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    int c  = u_coeff.extent_int(0);
    int p  = cellCub_->getNumPoints();
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<scalar_view>> J(d+2);
    for (int i = 0; i < d+2; ++i)
      J[i].resize(d+2);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j)
        J[i][j] = scalar_view("jac", c, fv, fv);
      J[d][i]   = scalar_view("jac", c, fp, fv);
      J[i][d]   = scalar_view("jac", c, fv, fp);
      J[d+1][i] = scalar_view("jac", c, fh, fv);
      J[i][d+1] = scalar_view("jac", c, fv, fh);
    }
    J[d][d]     = scalar_view("jac", c, fp, fp);
    J[d+1][d]   = scalar_view("jac", c, fh, fp);
    J[d][d+1]   = scalar_view("jac", c, fp, fh);
    J[d+1][d+1] = scalar_view("jac", c, fh, fh);

    // Split u_coeff into components.
    std::vector<scalar_view> U, Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // Evaluate problem coefficients at quadrature points.
    scalar_view nu("nu", c, p);
    scalar_view grav("grav", c, p);
    scalar_view kappa("kappa", c, p);
    computeCoefficients(nu,grav,kappa);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<scalar_view> valVel_vec(d);
    for (int i = 0; i < d; ++i) {
      valVel_vec[i] = scalar_view("valVel_vec", c, p);
      feVel_->evaluateValue(valVel_vec[i], U[i]);
    }
    // Evaluate/interpolate gradient of finite element fields on cells.
    std::vector<scalar_view> gradVel_vec(d);
    for (int i = 0; i < d; ++i) {
      gradVel_vec[i] = scalar_view("gradVel_vec", c, p, d);
      feVel_->evaluateGradient(gradVel_vec[i], U[i]);
    }
    scalar_view gradHeat_eval("gradHeat_eval", c, p, d);
    feThr_->evaluateGradient(gradHeat_eval, U[d+1]);

    // Assemble the velocity vector and its divergence.
    scalar_view valVel_eval("valVel_eval", c, p, d);
    std::vector<std::vector<scalar_view>> dVel_vec(d);
    std::vector<scalar_view> gradHeat_vec(d);
    for (int k = 0; k < d; ++k) {
      dVel_vec[k].resize(d);
      for (int l = 0; l < d; ++l) {
        dVel_vec[k][l] = scalar_view("dVel_vec", c, p);
        for (int i = 0; i < c; ++i) {
          for (int j = 0; j < p; ++j)
            (dVel_vec[k][l])(i,j) = (gradVel_vec[k])(i,j,l);
        }
      }
      gradHeat_vec[k] = scalar_view("gradHeat_vec", c, p);
      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          valVel_eval(i,j,k)  = (valVel_vec[k])(i,j);
          (gradHeat_vec[k])(i,j) = gradHeat_eval(i,j,k);
        }
      }
    }

    /**************************************************************************/
    /*** NAVIER STOKES ********************************************************/
    /**************************************************************************/
    std::vector<std::vector<scalar_view>> dVelPhi_vec(d);
    for (int i = 0; i < d; ++i) {
      // Compute nonlinear terms in the Navier-Stokes equations.
      dVelPhi_vec[i].resize(d);
      for (int j = 0; j < d; ++j) {
        dVelPhi_vec[i][j] = scalar_view("dVelPhi_vec", c, fv, p);
        fst::scalarMultiplyDataField(dVelPhi_vec[i][j], dVel_vec[i][j], feVel_->N());
      }
    }
    // Multiply velocity gradients with viscosity.
    scalar_view nuGradPhi_eval("nuGradPhi_eval", c, fv, p, d);
    fst::tensorMultiplyDataField(nuGradPhi_eval, nu, feVel_->gradN());
    // Compute nonlinear terms in the Navier-Stokes equations.
    scalar_view valVelDotgradPhi_eval("valVelDotgradPhi_eval", c, fv, p);
    fst::dotMultiplyDataField(valVelDotgradPhi_eval, valVel_eval, feVel_->gradN());
    // Negative pressure basis.
    scalar_view negPrsPhi("negPrsPhi", c, fp, p);
    Kokkos::deep_copy(negPrsPhi,fePrs_->N());
    rst::scale(negPrsPhi,static_cast<Real>(-1));
    // Compute gravity Jacobian
    scalar_view gravPhi("gravPhi", c, fh, p);
    fst::scalarMultiplyDataField(gravPhi, grav, feThr_->N());

    /**************************************************************************/
    /*** THERMAL **************************************************************/
    /**************************************************************************/
    scalar_view kappaGradPhi("kappaGradPhi", c, fh, p, d);
    scalar_view VelPhi("velPhi", c, fh, p);
    std::vector<scalar_view> dHeatPhi(d);
    // Compute kappa times gradient of basis.
    fst::tensorMultiplyDataField(kappaGradPhi, kappa, feThr_->gradN());
    // Compute scaled velocity.
    fst::dotMultiplyDataField(VelPhi, valVel_eval, feThr_->gradN());
    // Thermal-velocity Jacobians.
    for (int i = 0; i < d; ++i) {
      dHeatPhi[i] = scalar_view("dHeatPhi", c, fh, p);
      fst::scalarMultiplyDataField(dHeatPhi[i], gradHeat_vec[i], feThr_->N());
    }

    /*** Evaluate weak form of the Jacobian. ***/
    // X component of velocity equation.
    for (int i = 0; i < d; ++i) {
      // Velocity components
      for (int j = 0; j < d; ++j)
        fst::integrate(J[i][j],feVel_->NdetJ(),dVelPhi_vec[i][j],false);
      fst::integrate(J[i][i],nuGradPhi_eval,feVel_->gradNdetJ(),true);
      fst::integrate(J[i][i],feVel_->NdetJ(),valVelDotgradPhi_eval,true);
      // Pressure components
      fst::integrate(J[i][d],feVel_->DNDdetJ(i),negPrsPhi,false);
      fst::integrate(J[d][i],fePrs_->NdetJ(),feVel_->DND(i),false);
      rst::scale(J[d][i],static_cast<Real>(-1));
      // Thermal components
      fst::integrate(J[d+1][i],dHeatPhi[i],feVel_->NdetJ(),false);
    }
    // Gravitational component
    fst::integrate(J[d-1][d+1],feVel_->NdetJ(),gravPhi,false);
    // Thermal components
    fst::integrate(J[d+1][d+1],kappaGradPhi,feThr_->gradNdetJ(),false);
    fst::integrate(J[d+1][d+1],feThr_->NdetJ(),VelPhi,true);

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
            if (numCellsSide) {
              // Get U and Z coefficients on Robin boundary
              scalar_view u_coeff_bdry = getThermalBoundaryCoeff(U[d+1], i, j);
              scalar_view z_coeff_bdry = getThermalBoundaryCoeff(Z[d+1], i, j);
              // Evaluate U and Z on FE basis
              scalar_view valU_eval_bdry("valU_eval_bdry", numCellsSide, numCubPerSide);
              scalar_view valZ_eval_bdry("valZ_eval_bdry", numCellsSide, numCubPerSide);
              feThrBdry_[i][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
              feThrBdry_[i][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
              // Compute Stefan-Boltzmann residual
              scalar_view robinVal1("robinVal1", numCellsSide, numCubPerSide);
              scalar_view robinVal("robinVal", numCellsSide, fh, numCubPerSide);
              computeThermalRobin(robinVal1,valU_eval_bdry,valZ_eval_bdry,i,j,1,0);
              fst::scalarMultiplyDataField(robinVal,robinVal1,feThrBdry_[i][j]->N());
              // Evaluate boundary integral
              scalar_view robinJac("robinJac", numCellsSide, fh, fh);
              fst::integrate(robinJac,robinVal,feThrBdry_[i][j]->NdetJ(),false);
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < fh; ++l) {
                  for (int m = 0; m < fh; ++m) {
                    (J[d+1][d+1])(cidx,l,m) += robinJac(k,l,m);
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
                    for (int q = 0; q < d; ++q)
                      (J[n][q])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                    (J[n][n])(cidx,fvidx_[j][l],fvidx_[j][l]) = static_cast<Real>(1);
                  }
                }
                for (int m = 0; m < d; ++m) {
                  for (int n = 0; n < fp; ++n)
                    (J[m][d])(cidx,fvidx_[j][l],n) = static_cast<Real>(0);
                  for (int n = 0; n < fh; ++n)
                    (J[m][d+1])(cidx,fvidx_[j][l],n) = static_cast<Real>(0);
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
                for (int q=0; q < d; ++q)
                  (J[n][q])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                (J[n][n])(cidx,fvidx_[j][l],fvidx_[j][l]) = static_cast<Real>(1);
              }
            }
            for (int m=0; m < d; ++m) {
              for (int n=0; n < fp; ++n)
                (J[m][d])(cidx,fvidx_[j][l],n) = static_cast<Real>(0);
              for (int n = 0; n < fh; ++n)
                (J[m][d+1])(cidx,fvidx_[j][l],n) = static_cast<Real>(0);
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
                  for (int n = 0; n < d; ++n)
                    (J[d+1][n])(cidx,fhidx_[j][l],m) = static_cast<Real>(0);
                }
                for (int m = 0; m < fp; ++m)
                  (J[d+1][d])(cidx,fhidx_[j][l],m) = static_cast<Real>(0);
                for (int m = 0; m < fh; ++m)
                  (J[d+1][d+1])(cidx,fhidx_[j][l],m) = static_cast<Real>(0);
                (J[d+1][d+1])(cidx,fhidx_[j][l],fhidx_[j][l]) = static_cast<Real>(1);
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
                for (int q=0; q < d; ++q)
                  (J[d+1][q])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
              }
            }
            for (int m=0; m < fp; ++m)
              (J[d+1][d])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
            for (int m = 0; m < fh; ++m)
              (J[d+1][d+1])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
            (J[d+1][d+1])(cidx,fhidx_[j][l],fhidx_[j][l]) = static_cast<Real>(1);
          }
        }
        // Pressure pinning
        if (i==9) {
          int j = 2, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) {
            int cidx = bdryCellLocIds_[i][0][0];
            for (int m = 0; m < fv; ++m) {
              for (int n = 0; n < d; ++n)
                (J[d][n])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
            }
            for (int m = 0; m < fp; ++m)
              (J[d][d])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
            (J[d][d])(cidx,fpidx_[j][l],fpidx_[j][l]) = static_cast<Real>(1);
            for (int m = 0; m < fh; ++m)
              (J[d][d+1])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
          }
        }
      }
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }


  void Jacobian_2(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    int c  = u_coeff.extent_int(0);
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<scalar_view>> J(d+2);
    for (int i = 0; i < d+2; ++i)
      J[i].resize(d+2);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j)
        J[i][j] = scalar_view("jac", c, fv, fv);
      J[d][i]   = scalar_view("jac", c, fp, fv);
      J[i][d]   = scalar_view("jac", c, fv, fp);
      J[d+1][i] = scalar_view("jac", c, fh, fv);
      J[i][d+1] = scalar_view("jac", c, fv, fh);
    }
    J[d][d]     = scalar_view("jac", c, fp, fp);
    J[d+1][d]   = scalar_view("jac", c, fh, fp);
    J[d][d+1]   = scalar_view("jac", c, fp, fh);
    J[d+1][d+1] = scalar_view("jac", c, fh, fh);

    // Split u_coeff into components.
    std::vector<scalar_view> U, Z;
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
              scalar_view u_coeff_bdry = getThermalBoundaryCoeff(U[d+1], i, j);
              scalar_view z_coeff_bdry = getThermalBoundaryCoeff(Z[d+1], i, j);
              // Evaluate U and Z on FE basis
              scalar_view valU_eval_bdry("valU_eval_bdry", numCellsSide, numCubPerSide);
              scalar_view valZ_eval_bdry("valZ_eval_bdry", numCellsSide, numCubPerSide);
              feThrBdry_[i][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
              feThrBdry_[i][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
              // Compute Stefan-Boltzmann residual
              scalar_view robinVal1("robinVal1", numCellsSide, numCubPerSide);
              scalar_view robinVal("robinVal", numCellsSide, fh, numCubPerSide);
              computeThermalRobin(robinVal1,valU_eval_bdry,valZ_eval_bdry,i,j,1,1);
              fst::scalarMultiplyDataField(robinVal,robinVal1,feThrBdry_[i][j]->N());
              // Evaluate boundary integral
              scalar_view robinJac("robinJac", numCellsSide, fh, fh);
              fst::integrate(robinJac,robinVal,feThrBdry_[i][j]->NdetJ(),false);
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < fh; ++l) {
                  for (int m = 0; m < fh; ++m)
                    (J[d+1][d+1])(cidx,l,m) += robinJac(k,l,m);
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
                    for (int p = 0; p < d; ++p)
                      (J[n][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
                for (int m = 0; m < d; ++m) {
                  for (int n = 0; n < fp; ++n)
                    (J[m][d])(cidx,fvidx_[j][l],n) = static_cast<Real>(0);
                  for (int n = 0; n < fh; ++n)
                    (J[m][d+1])(cidx,fvidx_[j][l],n) = static_cast<Real>(0);
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
                for (int p=0; p < d; ++p)
                  (J[n][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
              }
            }
            for (int m=0; m < d; ++m) {
              for (int n=0; n < fp; ++n)
                (J[m][d])(cidx,fvidx_[j][l],n) = static_cast<Real>(0);
              for (int n = 0; n < fh; ++n)
                (J[m][d+1])(cidx,fvidx_[j][l],n) = static_cast<Real>(0);
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
                  for (int n = 0; n < d; ++n)
                    (J[d+1][n])(cidx,fhidx_[j][l],m) = static_cast<Real>(0);
                }
                for (int m = 0; m < fp; ++m)
                  (J[d+1][d])(cidx,fhidx_[j][l],m) = static_cast<Real>(0);
                for (int m = 0; m < fh; ++m)
                  (J[d+1][d+1])(cidx,fhidx_[j][l],m) = static_cast<Real>(0);
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
                for (int p=0; p < d; ++p)
                  (J[d+1][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
              }
            }
            for (int m=0; m < fp; ++m)
              (J[d+1][d])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
            for (int m = 0; m < fh; ++m)
              (J[d+1][d+1])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
          }
        }
        // Pressure pinning
        if (i==9) {
          int j = 2, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) {
            int cidx = bdryCellLocIds_[i][0][0];
            for (int m = 0; m < fv; ++m) {
              for (int n = 0; n < d; ++n)
                (J[d][n])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
            }
            for (int m = 0; m < fp; ++m)
              (J[d][d])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
            for (int m = 0; m < fh; ++m)
              (J[d][d+1])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
          }
        }
      }
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }

  void Hessian_11(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
//    throw Exception::NotImplemented("");
    // Retrieve dimensions.
    int c  = u_coeff.extent_int(0);
    int p  = cellCub_->getNumPoints();
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<scalar_view>> H(d+2);
    for (int i = 0; i < d+2; ++i)
      H[i].resize(d+2);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j)
        H[i][j] = scalar_view("hess", c, fv, fv);
    }
    for (int i = 0; i < d; ++i) {
      H[d][i]   = scalar_view("hess", c, fp, fv);
      H[i][d]   = scalar_view("hess", c, fv, fp);
      H[d+1][i] = scalar_view("hess", c, fh, fv);
      H[i][d+1] = scalar_view("hess", c, fv, fh);
    }
    H[d][d]     = scalar_view("hess", c, fp, fp);
    H[d+1][d]   = scalar_view("hess", c, fh, fp);
    H[d][d+1]   = scalar_view("hess", c, fp, fh);
    H[d+1][d+1] = scalar_view("hess", c, fh, fh);

    // Split l_coeff into components.
    std::vector<scalar_view> L;
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
                  (L[m])(cidx,fvidx_[j][l]) = static_cast<Real>(0);
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
              (L[m])(cidx,fvidx_[j][l]) = static_cast<Real>(0);
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
                (L[d+1])(cidx,fhidx_[j][l]) = static_cast<Real>(0);
              }
            }
          }
        }
        if (i==7) {
          int j = 0, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) { 
            int cidx = bdryCellLocIds_[i][0][0];
            (L[d+1])(cidx,fhidx_[j][l]) = static_cast<Real>(0);
          }
        }
        // Pressure pinning
        if (i==9) {
          int j = 2, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) { 
            int cidx = bdryCellLocIds_[i][0][0];
            (L[d])(cidx,fpidx_[j][l]) = static_cast<Real>(0);
          }
        }
      }
    }

    // Evaluate/interpolate finite element fields on cells.
    std::vector<scalar_view> valVel_vec(d);
    for (int i = 0; i < d; ++i) {
      valVel_vec[i] = scalar_view("valVel_vec", c, p);
      feVel_->evaluateValue(valVel_vec[i], L[i]);
    }
    scalar_view valHeat_eval("valHeat_eval", c, p);
    feThr_->evaluateValue(valHeat_eval, L[3]);

    // Compute nonlinear terms in the Navier-Stokes equations.
    std::vector<scalar_view> valVelPhi_vec(d);
    for (int i = 0; i < d; ++i) {
      valVelPhi_vec[i] = scalar_view("valVelPhi_vec", c, fv, p);
      fst::scalarMultiplyDataField(valVelPhi_vec[i], valVel_vec[i], feVel_->N());
    }

    // Compute nonlinear terms in the thermal equation.
    scalar_view valHeatVelPhi("valHeatVelPhi", c, fv, p);
    fst::scalarMultiplyDataField(valHeatVelPhi, valHeat_eval, feVel_->N());

    /*** Evaluate weak form of the Hessian. ***/
    for (int i = 0; i < d; ++i) {
      // velocity equation.
      for (int j = 0; j < d; ++j) {
        fst::integrate(H[i][j],valVelPhi_vec[j],feVel_->DNDdetJ(i),false);
        fst::integrate(H[i][j],feVel_->DNDdetJ(j),valVelPhi_vec[i],true);
      }
      fst::integrate(H[i][d+1],valHeatVelPhi,feThr_->DNDdetJ(i),false);
      // Thermal equation.
      fst::integrate(H[d+1][i],feThr_->DNDdetJ(i),valHeatVelPhi,false);
    }

    // Combine the Hessians.
    fieldHelper_->combineFieldCoeff(hess, H);
  }

  void Hessian_12(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_ThermalFluids::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(scalar_view & riesz) override {
    // Retrieve dimensions.
    int c  = feVel_->N().extent_int(0);
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<scalar_view>> J(d+2);
    for (int i = 0; i < d+2; ++i)
      J[i].resize(d+2);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j)
        J[i][j] = scalar_view("riesz1", c, fv, fv);
      J[d][i]   = scalar_view("riesz1", c, fp, fv);
      J[i][d]   = scalar_view("riesz1", c, fv, fp);
      J[d+1][i] = scalar_view("riesz1", c, fh, fv);
      J[i][d+1] = scalar_view("riesz1", c, fv, fh);
    }
    J[d][d]     = scalar_view("riesz1", c, fp, fp);
    J[d+1][d]   = scalar_view("riesz1", c, fh, fp);
    J[d][d+1]   = scalar_view("riesz1", c, fp, fh);
    J[d+1][d+1] = scalar_view("riesz1", c, fh, fh);

    for (int i = 0; i < d; ++i) {
      Kokkos::deep_copy(J[i][i],feVel_->stiffMat());
      rst::add(J[i][i],feVel_->massMat());
    }
    Kokkos::deep_copy(J[d][d],fePrs_->massMat());
    Kokkos::deep_copy(J[d+1][d+1],feThr_->stiffMat());
    rst::add(J[d+1][d+1],feThr_->massMat());

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  void RieszMap_2(scalar_view & riesz) override {
    // Retrieve dimensions.
    int c  = feVel_->N().extent_int(0);
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<scalar_view>> J(d+2);
    for (int i = 0; i < d+2; ++i)
      J[i].resize(d+2);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j)
        J[i][j] = scalar_view("reisz2", c, fv, fv);
      J[d][i]   = scalar_view("reisz2", c, fp, fv);
      J[i][d]   = scalar_view("reisz2", c, fv, fp);
      J[d+1][i] = scalar_view("reisz2", c, fh, fv);
      J[i][d+1] = scalar_view("reisz2", c, fv, fh);
    }
    J[d][d]     = scalar_view("reisz2", c, fp, fp);
    J[d+1][d]   = scalar_view("reisz2", c, fh, fp);
    J[d][d+1]   = scalar_view("reisz2", c, fp, fh);
    J[d+1][d+1] = scalar_view("reisz2", c, fh, fh);

    for (int i = 0; i < d; ++i)
      Kokkos::deep_copy(J[i][i],feVel_->massMat());
    Kokkos::deep_copy(J[d][d],fePrs_->massMat());
    Kokkos::deep_copy(J[d+1][d+1],feThr_->massMat());

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  std::vector<basis_ptr> getFields() override {
    return basisPtrs_;
  }

  void setCellNodes(const scalar_view &volCellNodes,
                    const std::vector<std::vector<scalar_view>> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds) override {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    feVel_ = ROL::makePtr<fe_type>(volCellNodes_,basisPtrVel_,cellCub_);
    fePrs_ = ROL::makePtr<fe_type>(volCellNodes_,basisPtrPrs_,cellCub_);
    feThr_ = ROL::makePtr<fe_type>(volCellNodes_,basisPtrThr_,cellCub_);
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
        if (bdryCellNodes[i][j] != scalar_view()) {
          feVelBdry_[i][j] = ROL::makePtr<fe_type>(bdryCellNodes[i][j],basisPtrVel_,bdryCub_,j);
          fePrsBdry_[i][j] = ROL::makePtr<fe_type>(bdryCellNodes[i][j],basisPtrPrs_,bdryCub_,j);
          feThrBdry_[i][j] = ROL::makePtr<fe_type>(bdryCellNodes[i][j],basisPtrThr_,bdryCub_,j);
        }
      }
    }
  }

  void setFieldPattern(const std::vector<std::vector<int>> & fieldPattern) override {
    fieldPattern_ = fieldPattern;
    fieldHelper_ = ROL::makePtr<FieldHelper<Real,DeviceType>>(numFields_, numDofs_, numFieldDofs_, fieldPattern_);
  }

  const ROL::Ptr<fe_type> getVelocityFE(void) const {
    return feVel_;
  }

  const ROL::Ptr<fe_type> getPressureFE(void) const {
    return fePrs_;
  }

  const ROL::Ptr<fe_type> getThermalFE(void) const {
    return feThr_;
  }

  const std::vector<std::vector<ROL::Ptr<fe_type>>> getVelocityBdryFE(void) const {
    return feVelBdry_;
  }

  const std::vector<std::vector<ROL::Ptr<fe_type>>> getPressureBdryFE(void) const {
    return fePrsBdry_;
  }

  const std::vector<std::vector<ROL::Ptr<fe_type>>> getThermalBdryFE(void) const {
    return feThrBdry_;
  }

  const std::vector<std::vector<std::vector<int>>> getBdryCellLocIds(void) const {
    return bdryCellLocIds_;
  }

  const ROL::Ptr<FieldHelper<Real,DeviceType>> getFieldHelper(void) const {
    return fieldHelper_;
  }

}; // PDE_ThermalFluids

#endif
