// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde_brinkman.hpp
    \brief Implements the local PDE interface for the Brinkman porosity optimization problem.

    Implementes the weak form of the Navier-Stokes-Brinkman equations
    \f[
       \begin{aligned}
         -\nabla\cdot(\mu\nabla u) + \rho (u\cdot\nabla) u + \nabla p &= -\mu K(z)^{-1}u \\
         \nabla\cdot u &=0
    \f]
*/

#ifndef PDE_BRINKMANK_HPP
#define PDE_BRINKMANK_HPP

#include "../../../../TOOLS/pdeK.hpp"
#include "../../../../TOOLS/feK.hpp"
#include "../../../../TOOLS/fieldhelperK.hpp"
#include "impermeabilityK.hpp"

#include "Intrepid2_HVOL_C0_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_C2_FEM.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real, class DeviceType>
class PDE_NavierStokes : public PDE<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  // Finite element basis information
  basis_ptr basisPtrVel_, basisPtrPrs_, basisPtrCtrl_;
  std::vector<basis_ptr> basisPtrs_, basisPtrsCtrl_;
  // Cell cubature information
  ROL::Ptr<Intrepid2::Cubature<DeviceType,Real,Real>> cellCub_, bdryCub_;
  // Cell node information
  scalar_view volCellNodes_;
  std::vector<std::vector<scalar_view>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<fe_type> feVel_, fePrs_, feCtrl_;
  std::vector<std::vector<ROL::Ptr<fe_type>>> feVelBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fvidx_, fpidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<scalar_view>> bdryCellDofValues_;
  // Field pattern, offsets, etc.
  std::vector<std::vector<int>> fieldPattern_;  // local Field/DOF pattern; set from DOF manager 
  int numFields_;                               // number of fields (equations in the PDE)
  int numDofs_;                                 // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                     // for each field, a counting offset
  std::vector<int> numFieldDofs_;               // for each field, number of degrees of freedom
  std::vector<std::vector<int>> fieldPatternCtrl_;  // local Field/DOF pattern; set from DOF manager 
  int numFieldsCtrl_;                               // number of fields (equations in the PDE)
  int numDofsCtrl_;                                 // total number of degrees of freedom for all (local) fields
  std::vector<int> offsetCtrl_;                     // for each field, a counting offset
  std::vector<int> numFieldDofsCtrl_;               // for each field, number of degrees of freedom
  
  // Problem parameters.
  Real Patm_, density_, viscosity_, inVelocity_;
  bool useStokes_, useNoSlip_, usePresOut_;

  ROL::Ptr<Impermeability<Real,DeviceType>> imp_;
  ROL::Ptr<FieldHelper<Real,DeviceType>> fieldHelper_;
  ROL::Ptr<FieldUtils::FieldInfo> fieldInfo_, fieldInfoCtrl_;

  // Extract velocity coefficients on boundary.
  scalar_view getBoundaryCoeff(const scalar_view cell_coeff, int sideSet, int cell) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][cell];
    const int numCellsSide = bdryCellLocId.size();
    const int f = basisPtrVel_->getCardinality();
    
    scalar_view bdry_coeff("bdry_coeff", numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        bdry_coeff(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

  Real DirichletFunc(const std::vector<Real> & coords, int sideset, int locSideId, int dir) const {
    const std::vector<Real> param = PDE<Real,DeviceType>::getParameter();
    Real val(0);
    // sideset == 0: Inflow
    // sideset == 1: No Slip
    // sideset == 2: Outflow
    // sideset == 3: No Normal
    if (sideset==0) val = (dir==1 ? inVelocity_ : static_cast<Real>(0));
    return val;
  }

  void computeDirichlet(void) {
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
        bdryCellDofValues_[i][j] = scalar_view("bdryCellDofValues", c, f, d);
        scalar_view coords("coords", c, f, d);
        if (c > 0) {
          feVel_->computeDofCoords(coords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<f; ++l) {
            std::vector<Real> dofpoint(d);
            //std::cout << "Sideset " << i << " LocalSide " << j << "  Cell " << k << "  Field " << l << "  Coord ";
            for (int m=0; m<d; ++m) {
              dofpoint[m] = coords(k, l, m);
              //std::cout << dofpoint[m] << "  ";
            }

            for (int m=0; m<d; ++m) {
              (bdryCellDofValues_[i][j])(k, l, m) = DirichletFunc(dofpoint, i, j, m);
              //std::cout << "  " << m << "-Value " << DirichletFunc(dofpoint, i, j, m);
            }
            //std::cout << std::endl;
          }
        }
      }
    }
  }

  Real viscosityFunc(const std::vector<Real> & coords) const {
    return viscosity_;
  }

  void computeViscosity(scalar_view &nu) const {
    int c = feVel_->gradN().extent_int(0);
    int p = feVel_->gradN().extent_int(2);
    int d = feVel_->gradN().extent_int(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (feVel_->cubPts())(i,j,k);
        }
        // Compute spatially dependent viscosity
        nu(i,j) = viscosityFunc(pt);
      }
    }
  }

  void multiplyByRadius(scalar_view &input, bool isField, bool useReciprocal = false) const {
    int c = feVel_->gradN().extent_int(0);
    int p = feVel_->gradN().extent_int(2);
    scalar_view r("r",c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        if (!useReciprocal) r(i,j) = (feVel_->cubPts())(i,j,0);
        else                r(i,j) = 1/(feVel_->cubPts())(i,j,0);
      }
    }
    scalar_view in(input);
    if (!isField)
      fst::scalarMultiplyDataData(input,r,in);
    else
      fst::scalarMultiplyDataField(input,r,in);
  }

public:
  PDE_NavierStokes(ROL::ParameterList &parlist) {
    // Finite element fields.
    int cubDegree     = parlist.sublist("Problem").get("Cubature Degree",4);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2);
    int probDim       = parlist.sublist("Problem").get("Problem Dimension",2);
    int basisDegCtrl  = parlist.sublist("Problem").get("Density Basis Degree",1);
    std::string elemType = parlist.sublist("Problem").get("Element Type","QUAD");
    TEUCHOS_TEST_FOR_EXCEPTION(probDim!=2, std::invalid_argument,
      ">>> rol/example/PDE-OPT/flow-opt/axisymmetric/models/brinkman/pde_brinkman.hpp: Problem dimension is not 2!");
    if (elemType == "TRI") {
      basisPtrVel_ = ROL::makePtr<Intrepid2::Basis_HGRAD_TRI_C2_FEM<DeviceType,Real,Real>>();
      basisPtrPrs_ = ROL::makePtr<Intrepid2::Basis_HGRAD_TRI_C1_FEM<DeviceType,Real,Real>>();
    }
    else {
      basisPtrVel_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<DeviceType,Real,Real>>();
      basisPtrPrs_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType,Real,Real>>();
    }
    basisPtrs_.clear();
    basisPtrs_.resize(probDim,basisPtrVel_); // Velocity component
    basisPtrs_.push_back(basisPtrPrs_);      // Pressure component
    shards::CellTopology cellType = basisPtrs_[0]->getBaseCellTopology(); // get the cell type from any basis
    if (basisDegCtrl == 1) {
      if (elemType == "TRI")
        basisPtrCtrl_ = ROL::makePtr<Intrepid2::Basis_HGRAD_TRI_C1_FEM<DeviceType,Real,Real>>();
      else
        basisPtrCtrl_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType,Real,Real>>();
    }
    else
      basisPtrCtrl_ = ROL::makePtr<Intrepid2::Basis_HVOL_C0_FEM<DeviceType,Real,Real>>(cellType);
    basisPtrsCtrl_.clear();
    basisPtrsCtrl_.push_back(basisPtrCtrl_); // Control component
    // Quadrature rules.
    Intrepid2::DefaultCubatureFactory cubFactory;                            // create cubature factory
    cellCub_ = cubFactory.create<DeviceType,Real,Real>(cellType, cubDegree); // create default cubature

    shards::CellTopology bdryCellType = cellType.getCellTopologyData(probDim-1, 0);
    bdryCub_ = cubFactory.create<DeviceType,Real,Real>(bdryCellType, bdryCubDegree);

    // Other problem parameters.
    Real dynvisco = parlist.sublist("Problem").get("Dynamic Viscosity", 0.84e-8); // kg/mm-s
    Real inRadius = parlist.sublist("Problem").get("Inlet Radius", 0.6875);       // mm
    Real Reynolds = parlist.sublist("Problem").get("Reynolds Number", 50.0);
    Patm_       = static_cast<Real>(101.325);
    density_    = parlist.sublist("Problem").get("Fluid Density", 8.988e-11);     // kg/mm^3
    viscosity_  = parlist.sublist("Problem").get("Effective Viscosity", 0.84e-8); // kg/mm-s
    inVelocity_ = Reynolds * dynvisco / (static_cast<Real>(2) * inRadius * density_);
    useStokes_  = parlist.sublist("Problem").get("Use Stokes", false);
    useNoSlip_  = parlist.sublist("Problem").get("Use No Slip", true);
    usePresOut_ = parlist.sublist("Problem").get("Use Pressure Outflow",false);
    imp_        = ROL::makePtr<Impermeability<Real,DeviceType>>(parlist);

    numDofs_ = 0;
    numFields_ = basisPtrs_.size();
    offset_.resize(numFields_);
    numFieldDofs_.resize(numFields_);
    for (int i=0; i<numFields_; ++i) {
      offset_[i] = i==0 ? 0 : offset_[i-1] + basisPtrs_[i-1]->getCardinality();
      numFieldDofs_[i] = basisPtrs_[i]->getCardinality();
      numDofs_ += numFieldDofs_[i];
    }
    numDofsCtrl_ = 0;
    numFieldsCtrl_ = basisPtrsCtrl_.size();
    offsetCtrl_.resize(numFieldsCtrl_);
    numFieldDofsCtrl_.resize(numFieldsCtrl_);
    for (int i=0; i<numFieldsCtrl_; ++i) {
      offsetCtrl_[i] = i==0 ? 0 : offsetCtrl_[i-1] + basisPtrsCtrl_[i-1]->getCardinality();
      numFieldDofsCtrl_[i] = basisPtrsCtrl_[i]->getCardinality();
      numDofsCtrl_ += numFieldDofsCtrl_[i];
    }
  }

  void residual(scalar_view & res,
                const scalar_view u_coeff,
                const scalar_view z_coeff = scalar_view(),
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c  = feVel_->gradN().extent_int(0);
    const int fv = feVel_->gradN().extent_int(1);
    const int fp = fePrs_->gradN().extent_int(1);
    const int p  = feVel_->gradN().extent_int(2);
    const int d  = feVel_->gradN().extent_int(3);
    // Initialize residuals.
    std::vector<scalar_view> R;
    R.resize(numFields_);
    for (int i = 0; i < d; ++i)
      R[i] = scalar_view("res", c, fv);
    R[d] = scalar_view("res", c, fp);
    // Split coefficients into components.
    std::vector<scalar_view> U, Z;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Evaluate/interpolate finite element fields on cells.
    std::vector<scalar_view> nuGradVel(d), valVel(d), gradVel(d), valVelDotgradVel(d), alphaVel(d);
    scalar_view nu("nu", c, p);
    scalar_view valVel0("valVel0", c, p, d);
    scalar_view valPres("valPres", c, p);
    scalar_view divVel("divVel", c, p);
    scalar_view valCtrl("valCtrl", c, p);
    scalar_view alpha("alpha", c, p);
    scalar_view nuVel0("nuVel0", c, p);
    computeViscosity(nu);
    feCtrl_->evaluateValue(valCtrl, Z[0]);
    imp_->compute(alpha, valCtrl, feVel_->cubPts(), 0);
    for (int i = 0; i < d; ++i) {
      nuGradVel[i]        = scalar_view("nuGradVel", c, p, d);
      valVel[i]           = scalar_view("valVel", c, p);
      gradVel[i]          = scalar_view("gradVel", c, p, d);
      valVelDotgradVel[i] = scalar_view("valVelDotgradVel", c, p);
      alphaVel[i]         = scalar_view("alphaVel", c, p);
      // Evaluate on FE basis
      feVel_->evaluateValue(valVel[i], U[i]);
      feVel_->evaluateGradient(gradVel[i], U[i]);
      // Multiply velocity gradients with viscosity.
      fst::tensorMultiplyDataData(nuGradVel[i], nu, gradVel[i]);
      // Multiply velocity with alpha
      fst::scalarMultiplyDataData(alphaVel[i], alpha, valVel[i]);
      // Assemble the velocity vector and its divergence.
      for (int j = 0; j < c; ++j) {
        for (int k = 0; k < p; ++k) {
          valVel0(j,k,i) = (valVel[i])(j,k);
          divVel(j,k)   += (gradVel[i])(j,k,i);
        }
      }
    }
    fst::scalarMultiplyDataData(nuVel0, nu, valVel[0]);
    // Compute nonlinear terms in the Navier-Stokes equations.
    for (int i = 0; i < d; ++i) {
      fst::dotMultiplyDataData(valVelDotgradVel[i], valVel0, gradVel[i]);
      rst::scale(valVelDotgradVel[i], density_);
    }
    // Negative pressure
    fePrs_->evaluateValue(valPres, U[d]);
    rst::scale(valPres, static_cast<Real>(-1));
    /*** Evaluate weak form of the residual. ***/
    // Integrate velocity equation.
    fst::integrate(R[0],valPres,feVel_->NdetJ(),false);
    multiplyByRadius(nuVel0,false,true);
    fst::integrate(R[0],nuVel0,feVel_->NdetJ(),true);
    multiplyByRadius(valPres,false);
    for (int i = 0; i < d; ++i) {
      multiplyByRadius(nuGradVel[i],false);
      fst::integrate(R[i],nuGradVel[i],feVel_->gradNdetJ(),true);
      if (!useStokes_) {
        multiplyByRadius(valVelDotgradVel[i],false);
        fst::integrate(R[i],valVelDotgradVel[i],feVel_->NdetJ(),true);
      }
      fst::integrate(R[i],valPres,feVel_->DNDdetJ(i),true);
      multiplyByRadius(alphaVel[i],false);
      fst::integrate(R[i],alphaVel[i],feVel_->NdetJ(),true);
    }
    // Integrate pressure equation.
    multiplyByRadius(divVel,false);
    fst::integrate(R[d],divVel,fePrs_->NdetJ(),false);
    fst::integrate(R[d],valVel[0],fePrs_->NdetJ(),true);
    rst::scale(R[d], static_cast<Real>(-1));
    // Boundary conditions.
    applyDirichletToResidual(R,U);
    // Combine the residuals.
    FieldUtils::combineFieldCoeff<Real>(res, R, fieldInfo_);
  }

  void Jacobian_1(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c  = feVel_->gradN().extent_int(0);
    const int fv = feVel_->gradN().extent_int(1);
    const int fp = fePrs_->gradN().extent_int(1);
    const int p  = feVel_->gradN().extent_int(2);
    const int d  = feVel_->gradN().extent_int(3);
    // Initialize Jacobians.
    std::vector<std::vector<scalar_view>> J(numFields_);
    J[d].resize(numFields_);
    for (int i = 0; i < d; ++i) {
      J[i].resize(numFields_);
      for (int j = 0; j < d; ++j)
        J[i][j] = scalar_view("jac", c, fv, fv);
      J[i][d] = scalar_view("jac", c, fv, fp);
      J[d][i] = scalar_view("jac", c, fp, fv);
    }
    J[d][d] = scalar_view("jac", c, fp, fp);
    // Split coefficients into components.
    std::vector<scalar_view> U, Z;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Evaluate/interpolate finite element fields on cells.
    std::vector<scalar_view> valVel(d), gradVel(d);
    std::vector<std::vector<scalar_view>> dVel(d), dVelPhi(d);
    scalar_view nu("nu", c, p);
    scalar_view nuGradPhi("nuGradPhi", c, fv, p, d);
    scalar_view valCtrl("valCtrl", c, p);
    scalar_view alpha("alpha", c, p);
    scalar_view alphaPhi("alphaPhi", c, fv, p);
    scalar_view valVel0("valVel0", c, p, d);
    scalar_view valVelDotgradPhi("valVelDotgradPhi", c, fv, p);
    scalar_view nuNVel0("nuNVel0", c, fv, p);
    scalar_view posPrsN("posPrsN", c, fp, p);
    scalar_view negPrsN("negPrsN", c, fp, p);
    computeViscosity(nu);
    feCtrl_->evaluateValue(valCtrl, Z[0]);
    imp_->compute(alpha, valCtrl, feVel_->cubPts(), 0);
    for (int i = 0; i < d; ++i) {
      valVel[i]           = scalar_view("valVel", c, p);
      gradVel[i]          = scalar_view("gradVel", c, p, d);
      feVel_->evaluateValue(valVel[i], U[i]);
      feVel_->evaluateGradient(gradVel[i], U[i]);
      for (int j = 0; j < c; ++j) {
        for (int k = 0; k < p; ++k)
          valVel0(j,k,i) = (valVel[i])(j,k);
      }
      dVel[i].resize(d);
      dVelPhi[i].resize(d);
      for (int j = 0; j < d; ++j) {
        dVel[i][j]    = scalar_view("dVel", c, p);
        dVelPhi[i][j] = scalar_view("dVelPhi", c, fv, p);
        for (int k = 0; k < c; ++k) {
          for (int l = 0; l < p; ++l)
            (dVel[i][j])(k,l) = (gradVel[i])(k,l,j);
        }
      }
    }
    fst::scalarMultiplyDataField(nuNVel0, nu, feVel_->N());
    // Multiply velocity gradients with viscosity.
    fst::tensorMultiplyDataField(nuGradPhi, nu, feVel_->gradN());
    // Compute nonlinear terms in the Navier-Stokes equations.
    fst::dotMultiplyDataField(valVelDotgradPhi, valVel0, feVel_->gradN());
    rst::scale(valVelDotgradPhi, density_);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        fst::scalarMultiplyDataField(dVelPhi[i][j], dVel[i][j], feVel_->N());
        rst::scale(dVelPhi[i][j], density_);
      }
    }
    // Multiply velocity with alpha
    fst::scalarMultiplyDataField(alphaPhi, alpha, feVel_->N());
    // Negative pressure basis.
    Kokkos::deep_copy(posPrsN,fePrs_->NdetJ());
    Kokkos::deep_copy(negPrsN,fePrs_->N());
    rst::scale(negPrsN, static_cast<Real>(-1));
    /*** Evaluate weak form of the Jacobian. ***/
    fst::integrate(J[0][d],feVel_->NdetJ(),negPrsN,false);
    multiplyByRadius(nuNVel0,true,true);
    fst::integrate(J[0][0],nuNVel0,feVel_->NdetJ(),false);
    fst::integrate(J[d][0],fePrs_->N(),feVel_->NdetJ(),false);
    multiplyByRadius(nuGradPhi,true);
    multiplyByRadius(valVelDotgradPhi,true);
    multiplyByRadius(alphaPhi,true);
    multiplyByRadius(negPrsN,true);
    multiplyByRadius(posPrsN,true);
    for (int i = 0; i < d; ++i) {
      // Velocity equation.
      fst::integrate(J[i][i],nuGradPhi,feVel_->gradNdetJ(),true);
      fst::integrate(J[i][i],alphaPhi,feVel_->NdetJ(),true);
      if (!useStokes_) {
        fst::integrate(J[i][i],feVel_->NdetJ(),valVelDotgradPhi,true);
        for (int j = 0; j < d; ++j) {
          multiplyByRadius(dVelPhi[i][j],true);
          fst::integrate(J[i][j],feVel_->NdetJ(),dVelPhi[i][j],true);
        }
      }
      fst::integrate(J[i][d],feVel_->DNDdetJ(i),negPrsN,true);
      // Pressure equation.
      fst::integrate(J[d][i],posPrsN,feVel_->DND(i),true);
      rst::scale(J[d][i], static_cast<Real>(-1));
    }
    // APPLY DIRICHLET CONDITIONS
    applyDirichletToJacobian1(J);
    // Combine the jacobians.
    FieldUtils::combineFieldCoeff(jac, J, fieldInfo_, fieldInfo_);
  }


  void Jacobian_2(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c  = feVel_->gradN().extent_int(0);
    const int fv = feVel_->gradN().extent_int(1);
    const int fp = fePrs_->gradN().extent_int(1);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int p  = feVel_->gradN().extent_int(2);
    const int d  = feVel_->gradN().extent_int(3);
    // Initialize Jacobians.
    std::vector<std::vector<scalar_view>> J(numFields_);
    for (int i = 0; i < d; ++i) {
      J[i].resize(1);
      J[i][0] = scalar_view("jac", c, fv, fc);
    }
    J[d].resize(1);
    J[d][0] = scalar_view("jac", c, fp, fc);
    // Split coefficients into components.
    std::vector<scalar_view> U, Z;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Evaluate on FE basis.
    scalar_view valCtrl("valCtrl", c, p);
    scalar_view alpha("alpha", c, p);
    scalar_view alphaPhi("alphaPhi", c, fc, p);
    scalar_view valVel("valVel", c, p);
    scalar_view alphaU("alphaU", c, fc, p);
    feCtrl_->evaluateValue(valCtrl, Z[0]);
    imp_->compute(alpha, valCtrl, feVel_->cubPts(), 1);
    fst::scalarMultiplyDataField(alphaPhi, alpha, feCtrl_->N());
    for (int i = 0; i < d; ++i) {
      Kokkos::deep_copy(alphaU,static_cast<Real>(0));
      Kokkos::deep_copy(valVel,static_cast<Real>(0));
      feVel_->evaluateValue(valVel, U[i]);
      // Multiply velocity with alpha
      fst::scalarMultiplyDataField(alphaU, valVel, alphaPhi);
      // Integrate
      multiplyByRadius(alphaU,true);
      fst::integrate(J[i][0],feVel_->NdetJ(),alphaU,false);
    }
    // APPLY DIRICHLET CONDITIONS
    applyDirichletToJacobian2(J);
    // Combine the jacobians.
    FieldUtils::combineFieldCoeff<Real>(jac, J, fieldInfo_, fieldInfoCtrl_);
  }

  void Jacobian_3(std::vector<scalar_view> & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> Jacobian_3 is zero!");
  }

  void Hessian_11(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if (!useStokes_) {
      // Retrieve dimensions.
      const int c  = feVel_->gradN().extent_int(0);
      const int fv = feVel_->gradN().extent_int(1);
      const int fp = fePrs_->gradN().extent_int(1);
      const int p  = feVel_->gradN().extent_int(2);
      const int d  = feVel_->gradN().extent_int(3);
      // Initialize hessians.
      std::vector<std::vector<scalar_view>> H;
      H.resize(numFields_);
      H[d].resize(numFields_);
      for (int i = 0; i < d; ++i) {
        H[i].resize(numFields_);
        for (int j = 0; j < d; ++j)
          H[i][j] = scalar_view("hess", c, fv, fv);
        H[i][d] = scalar_view("hess", c, fv, fp);
        H[d][i] = scalar_view("hess", c, fp, fv);
      }
      H[d][d] = scalar_view("hess", c, fp, fp);
      // Split coefficients into components.
      std::vector<scalar_view> L;
      FieldUtils::splitFieldCoeff<Real>(L, l_coeff, fieldInfo_);
      // Apply Dirichlet conditions
      applyDirichletToMultiplier(L);
      // Evaluate/interpolate finite element fields on cells.
      std::vector<scalar_view> Lval(d), LPhi(d);
      for (int i = 0; i < d; ++i) {
        Lval[i] = scalar_view("Lval", c, p);
        LPhi[i] = scalar_view("LPhi", c, fv, p);
        feVel_->evaluateValue(Lval[i], L[i]);
        fst::scalarMultiplyDataField(LPhi[i], Lval[i], feVel_->N());
        rst::scale(LPhi[i], density_);
        multiplyByRadius(LPhi[i],true);
      }
      /*** Evaluate weak form of the Hessian. ***/
      for (int i = 0; i < d; ++i) {
        fst::integrate(H[i][i],LPhi[i],feVel_->DNDdetJ(i),false);
        fst::integrate(H[i][i],feVel_->DNDdetJ(i),LPhi[i],true);
        for (int j = 0; j < d; ++j) {
          fst::integrate(H[i][j],feVel_->DNDdetJ(j),LPhi[i],false);
          fst::integrate(H[i][j],LPhi[j],feVel_->DNDdetJ(i),true);
        }
      }
      // Combine the Hessians.
      FieldUtils::combineFieldCoeff<Real>(hess, H, fieldInfo_, fieldInfo_);
    }
    else {
      throw Exception::Zero(">>> Hessian_11 is zero!");
    }
  }

  void Hessian_12(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c  = feVel_->gradN().extent_int(0);
    const int fv = feVel_->gradN().extent_int(1);
    const int fp = fePrs_->gradN().extent_int(1);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int p  = feVel_->gradN().extent_int(2);
    const int d  = feVel_->gradN().extent_int(3);
    // Initialize hessians.
    std::vector<std::vector<scalar_view>> H(1);
    H[0].resize(numFields_);
    for (int i = 0; i < d; ++i)
      H[0][i] = scalar_view("hess", c, fc, fv);
    H[0][d] = scalar_view("hess", c, fc, fp);
    // Split coefficients into components.
    std::vector<scalar_view> L, Z;
    FieldUtils::splitFieldCoeff<Real>(L, l_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Apply Dirichlet conditions
    applyDirichletToMultiplier(L);
    // Evaluate/interpolate finite element fields on cells.
    scalar_view Z0("Z0", c, p);
    scalar_view alpha("alpha", c, p);
    scalar_view Lval("Lval", c, p);
    scalar_view alphaL("alphaL", c, fc, p);
    feCtrl_->evaluateValue(Z0, Z[0]);
    imp_->compute(alpha, Z0, feVel_->cubPts(), 1);
    /*** Evaluate weak form of the Hessian. ***/
    // Velocity equation.
    for (int i = 0; i < d; ++i) {
      // Multiply velocity with alpha
      Kokkos::deep_copy(Lval, static_cast<Real>(0));
      Kokkos::deep_copy(alphaL, static_cast<Real>(0));
      feVel_->evaluateValue(Lval, L[i]);
      for (int j = 0; j < c; ++j) {
        for (int k = 0; k < fc; ++k) {
          for (int l = 0; l < p; ++l) {
            alphaL(j,k,l) = alpha(j,l) * (feCtrl_->N())(j,k,l) * Lval(j,l);
          }
        }
      }
      multiplyByRadius(alphaL,true);
      fst::integrate(H[0][i],alphaL,feVel_->NdetJ(),false);
    }
    // Combine the Hessians.
    FieldUtils::combineFieldCoeff<Real>(hess, H, fieldInfoCtrl_, fieldInfo_);
  }

  void Hessian_13(std::vector<scalar_view> & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> Hessian_13 is zero!");
  }

  void Hessian_21(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c  = feVel_->gradN().extent_int(0);
    const int fv = feVel_->gradN().extent_int(1);
    const int fp = fePrs_->gradN().extent_int(1);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int p  = feVel_->gradN().extent_int(2);
    const int d  = feVel_->gradN().extent_int(3);
    // Initialize hessians.
    std::vector<std::vector<scalar_view>> H(numFields_);
    for (int i = 0; i < d; ++i) {
      H[i].resize(1);
      H[i][0] = scalar_view("hess", c, fv, fc);
    }
    H[d].resize(1);
    H[d][0] = scalar_view("hess", c, fp, fc);
    // Split coefficients into components.
    std::vector<scalar_view> L, Z;
    FieldUtils::splitFieldCoeff<Real>(L, l_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Apply Dirichlet conditions
    applyDirichletToMultiplier(L);
    // Evaluate/interpolate finite element fields on cells.
    scalar_view Z0("Z0", c, p);
    scalar_view alpha("alpha", c, p);
    scalar_view Lval("Lval", c, p);
    scalar_view alphaL("alphaL", c, fc, p);
    feCtrl_->evaluateValue(Z0, Z[0]);
    imp_->compute(alpha, Z0, feVel_->cubPts(), 1);
    /*** Evaluate weak form of the Hessian. ***/
    // Velocity equation.
    for (int i = 0; i < d; ++i) {
      // Multiply velocity with alpha
      Kokkos::deep_copy(Lval, static_cast<Real>(0));
      Kokkos::deep_copy(alphaL, static_cast<Real>(0));
      feVel_->evaluateValue(Lval, L[i]);
      for (int j = 0; j < c; ++j) {
        for (int k = 0; k < fc; ++k) {
          for (int l = 0; l < p; ++l) {
            alphaL(j,k,l) = alpha(j,l) * (feCtrl_->N())(j,k,l) * Lval(j,l);
          }
        }
      }
      multiplyByRadius(alphaL,true);
      fst::integrate(H[i][0],feVel_->NdetJ(),alphaL,false);
    }
    // Combine the Hessians.
    FieldUtils::combineFieldCoeff<Real>(hess, H, fieldInfo_, fieldInfoCtrl_);
  }

  void Hessian_22(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c  = feVel_->gradN().extent_int(0);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int p  = feVel_->gradN().extent_int(2);
    const int d  = feVel_->gradN().extent_int(3);
    // Initialize hessians.
    std::vector<std::vector<scalar_view>> H(1);
    H[0].resize(1);
    H[0][0] = scalar_view("hess", c, fc, fc);
    // Split coefficients into components.
    std::vector<scalar_view> L, U, Z;
    FieldUtils::splitFieldCoeff<Real>(L, l_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Apply Dirichlet conditions
    applyDirichletToMultiplier(L);
    // Evaluate/interpolate finite element fields on cells.
    std::vector<scalar_view> Lval(d), Uval(d);
    scalar_view Z0("Z0", c, p);
    scalar_view alpha("alpha", c, p);
    scalar_view alphaUL("alphaUL", c, fc, p);
    for (int i = 0; i < d; ++i) {
      Lval[i] = scalar_view("Lval", c, p);
      Uval[i] = scalar_view("Uval", c, p);
      feVel_->evaluateValue(Lval[i], L[i]);
      feVel_->evaluateValue(Uval[i], U[i]);
    }
    feCtrl_->evaluateValue(Z0, Z[0]);
    imp_->compute(alpha, Z0, feVel_->cubPts(), 2);
    // Multiply velocity with alpha
    Real dot(0);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < fc; ++j) {
        for (int k = 0; k < p; ++k) {
          dot = static_cast<Real>(0);
          for (int l = 0; l < d; ++l) {
            dot += (Uval[l])(i,k) * (Lval[l])(i,k);
          }
          alphaUL(i,j,k) = alpha(i,k) * (feCtrl_->N())(i,j,k) * dot;
        }
      }
    }
    /*** Evaluate weak form of the Hessian. ***/
    // Velocity equation.
    multiplyByRadius(alphaUL,true);
    fst::integrate(H[0][0],alphaUL,feCtrl_->NdetJ(),false);
    // Combine hessians.
    FieldUtils::combineFieldCoeff<Real>(hess, H, fieldInfoCtrl_, fieldInfoCtrl_);
  }

  void Hessian_23(std::vector<scalar_view> & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> Hessian_23 is zero!");
  }

  void Hessian_31(std::vector<scalar_view> & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> Hessian_31 is zero!");
  }

  void Hessian_32(std::vector<scalar_view> & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> Hessian_32 is zero!");
  }

  void Hessian_33(std::vector<std::vector<scalar_view>> & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> Hessian_33 is zero!");
  }

  void RieszMap_1(scalar_view & riesz) override {
    // Retrieve dimensions.
    const int c  = feVel_->gradN().extent_int(0);
    const int fv = feVel_->gradN().extent_int(1);
    const int fp = fePrs_->gradN().extent_int(1);
    const int p  = feVel_->gradN().extent_int(2);
    const int d  = feVel_->gradN().extent_int(3);
    // Initialize Riesz maps.
    std::vector<std::vector<scalar_view>> J(numFields_);
    J[d].resize(numFields_);
    for (int i = 0; i < d; ++i) {
      J[i].resize(numFields_);
      for (int j = 0; j < d; ++j)
        J[i][j] = scalar_view("riesz", c, fv, fv);
      J[i][d] = scalar_view("riesz", c, fv, fp);
      J[d][i] = scalar_view("riesz", c, fp, fv);
    }
    J[d][d] = scalar_view("riesz", c, fp, fp);
    // Evaluate/interpolate finite element fields on cells.
    scalar_view valVel("valVel", c, fv, p);
    scalar_view gradVel("gradVel", c, fv, p, d);
    scalar_view valPrs("valPrs", c, fp, p);
    Kokkos::deep_copy(valVel, feVel_->N());
    Kokkos::deep_copy(gradVel, feVel_->gradN());
    Kokkos::deep_copy(valPrs, fePrs_->N());
    /*** Evaluate weak form of the Jacobian. ***/
    multiplyByRadius(gradVel,true);
    multiplyByRadius(valVel,true);
    multiplyByRadius(valPrs,true);
    for (int i = 0; i < d; ++i) {
      fst::integrate(J[i][i],gradVel,feVel_->gradNdetJ(),false);
      fst::integrate(J[i][i],valVel,feVel_->NdetJ(),true);
    }
    fst::integrate(J[d][d],valPrs,fePrs_->NdetJ(),false);
    // Combine Riesz maps.
    FieldUtils::combineFieldCoeff<Real>(riesz, J, fieldInfo_, fieldInfo_);
  }

  void RieszMap_2(scalar_view & riesz) override {
    // Retrieve dimensions.
    const int c  = feVel_->gradN().extent_int(0);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int p  = feCtrl_->gradN().extent_int(2);
    // Initialize residuals.
    std::vector<std::vector<scalar_view>> J(1);
    J[0].resize(1);
    J[0][0] = scalar_view("riesz", c, fc, fc);
    // Evaluate on FE basis.
    scalar_view valPsi("valPsi", c, fc, p);
    Kokkos::deep_copy(valPsi, feCtrl_->N());
    /*** Evaluate weak form of the residual. ***/
    multiplyByRadius(valPsi,true);
    fst::integrate(J[0][0],valPsi,feCtrl_->NdetJ(),false);
    // Combine Riesz maps.
    FieldUtils::combineFieldCoeff<Real>(riesz, J, fieldInfoCtrl_, fieldInfoCtrl_);
  }

  std::vector<basis_ptr> getFields() override {
    return basisPtrs_;
  }

  std::vector<basis_ptr> getFields2() override {
    return basisPtrsCtrl_;
  }

  void setCellNodes(const scalar_view &volCellNodes,
                    const std::vector<std::vector<scalar_view>> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds) override {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    feVel_  = ROL::makePtr<fe_type>(volCellNodes_,basisPtrVel_,cellCub_);
    fePrs_  = ROL::makePtr<fe_type>(volCellNodes_,basisPtrPrs_,cellCub_);
    feCtrl_ = ROL::makePtr<fe_type>(volCellNodes_,basisPtrCtrl_,cellCub_,false);
    fvidx_  = feVel_->getBoundaryDofs();
    fpidx_  = fePrs_->getBoundaryDofs();
    // Construct control boundary FE
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      feVelBdry_.resize(numSideSets);
      for (int i = 0; i < numSideSets; ++i) {
        int numLocSides = bdryCellNodes[i].size();
        feVelBdry_[i].resize(numLocSides);
        for (int j = 0; j < numLocSides; ++j) {
          if (bdryCellNodes[i][j] != scalar_view()) {
            feVelBdry_[i][j] = ROL::makePtr<fe_type>(bdryCellNodes[i][j],basisPtrVel_,bdryCub_,j);
          }
        }
      }
    }
    computeDirichlet();
  }

  void setFieldPattern(const std::vector<std::vector<int>> &fieldPattern,
                       const std::vector<std::vector<int>> &fieldPattern2) override {
    fieldPattern_ = fieldPattern;
    fieldInfo_    = ROL::makePtr<FieldUtils::FieldInfo>(numFields_,numDofs_,numFieldDofs_,fieldPattern_);
    fieldPatternCtrl_ = fieldPattern2;
    fieldInfoCtrl_    = ROL::makePtr<FieldUtils::FieldInfo>(numFieldsCtrl_,numDofsCtrl_,numFieldDofsCtrl_,fieldPatternCtrl_);
    fieldHelper_ = ROL::makePtr<FieldHelper<Real,DeviceType>>(numFields_, numDofs_, numFieldDofs_, fieldPattern_);
  }

  const ROL::Ptr<fe_type> getVelocityFE(void) const {
    return feVel_;
  }

  const ROL::Ptr<fe_type> getPressureFE(void) const {
    return fePrs_;
  }

  const ROL::Ptr<fe_type> getControlFE(void) const {
    return feCtrl_;
  }

  const std::vector<ROL::Ptr<fe_type>> getVelocityBdryFE(const int sideset = -1) const {
    int side = sideset;
    if ( sideset < 0 || sideset > 4 ) {
      //side = (useDirichletControl_) ? 6 : 10;
      side = 0;
    }
    return feVelBdry_[side];
  }

  const std::vector<std::vector<int>> getBdryCellLocIds(const int sideset = -1) const {
    int side = sideset;
    if ( sideset < 0 || sideset > 4 ) {
      //side = (useDirichletControl_) ? 6 : 10;
      side = 0;
    }
    return bdryCellLocIds_[side];
  }

  const ROL::Ptr<FieldHelper<Real,DeviceType>> getFieldHelper(void) const {
    return fieldHelper_;
  }

  const ROL::Ptr<const FieldUtils::FieldInfo> getStateFieldInfo(void) const {
    return fieldInfo_;
  }

  const ROL::Ptr<const FieldUtils::FieldInfo> getControlFieldInfo(void) const {
    return fieldInfoCtrl_;
  }

  void printData(std::string tag,
                 const scalar_view u_coeff,
                 const scalar_view z_coeff = scalar_view(),
                 const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c = fePrs_->gradN().extent_int(0);
    const int p = fePrs_->gradN().extent_int(2);
    const int d = fePrs_->gradN().extent_int(3);
    // Evaluate/interpolate finite element fields on cells.
    scalar_view valCtrl("valCtrl", c, p);
    scalar_view alpha("alpha", c, p);
    feCtrl_->evaluateValue(valCtrl, z_coeff);
    imp_->compute(alpha, valCtrl, fePrs_->cubPts(), 0);
    // Print to permeability file
    std::stringstream namePer;
    namePer << "permeability" << tag << ".txt";
    std::ofstream filePer;
    filePer.open(namePer.str());
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        filePer << std::scientific << std::setprecision(15);
        for (int k = 0; k < d; ++k) filePer << std::setw(25) << (fePrs_->cubPts())(i,j,k);
        filePer << std::setw(25) << viscosity_/alpha(i,j);
        filePer << std::endl;
      }
    }
    filePer.close();
  }

private:

  void applyDirichletToResidual(std::vector<scalar_view> &R,
                          const std::vector<scalar_view> &U) const {
    const int d = feVel_->gradN().extent_int(3);
    const int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      // APPLY DIRICHLET CONDITIONS
      for (int i = 0; i < numSideSets; ++i) {
        if (i==0) /* in flow */ {
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            const int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              const int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < d; ++m) {
                  (R[m])(cidx,fvidx_[j][l]) = (U[m])(cidx,fvidx_[j][l]) - (bdryCellDofValues_[i][j])(k,fvidx_[j][l],m);
                }
              }
            }
          }
        }
        else if (i==1 && useNoSlip_) /* no slip */ {
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            const int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              const int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < d; ++m) {
                  (R[m])(cidx,fvidx_[j][l]) = (U[m])(cidx,fvidx_[j][l]);
                }
              }
            }
          }
        }
        else if (i==2 && usePresOut_) /* out flow */ {
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            const int numBdryDofs = fpidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              const int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                (R[d])(cidx,fpidx_[j][l]) = (U[d])(cidx,fpidx_[j][l]) - Patm_;
              }
            }
          }
        }
        else if (i==3) /* no normal */ { 
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            const int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              const int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                (R[0])(cidx,fvidx_[j][l]) = (U[0])(cidx,fvidx_[j][l]);
              }
            }
          }
        }
      }
    }
  }

  void applyDirichletToJacobian1(std::vector<std::vector<scalar_view>> &J) const {
    const int fv = feVel_->gradN().extent_int(1);
    const int fp = fePrs_->gradN().extent_int(1);
    const int d  = feVel_->gradN().extent_int(3);
    const int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==0) /* in flow */ {
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            const int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              const int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < fv; ++m) {
                  for (int n=0; n < d; ++n) {
                    for (int p=0; p < d; ++p) {
                      (J[n][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                    }
                    (J[n][n])(cidx,fvidx_[j][l],fvidx_[j][l]) = static_cast<Real>(1);
                  }
                }
                for (int m=0; m < fp; ++m) {
                  for (int n=0; n < d; ++n) {
                    (J[n][d])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
              }
            }
          }
        }
        else if (i==1 && useNoSlip_) /* no slip */ {
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            const int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              const int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < fv; ++m) {
                  for (int n=0; n < d; ++n) {
                    for (int p=0; p < d; ++p) {
                      (J[n][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                    }
                    (J[n][n])(cidx,fvidx_[j][l],fvidx_[j][l]) = static_cast<Real>(1);
                  }
                }
                for (int m=0; m < fp; ++m) {
                  for (int n=0; n < d; ++n) {
                    (J[n][d])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
              }
            }
          }
        }
        else if (i==2 && usePresOut_) /* out flow */ {
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            const int numBdryDofs = fpidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              const int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m = 0; m < fv; ++m) {
                  for (int p = 0; p < d; ++p) {
                    (J[d][p])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
                for (int m = 0; m < fp; ++m) {
                  (J[d][d])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
                }
                (J[d][d])(cidx,fpidx_[j][l],fpidx_[j][l]) = static_cast<Real>(1);
              }
            }
          }
        }
        else if ( i==3 ) /* no normal */ {
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            const int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              const int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < fv; ++m) {
                  for (int p=0; p < d; ++p) {
                    (J[0][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                  (J[0][0])(cidx,fvidx_[j][l],fvidx_[j][l]) = static_cast<Real>(1);
                }
                for (int m=0; m < fp; ++m) {
                  (J[0][d])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                }
              }
            }
          }
        }
      }
    }
  }

  void applyDirichletToJacobian2(std::vector<std::vector<scalar_view>> &J) const {
    const int fc = feCtrl_->gradN().extent_int(1);
    const int d  = feVel_->gradN().extent_int(3);
    const int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==0) /* in flow */ {
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            const int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              const int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int n=0; n < d; ++n) {
                  for (int m=0; m < fc; ++m) {
                    (J[n][0])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
              }
            }
          }
        }
        else if (i==1 && useNoSlip_) /* no slip */ {
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            const int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              const int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int n=0; n < d; ++n) {
                  for (int m=0; m < fc; ++m) {
                    (J[n][0])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
              }
            }
          }
        }
        else if (i==2 && usePresOut_) /* out flow */ {}
        else if ( i==3 ) /* no normal */ {
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            const int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              const int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < fc; ++m) {
                  (J[0][0])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                }
              }
            }
          }
        }
      }
    }
  }

  void applyDirichletToMultiplier(std::vector<scalar_view> &L) const {
    const int d = feVel_->gradN().extent_int(3);
    const int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==0) /* in flow */ {
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            const int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              const int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < d; ++m) {
                  (L[m])(cidx,fvidx_[j][l]) = static_cast<Real>(0);
                }
              }
            }
          }
        }
        else if (i==1 && useNoSlip_) /* no slip */ {
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            const int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              const int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < d; ++m) {
                  (L[m])(cidx,fvidx_[j][l]) = static_cast<Real>(0);
                }
              }
            }
          }
        }
        else if (i==2 && usePresOut_) /* out flow */ {
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            const int numBdryDofs = fpidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              const int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                (L[d])(cidx,fpidx_[j][l]) = static_cast<Real>(0);
              }
            }
          }
        }
        else if ( i==3 ) /* no normal */ {
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            const int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              const int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                (L[0])(cidx,fvidx_[j][l]) = static_cast<Real>(0);
              }
            }
          }
        }
      }
    }
  }

}; // PDE_NavierStokes

#endif
