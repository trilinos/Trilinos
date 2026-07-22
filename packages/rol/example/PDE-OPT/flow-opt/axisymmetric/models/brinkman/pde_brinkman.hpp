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

#ifndef PDE_BRINKMAN_HPP
#define PDE_BRINKMAN_HPP

#include "../../../../TOOLS/pde.hpp"
#include "../../../../TOOLS/fe.hpp"
#include "../../../../TOOLS/fieldhelper.hpp"
#include "impermeability.hpp"

#include "../../../../TOOLS/Intrepid_HGRAD_C0_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid_HGRAD_TRI_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real>
class PDE_NavierStokes : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>> basisPtrVel_;
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>> basisPtrPrs_;
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>> basisPtrCtrl_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrs_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrsCtrl_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real>> cellCub_;
  ROL::Ptr<Intrepid::Cubature<Real>> bdryCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real>> volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real>> feVel_;
  ROL::Ptr<FE<Real>> fePrs_;
  ROL::Ptr<FE<Real>> feCtrl_;
  std::vector<std::vector<ROL::Ptr<FE<Real>>>> feVelBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fvidx_;
  std::vector<std::vector<int>> fpidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellDofValues_;
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

  ROL::Ptr<Impermeability<Real>> imp_;
  ROL::Ptr<FieldHelper<Real>> fieldHelper_;
  ROL::Ptr<FieldUtils::FieldInfo> fieldInfo_, fieldInfoCtrl_;

  // Extract velocity coefficients on boundary.
  ROL::Ptr<Intrepid::FieldContainer<Real>> getBoundaryCoeff(
      const Intrepid::FieldContainer<Real> & cell_coeff,
      int sideSet, int cell) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][cell];
    const int numCellsSide = bdryCellLocId.size();
    const int f = basisPtrVel_->getCardinality();
    
    ROL::Ptr<Intrepid::FieldContainer<Real >> bdry_coeff = 
      ROL::makePtr<Intrepid::FieldContainer<Real >>(numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        (*bdry_coeff)(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

  Real DirichletFunc(const std::vector<Real> & coords, int sideset, int locSideId, int dir) const {
    const std::vector<Real> param = PDE<Real>::getParameter();
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
        bdryCellDofValues_[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, d);
        ROL::Ptr<Intrepid::FieldContainer<Real>> coords =
          ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, d);
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
              (*bdryCellDofValues_[i][j])(k, l, m) = DirichletFunc(dofpoint, i, j, m);
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

  void computeViscosity(ROL::Ptr<Intrepid::FieldContainer<Real>> &nu) const {
    int c = feVel_->gradN()->dimension(0);
    int p = feVel_->gradN()->dimension(2);
    int d = feVel_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*feVel_->cubPts())(i,j,k);
        }
        // Compute spatially dependent viscosity
        (*nu)(i,j) = viscosityFunc(pt);
      }
    }
  }

  void multiplyByRadius(ROL::Ptr<Intrepid::FieldContainer<Real>> &input, bool isField, bool useReciprocal = false) const {
    int c = feVel_->gradN()->dimension(0);
    int p = feVel_->gradN()->dimension(2);
    Intrepid::FieldContainer<Real> r(c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        if (!useReciprocal) r(i,j) = (*feVel_->cubPts())(i,j,0);
        else                r(i,j) = 1/(*feVel_->cubPts())(i,j,0);
      }
    }
    Intrepid::FieldContainer<Real> in = *input;
    if (!isField)
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*input,r,in);
    else
      Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*input,r,in);
  }

public:
  PDE_NavierStokes(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int cubDegree     = parlist.sublist("Problem").get("Cubature Degree",4);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2);
    int probDim       = parlist.sublist("Problem").get("Problem Dimension",2);
    int basisDegCtrl  = parlist.sublist("Problem").get("Density Basis Degree",1);
    std::string elemType = parlist.sublist("Problem").get("Element Type","QUAD");
    TEUCHOS_TEST_FOR_EXCEPTION(probDim!=2, std::invalid_argument,
      ">>> rol/example/PDE-OPT/flow-opt/axisymmetric/models/brinkman/pde_brinkman.hpp: Problem dimension is not 2!");
    if (elemType == "TRI") {
      basisPtrVel_ = ROL::makePtr<Intrepid::Basis_HGRAD_TRI_C2_FEM<Real, Intrepid::FieldContainer<Real>>>();
      basisPtrPrs_ = ROL::makePtr<Intrepid::Basis_HGRAD_TRI_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      if (basisDegCtrl == 1)
        basisPtrCtrl_ = ROL::makePtr<Intrepid::Basis_HGRAD_TRI_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      else
        basisPtrCtrl_ = ROL::makePtr<Intrepid::Basis_HGRAD_C0_FEM<Real, Intrepid::FieldContainer<Real>>>();
    }
    else {
      basisPtrVel_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real,Intrepid::FieldContainer<Real>>>();
      basisPtrPrs_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real,Intrepid::FieldContainer<Real>>>();
      if (basisDegCtrl == 1)
        basisPtrCtrl_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real,Intrepid::FieldContainer<Real>>>();
      else
        basisPtrCtrl_ = ROL::makePtr<Intrepid::Basis_HGRAD_C0_FEM<Real,Intrepid::FieldContainer<Real>>>();
    }
    basisPtrs_.clear(); basisPtrsCtrl_.clear();
    for (int i=0; i<probDim; ++i) basisPtrs_.push_back(basisPtrVel_); // Velocity component
    basisPtrs_.push_back(basisPtrPrs_);                               // Pressure component
    basisPtrsCtrl_.push_back(basisPtrCtrl_);                          // Control component
    // Quadrature rules.
    shards::CellTopology cellType = basisPtrs_[0]->getBaseCellTopology(); // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                    // create cubature factory
    cellCub_ = cubFactory.create(cellType, cubDegree);                    // create default cubature

    shards::CellTopology bdryCellType = cellType.getCellTopologyData(probDim-1, 0);
    bdryCub_ = cubFactory.create(bdryCellType, bdryCubDegree);

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
    imp_        = ROL::makePtr<Impermeability<Real>>(parlist);

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

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real>> & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int p  = feVel_->gradN()->dimension(2);
    const int d  = feVel_->gradN()->dimension(3);
    // Initialize residuals.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> R;
    R.resize(numFields_);
    for (int i = 0; i < d; ++i) {
      R[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    R[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    // Split coefficients into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U, Z;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Evaluate/interpolate finite element fields on cells.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> nuGradVel(d), valVel(d), gradVel(d), valVelDotgradVel(d), alphaVel(d);
    ROL::Ptr<Intrepid::FieldContainer<Real>> nu, valVel0, valPres, divVel, valCtrl, alpha, nuVel0;
    nu      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    valVel0 = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    valPres = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    divVel  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    valCtrl = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alpha   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    nuVel0  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    computeViscosity(nu);
    feCtrl_->evaluateValue(valCtrl, Z[0]);
    imp_->compute(alpha, valCtrl, feVel_->cubPts(), 0);
    for (int i = 0; i < d; ++i) {
      nuGradVel[i]        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      valVel[i]           = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      gradVel[i]          = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      valVelDotgradVel[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      alphaVel[i]         = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      // Evaluate on FE basis
      feVel_->evaluateValue(valVel[i], U[i]);
      feVel_->evaluateGradient(gradVel[i], U[i]);
      // Multiply velocity gradients with viscosity.
      Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(*nuGradVel[i], *nu, *gradVel[i]);
      // Multiply velocity with alpha
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*alphaVel[i], *alpha, *valVel[i]);
      // Assemble the velocity vector and its divergence.
      for (int j = 0; j < c; ++j) {
        for (int k = 0; k < p; ++k) {
          (*valVel0)(j,k,i) = (*valVel[i])(j,k);
          (*divVel)(j,k)   += (*gradVel[i])(j,k,i);
        }
      }
    }
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*nuVel0, *nu, *valVel[0]);
    // Compute nonlinear terms in the Navier-Stokes equations.
    for (int i = 0; i < d; ++i) {
      Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(*valVelDotgradVel[i], *valVel0, *gradVel[i]);
      Intrepid::RealSpaceTools<Real>::scale(*valVelDotgradVel[i], density_);
    }
    // Negative pressure
    fePrs_->evaluateValue(valPres, U[d]);
    Intrepid::RealSpaceTools<Real>::scale(*valPres, static_cast<Real>(-1));
    /*** Evaluate weak form of the residual. ***/
    // Integrate velocity equation.
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[0],
                                                  *valPres,                 // p
                                                  *(feVel_->NdetJ()),       // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    multiplyByRadius(nuVel0,false,true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[0],
                                                  *nuVel0,                  // nu UX
                                                  *(feVel_->NdetJ()),       // Phi
                                                  Intrepid::COMP_CPP,
                                                  true);
    multiplyByRadius(valPres,false);
    for (int i = 0; i < d; ++i) {
      multiplyByRadius(nuGradVel[i],false);
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *nuGradVel[i],            // nu gradUX
                                                    *(feVel_->gradNdetJ()),   // gradPhi
                                                    Intrepid::COMP_CPP,
                                                    true);
      if (!useStokes_) {
        multiplyByRadius(valVelDotgradVel[i],false);
        Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                      *valVelDotgradVel[i],   // (U . gradUX)
                                                      *(feVel_->NdetJ()),     // Phi
                                                      Intrepid::COMP_CPP,
                                                      true);
      }
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *valPres,                 // p
                                                    *(feVel_->DNDdetJ(i)),    // dPhi/dx
                                                    Intrepid::COMP_CPP,
                                                    true);
      multiplyByRadius(alphaVel[i],false);
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *alphaVel[i],             // alpha UX
                                                    *(feVel_->NdetJ()),       // Phi
                                                    Intrepid::COMP_CPP,
                                                    true);
    }
    // Integrate pressure equation.
    multiplyByRadius(divVel,false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[d],
                                                  *divVel,                  // divU
                                                  *(fePrs_->NdetJ()),       // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[d],
                                                  *valVel[0],               // UX
                                                  *(fePrs_->NdetJ()),       // Phi
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::RealSpaceTools<Real>::scale(*R[d], static_cast<Real>(-1));
    // Boundary conditions.
    applyDirichletToResidual(R,U);
    // Combine the residuals.
    FieldUtils::combineFieldCoeff<Real>(res, R, fieldInfo_);
  }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int p  = feVel_->gradN()->dimension(2);
    const int d  = feVel_->gradN()->dimension(3);
    // Initialize Jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J(numFields_);
    J[d].resize(numFields_);
    for (int i = 0; i < d; ++i) {
      J[i].resize(numFields_);
      for (int j = 0; j < d; ++j) {
        J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fv);
      }
      J[i][d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fp);
      J[d][i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fv);
    }
    J[d][d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fp);
    // Split coefficients into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U, Z;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Evaluate/interpolate finite element fields on cells.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> valVel(d), gradVel(d);
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> dVel(d), dVelPhi(d);
    ROL::Ptr<Intrepid::FieldContainer<Real>> nu, nuGradPhi, alpha, alphaPhi, valCtrl;
    ROL::Ptr<Intrepid::FieldContainer<Real>> valVel0, valVelDotgradPhi, nuNVel0, posPrsN, negPrsN;
    nu                = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    nuGradPhi         = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p, d);
    valCtrl           = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alpha             = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alphaPhi          = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p);
    valVel0           = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    valVelDotgradPhi  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p);
    nuNVel0           = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p);
    posPrsN           = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, p);
    negPrsN           = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, p);
    computeViscosity(nu);
    feCtrl_->evaluateValue(valCtrl, Z[0]);
    imp_->compute(alpha, valCtrl, feVel_->cubPts(), 0);
    for (int i = 0; i < d; ++i) {
      valVel[i]           = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      gradVel[i]          = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      feVel_->evaluateValue(valVel[i], U[i]);
      feVel_->evaluateGradient(gradVel[i], U[i]);
      for (int j = 0; j < c; ++j) {
        for (int k = 0; k < p; ++k) {
          (*valVel0)(j,k,i) = (*valVel[i])(j,k);
        }
      }
      dVel[i].resize(d);
      dVelPhi[i].resize(d);
      for (int j = 0; j < d; ++j) {
        dVel[i][j]    = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
        dVelPhi[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p);
        for (int k = 0; k < c; ++k) {
          for (int l = 0; l < p; ++l) {
            (*dVel[i][j])(k,l) = (*gradVel[i])(k,l,j);
          }
        }
      }
    }
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*nuNVel0, *nu, *(feVel_->N()));
    // Multiply velocity gradients with viscosity.
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(*nuGradPhi, *nu, *(feVel_->gradN()));
    // Compute nonlinear terms in the Navier-Stokes equations.
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*valVelDotgradPhi, *valVel0, *(feVel_->gradN()));
    Intrepid::RealSpaceTools<Real>::scale(*valVelDotgradPhi, density_);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*dVelPhi[i][j], *dVel[i][j], *(feVel_->N()));
        Intrepid::RealSpaceTools<Real>::scale(*dVelPhi[i][j], density_);
      }
    }
    // Multiply velocity with alpha
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*alphaPhi, *alpha, *(feVel_->N()));
    // Negative pressure basis.
    *posPrsN = *(fePrs_->NdetJ());
    *negPrsN = *(fePrs_->N());
    Intrepid::RealSpaceTools<Real>::scale(*negPrsN, static_cast<Real>(-1));
    /*** Evaluate weak form of the Jacobian. ***/
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][d],
                                                  *(feVel_->NdetJ()),       // Phi
                                                  *negPrsN,                 // -Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    multiplyByRadius(nuNVel0,true,true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0],
                                                  *nuNVel0,                 // nu Phi
                                                  *(feVel_->NdetJ()),       // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[d][0],
                                                  *(fePrs_->N()),           // -Phi
                                                  *(feVel_->NdetJ()),       // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    multiplyByRadius(nuGradPhi,true);
    multiplyByRadius(valVelDotgradPhi,true);
    multiplyByRadius(alphaPhi,true);
    multiplyByRadius(negPrsN,true);
    multiplyByRadius(posPrsN,true);
    for (int i = 0; i < d; ++i) {
      // Velocity equation.
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][i],
                                                    *nuGradPhi,             // nu gradPhi
                                                    *(feVel_->gradNdetJ()), // gradPhi
                                                    Intrepid::COMP_CPP,
                                                    true);
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][i],
                                                    *alphaPhi,              // alpha Phi
                                                    *(feVel_->NdetJ()),     // Phi
                                                    Intrepid::COMP_CPP,
                                                    true);
      if (!useStokes_) {
        Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][i],
                                                      *(feVel_->NdetJ()),   // Phi
                                                      *valVelDotgradPhi,    // (U.gradPhiX)
                                                      Intrepid::COMP_CPP,
                                                      true);
        for (int j = 0; j < d; ++j) {
          multiplyByRadius(dVelPhi[i][j],true);
          Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][j],
                                                        *(feVel_->NdetJ()), // Phi
                                                        *dVelPhi[i][j],     // (Phi.gradU)
                                                        Intrepid::COMP_CPP,
                                                        true);
        }
      }
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][d],
                                                    *(feVel_->DNDdetJ(i)),  // dPhi/dx
                                                    *negPrsN,               // -Phi
                                                    Intrepid::COMP_CPP,
                                                    true);
      // Pressure equation.
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[d][i],
                                                    *posPrsN,               // Phi
                                                    *(feVel_->DND(i)),      // dPhi/dx
                                                    Intrepid::COMP_CPP,
                                                    true);
      Intrepid::RealSpaceTools<Real>::scale(*J[d][i], static_cast<Real>(-1));
    }
    // APPLY DIRICHLET CONDITIONS
    applyDirichletToJacobian1(J);
    // Combine the jacobians.
    FieldUtils::combineFieldCoeff<Real>(jac, J, fieldInfo_, fieldInfo_);
  }


  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int fc = feCtrl_->gradN()->dimension(1);
    const int p  = feVel_->gradN()->dimension(2);
    const int d  = feVel_->gradN()->dimension(3);
    // Initialize Jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J(numFields_);
    for (int i = 0; i < d; ++i) {
      J[i].resize(1,ROL::nullPtr);
      J[i][0] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fc);
    }
    J[d].resize(1,ROL::nullPtr);
    J[d][0] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fc);
    // Split coefficients into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U, Z;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Evaluate on FE basis.
    ROL::Ptr<Intrepid::FieldContainer<Real>> valCtrl, alpha, alphaPhi, valVel, alphaU;
    valCtrl  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alpha    = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alphaPhi = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc, p);
    valVel   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alphaU   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc, p);
    feCtrl_->evaluateValue(valCtrl, Z[0]);
    imp_->compute(alpha, valCtrl, feVel_->cubPts(), 1);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*alphaPhi, *alpha, *feCtrl_->N());
    for (int i = 0; i < d; ++i) {
      alphaU->initialize();
      valVel->initialize();
      feVel_->evaluateValue(valVel, U[i]);
      // Multiply velocity with alpha
      Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*alphaU, *valVel, *alphaPhi);
      // Integrate
      multiplyByRadius(alphaU,true);
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][0],
                                                    *feVel_->NdetJ(),  // Phi
                                                    *alphaU,           // alpha' Phi U
                                                    Intrepid::COMP_CPP,
                                                    false);
    }
    // APPLY DIRICHLET CONDITIONS
    applyDirichletToJacobian2(J);
    // Combine the jacobians.
    FieldUtils::combineFieldCoeff<Real>(jac, J, fieldInfo_, fieldInfoCtrl_);
  }

  void Jacobian_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> Jacobian_3 is zero!");
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if (!useStokes_) {
      // Retrieve dimensions.
      const int c  = feVel_->gradN()->dimension(0);
      const int fv = feVel_->gradN()->dimension(1);
      const int fp = fePrs_->gradN()->dimension(1);
      const int p  = feVel_->gradN()->dimension(2);
      const int d  = feVel_->gradN()->dimension(3);
      // Initialize hessians.
      std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> H;
      H.resize(numFields_);
      H[d].resize(numFields_);
      for (int i = 0; i < d; ++i) {
        H[i].resize(numFields_);
        for (int j = 0; j < d; ++j) {
          H[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fv);
        }
        H[i][d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fp);
        H[d][i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fv);
      }
      H[d][d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fp);
      // Split coefficients into components.
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> L;
      FieldUtils::splitFieldCoeff<Real>(L, l_coeff, fieldInfo_);
      // Apply Dirichlet conditions
      applyDirichletToMultiplier(L);
      // Evaluate/interpolate finite element fields on cells.
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> Lval(d), LPhi(d);
      for (int i = 0; i < d; ++i) {
        Lval[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
        LPhi[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p);
        feVel_->evaluateValue(Lval[i], L[i]);
        Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*LPhi[i], *Lval[i], *(feVel_->N()));
        Intrepid::RealSpaceTools<Real>::scale(*LPhi[i], density_);
        multiplyByRadius(LPhi[i],true);
      }
      /*** Evaluate weak form of the Hessian. ***/
      for (int i = 0; i < d; ++i) {
        Intrepid::FunctionSpaceTools::integrate<Real>(*H[i][i],
                                                      *LPhi[i],                // L Phi
                                                      *(feVel_->DNDdetJ(i)),   // dPhi/dx
                                                      Intrepid::COMP_CPP,
                                                      false);
        Intrepid::FunctionSpaceTools::integrate<Real>(*H[i][i],
                                                      *(feVel_->DNDdetJ(i)),   // dPhi/dx
                                                      *LPhi[i],                // L Phi
                                                      Intrepid::COMP_CPP,
                                                      true);
        for (int j = 0; j < d; ++j) {
          Intrepid::FunctionSpaceTools::integrate<Real>(*H[i][j],
                                                        *(feVel_->DNDdetJ(j)),   // dPhi/dy
                                                        *LPhi[i],                // L Phi
                                                        Intrepid::COMP_CPP,
                                                        false);
          Intrepid::FunctionSpaceTools::integrate<Real>(*H[i][j],
                                                        *LPhi[j],                // L Phi
                                                        *(feVel_->DNDdetJ(i)),   // dPhi/dx
                                                        Intrepid::COMP_CPP,
                                                        true);
        }
      }
      // Combine the Hessians.
      FieldUtils::combineFieldCoeff<Real>(hess, H, fieldInfo_, fieldInfo_);
    }
    else {
      throw Exception::Zero(">>> Hessian_11 is zero!");
    }
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int fc = feCtrl_->gradN()->dimension(1);
    const int p  = feVel_->gradN()->dimension(2);
    const int d  = feVel_->gradN()->dimension(3);
    // Initialize hessians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> H(1);
    H[0].resize(numFields_);
    for (int i = 0; i < d; ++i) {
      H[0][i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc, fv);
    }
    H[0][d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc, fp);
    // Split coefficients into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> L, Z;
    FieldUtils::splitFieldCoeff<Real>(L, l_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Apply Dirichlet conditions
    applyDirichletToMultiplier(L);
    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> Z0, alpha, Lval, alphaL;
    Z0     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alpha  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Lval   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alphaL = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc, p);
    feCtrl_->evaluateValue(Z0, Z[0]);
    imp_->compute(alpha, Z0, feVel_->cubPts(), 1);
    /*** Evaluate weak form of the Hessian. ***/
    // Velocity equation.
    for (int i = 0; i < d; ++i) {
      // Multiply velocity with alpha
      Lval->initialize(); alphaL->initialize();
      feVel_->evaluateValue(Lval, L[i]);
      for (int j = 0; j < c; ++j) {
        for (int k = 0; k < fc; ++k) {
          for (int l = 0; l < p; ++l) {
            (*alphaL)(j,k,l) = (*alpha)(j,l) * (*feCtrl_->N())(j,k,l) * (*Lval)(j,l);
          }
        }
      }
      multiplyByRadius(alphaL,true);
      Intrepid::FunctionSpaceTools::integrate<Real>(*H[0][i],
                                                    *alphaL,               // alpha' L Phi
                                                    *(feVel_->NdetJ()),    // Phi
                                                    Intrepid::COMP_CPP,
                                                    false);
    }
    // Combine the Hessians.
    FieldUtils::combineFieldCoeff<Real>(hess, H, fieldInfoCtrl_, fieldInfo_);
  }

  void Hessian_13(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> Hessian_13 is zero!");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int fc = feCtrl_->gradN()->dimension(1);
    const int p  = feVel_->gradN()->dimension(2);
    const int d  = feVel_->gradN()->dimension(3);
    // Initialize hessians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> H(numFields_);
    for (int i = 0; i < d; ++i) {
      H[i].resize(1,ROL::nullPtr);
      H[i][0] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fc);
    }
    H[d].resize(1,ROL::nullPtr);
    H[d][0] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fc);
    // Split coefficients into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> L, Z;
    FieldUtils::splitFieldCoeff<Real>(L, l_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Apply Dirichlet conditions
    applyDirichletToMultiplier(L);
    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> Z0, alpha, Lval, alphaL;
    Z0     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alpha  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Lval   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alphaL = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc, p);
    feCtrl_->evaluateValue(Z0, Z[0]);
    imp_->compute(alpha, Z0, feVel_->cubPts(), 1);
    /*** Evaluate weak form of the Hessian. ***/
    // Velocity equation.
    for (int i = 0; i < d; ++i) {
      // Multiply velocity with alpha
      Lval->initialize(); alphaL->initialize();
      feVel_->evaluateValue(Lval, L[i]);
      for (int j = 0; j < c; ++j) {
        for (int k = 0; k < fc; ++k) {
          for (int l = 0; l < p; ++l) {
            (*alphaL)(j,k,l) = (*alpha)(j,l) * (*feCtrl_->N())(j,k,l) * (*Lval)(j,l);
          }
        }
      }
      multiplyByRadius(alphaL,true);
      Intrepid::FunctionSpaceTools::integrate<Real>(*H[i][0],
                                                    *(feVel_->NdetJ()),    // Phi
                                                    *alphaL,               // alpha' L Phi
                                                    Intrepid::COMP_CPP,
                                                    false);
    }
    // Combine the Hessians.
    FieldUtils::combineFieldCoeff<Real>(hess, H, fieldInfo_, fieldInfoCtrl_);
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c  = feVel_->gradN()->dimension(0);
    const int fc = feCtrl_->gradN()->dimension(1);
    const int p  = feVel_->gradN()->dimension(2);
    const int d  = feVel_->gradN()->dimension(3);
    // Initialize hessians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> H(1);
    H[0].resize(1,ROL::nullPtr);
    H[0][0] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc, fc);
    // Split coefficients into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> L, U, Z;
    FieldUtils::splitFieldCoeff<Real>(L, l_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Apply Dirichlet conditions
    applyDirichletToMultiplier(L);
    // Evaluate/interpolate finite element fields on cells.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> Lval(d), Uval(d);
    ROL::Ptr<Intrepid::FieldContainer<Real>> Z0, alpha, alphaUL;
    Z0      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alpha   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alphaUL = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc, p);
    for (int i = 0; i < d; ++i) {
      Lval[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      Uval[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
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
            dot += (*Uval[l])(i,k) * (*Lval[l])(i,k);
          }
          (*alphaUL)(i,j,k) = (*alpha)(i,k) * (*feCtrl_->N())(i,j,k) * dot;
        }
      }
    }
    /*** Evaluate weak form of the Hessian. ***/
    // Velocity equation.
    multiplyByRadius(alphaUL,true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*H[0][0],
                                                  *alphaUL,              // alpha'' U L Phi
                                                  *(feCtrl_->NdetJ()),   // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    // Combine hessians.
    FieldUtils::combineFieldCoeff<Real>(hess, H, fieldInfoCtrl_, fieldInfoCtrl_);
  }

  void Hessian_23(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> Hessian_23 is zero!");
  }

  void Hessian_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> Hessian_31 is zero!");
  }

  void Hessian_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> Hessian_32 is zero!");
  }

  void Hessian_33(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> Hessian_33 is zero!");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) override {
    // Retrieve dimensions.
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int p  = feVel_->gradN()->dimension(2);
    const int d  = feVel_->gradN()->dimension(3);
    // Initialize Riesz maps.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J(numFields_);
    J[d].resize(numFields_);
    for (int i = 0; i < d; ++i) {
      J[i].resize(numFields_);
      for (int j = 0; j < d; ++j) {
        J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fv);
      }
      J[i][d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fp);
      J[d][i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fv);
    }
    J[d][d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fp);
    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> valVel, gradVel, valPrs;
    valVel  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p);
    gradVel = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p, d);
    valPrs  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, p);
    *valVel  = *feVel_->N();
    *gradVel = *feVel_->gradN();
    *valPrs  = *fePrs_->N();
    /*** Evaluate weak form of the Jacobian. ***/
    multiplyByRadius(gradVel,true);
    multiplyByRadius(valVel,true);
    multiplyByRadius(valPrs,true);
    for (int i = 0; i < d; ++i) {
      Intrepid::FunctionSpaceTools::integrate<Real>(*(J[i][i]),
                                                    *gradVel,             // alpha dPhi/dx
                                                    *feVel_->gradNdetJ(), // dPhi/dx
                                                    Intrepid::COMP_CPP,
                                                    false);
      Intrepid::FunctionSpaceTools::integrate<Real>(*(J[i][i]),
                                                    *valVel,              // alpha dPhi/dx
                                                    *feVel_->NdetJ(),     // dPhi/dx
                                                    Intrepid::COMP_CPP,
                                                    true);
    }
    Intrepid::FunctionSpaceTools::integrate<Real>(*(J[d][d]),
                                                  *valPrs,              // alpha dPhi/dx
                                                  *fePrs_->NdetJ(),     // dPhi/dx
                                                  Intrepid::COMP_CPP,
                                                  false);
    // Combine Riesz maps.
    FieldUtils::combineFieldCoeff<Real>(riesz, J, fieldInfo_, fieldInfo_);
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) override {
    // Retrieve dimensions.
    const int c  = feVel_->gradN()->dimension(0);
    const int fc = feCtrl_->gradN()->dimension(1);
    const int p  = feCtrl_->gradN()->dimension(2);
    // Initialize residuals.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J(1);
    J[0].resize(1,ROL::nullPtr);
    J[0][0] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc, fc);
    // Evaluate on FE basis.
    ROL::Ptr<Intrepid::FieldContainer<Real>> valPsi;
    valPsi  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc, p);
    *valPsi = *feCtrl_->N();
    /*** Evaluate weak form of the residual. ***/
    multiplyByRadius(valPsi,true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*(J[0][0]),
                                                  *valPsi,           // Psi
                                                  *feCtrl_->NdetJ(), // Psi
                                                  Intrepid::COMP_CPP,
                                                  false);
    // Combine Riesz maps.
    FieldUtils::combineFieldCoeff<Real>(riesz, J, fieldInfoCtrl_, fieldInfoCtrl_);
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields() override {
    return basisPtrs_;
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields2() override {
    return basisPtrsCtrl_;
  }

  void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real>> &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds) override {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    feVel_  = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrVel_,cellCub_);
    fePrs_  = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrPrs_,cellCub_);
    feCtrl_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrCtrl_,cellCub_,false);
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
          if (bdryCellNodes[i][j] != ROL::nullPtr) {
            feVelBdry_[i][j] = ROL::makePtr<FE<Real>>(bdryCellNodes[i][j],basisPtrVel_,bdryCub_,j);
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
    fieldHelper_ = ROL::makePtr<FieldHelper<Real>>(numFields_, numDofs_, numFieldDofs_, fieldPattern_);
  }

  const ROL::Ptr<FE<Real>> getVelocityFE(void) const {
    return feVel_;
  }

  const ROL::Ptr<FE<Real>> getPressureFE(void) const {
    return fePrs_;
  }

  const ROL::Ptr<FE<Real>> getControlFE(void) const {
    return feCtrl_;
  }

  const std::vector<ROL::Ptr<FE<Real>>> getVelocityBdryFE(const int sideset = -1) const {
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

  const ROL::Ptr<FieldHelper<Real>> getFieldHelper(void) const {
    return fieldHelper_;
  }

  const ROL::Ptr<const FieldUtils::FieldInfo> getStateFieldInfo(void) const {
    return fieldInfo_;
  }

  const ROL::Ptr<const FieldUtils::FieldInfo> getControlFieldInfo(void) const {
    return fieldInfoCtrl_;
  }

  void printData(std::string tag,
                 const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                 const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                 const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c = fePrs_->gradN()->dimension(0);
    const int p = fePrs_->gradN()->dimension(2);
    const int d = fePrs_->gradN()->dimension(3);
    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> valCtrl, alpha;
    valCtrl = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alpha   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
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
        for (int k = 0; k < d; ++k) filePer << std::setw(25) << (*fePrs_->cubPts())(i,j,k);
        filePer << std::setw(25) << viscosity_/(*alpha)(i,j);
        filePer << std::endl;
      }
    }
    filePer.close();
  }

private:

  void applyDirichletToResidual(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> &R,
                          const std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> &U) const {
    const int d = feVel_->gradN()->dimension(3);
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
                  (*R[m])(cidx,fvidx_[j][l]) = (*U[m])(cidx,fvidx_[j][l]) - (*bdryCellDofValues_[i][j])(k,fvidx_[j][l],m);
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
                  (*R[m])(cidx,fvidx_[j][l]) = (*U[m])(cidx,fvidx_[j][l]);
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
                (*R[d])(cidx,fpidx_[j][l]) = (*U[d])(cidx,fpidx_[j][l]) - Patm_;
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
                (*R[0])(cidx,fvidx_[j][l]) = (*U[0])(cidx,fvidx_[j][l]);
              }
            }
          }
        }
      }
    }
  }

  void applyDirichletToJacobian1(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> &J) const {
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int d  = feVel_->gradN()->dimension(3);
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
                      (*J[n][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                    }
                    (*J[n][n])(cidx,fvidx_[j][l],fvidx_[j][l]) = static_cast<Real>(1);
                  }
                }
                for (int m=0; m < fp; ++m) {
                  for (int n=0; n < d; ++n) {
                    (*J[n][d])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
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
                      (*J[n][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                    }
                    (*J[n][n])(cidx,fvidx_[j][l],fvidx_[j][l]) = static_cast<Real>(1);
                  }
                }
                for (int m=0; m < fp; ++m) {
                  for (int n=0; n < d; ++n) {
                    (*J[n][d])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
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
                    (*J[d][p])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
                for (int m = 0; m < fp; ++m) {
                  (*J[d][d])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
                }
                (*J[d][d])(cidx,fpidx_[j][l],fpidx_[j][l]) = static_cast<Real>(1);
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
                    (*J[0][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                  (*J[0][0])(cidx,fvidx_[j][l],fvidx_[j][l]) = static_cast<Real>(1);
                }
                for (int m=0; m < fp; ++m) {
                  (*J[0][d])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                }
              }
            }
          }
        }
      }
    }
  }

  void applyDirichletToJacobian2(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> &J) const {
    const int fc = feCtrl_->gradN()->dimension(1);
    const int d  = feVel_->gradN()->dimension(3);
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
                    (*J[n][0])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
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
                    (*J[n][0])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
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
                  (*J[0][0])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                }
              }
            }
          }
        }
      }
    }
  }

  void applyDirichletToMultiplier(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> &L) const {
    const int d = feVel_->gradN()->dimension(3);
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
                  (*L[m])(cidx,fvidx_[j][l]) = static_cast<Real>(0);
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
                  (*L[m])(cidx,fvidx_[j][l]) = static_cast<Real>(0);
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
                (*L[d])(cidx,fpidx_[j][l]) = static_cast<Real>(0);
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
                (*L[0])(cidx,fvidx_[j][l]) = static_cast<Real>(0);
              }
            }
          }
        }
      }
    }
  }

}; // PDE_NavierStokes

#endif
