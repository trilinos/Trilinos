// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde_darcy.hpp
    \brief Implements the local PDE interface for the Darcy porosity optimization problem.
*/

#ifndef PDE_DARCY_HPP
#define PDE_DARCY_HPP

#include "../../../../TOOLS/pde.hpp"
#include "../../../../TOOLS/fe.hpp"
#include "../../../../TOOLS/fieldhelper.hpp"
#include "permeability.hpp"

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
class PDE_Darcy : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>> basisPtrPrs_, basisPtrCtrl_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrs_, basisPtrsCtrl_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real>> cellCub_, bdryCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real>> volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real>> fePrs_, feCtrl_;
  std::vector<std::vector<ROL::Ptr<FE<Real>>>> fePrsBdry_, feCtrlBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fpidx_;
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
  Real Patm_, inVelocity_, dynvisco_;
  bool useNoSlip_;

  ROL::Ptr<Permeability<Real>> perm_;
  ROL::Ptr<FieldUtils::FieldInfo> fieldInfo_, fieldInfoCtrl_;

  void multiplyByRadius(ROL::Ptr<Intrepid::FieldContainer<Real>> &input,
                        const ROL::Ptr<const Intrepid::FieldContainer<Real>> &cubPts,
                        bool isField, bool useReciprocal = false) const {
    int c = cubPts->dimension(0);
    int p = cubPts->dimension(1);
    Intrepid::FieldContainer<Real> r(c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        if (!useReciprocal) r(i,j) = (*cubPts)(i,j,0);
        else                r(i,j) = static_cast<Real>(1)/(*cubPts)(i,j,0);
      }
    }
    Intrepid::FieldContainer<Real> in = *input;
    if (!isField)
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*input,r,in);
    else
      Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*input,r,in);
  }

public:
  PDE_Darcy(Teuchos::ParameterList &parlist) : Patm_(101.325 /* kg/mm-s^2 */) {
    // Finite element fields.
    int cubDegree     = parlist.sublist("Problem").get("Cubature Degree",4);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2);
    int basisDegPres  = parlist.sublist("Problem").get("Pressure Basis Degree",1);
    int basisDegCtrl  = parlist.sublist("Problem").get("Density Basis Degree",1);
    std::string elemType = parlist.sublist("Problem").get("Element Type","QUAD");
    if (elemType == "TRI") {
      if (basisDegPres == 2)
        basisPtrPrs_ = ROL::makePtr<Intrepid::Basis_HGRAD_TRI_C2_FEM<Real, Intrepid::FieldContainer<Real>>>();
      else 
        basisPtrPrs_ = ROL::makePtr<Intrepid::Basis_HGRAD_TRI_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      if (basisDegCtrl == 1)
        basisPtrCtrl_ = ROL::makePtr<Intrepid::Basis_HGRAD_TRI_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      else
        basisPtrCtrl_ = ROL::makePtr<Intrepid::Basis_HGRAD_C0_FEM<Real, Intrepid::FieldContainer<Real>>>();
    }
    else {
      if (basisDegPres == 2)
        basisPtrPrs_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real>>>();
      else 
        basisPtrPrs_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      if (basisDegCtrl == 1)
        basisPtrCtrl_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      else
        basisPtrCtrl_ = ROL::makePtr<Intrepid::Basis_HGRAD_C0_FEM<Real, Intrepid::FieldContainer<Real>>>();
    }
    basisPtrs_.push_back(basisPtrPrs_);       // Pressure component
    basisPtrsCtrl_.push_back(basisPtrCtrl_);  // Control component
    // Quadrature rules.
    shards::CellTopology cellType = basisPtrs_[0]->getBaseCellTopology();        // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                           // create cubature factory
    cellCub_ = cubFactory.create(cellType, cubDegree);                           // create default cubature

    shards::CellTopology bdryCellType = cellType.getCellTopologyData(1, 0);
    bdryCub_ = cubFactory.create(bdryCellType, bdryCubDegree);

    // Other problem parameters.
    dynvisco_     = parlist.sublist("Problem").get("Dynamic Viscosity", 0.84e-8); // kg/mm-s
    Real density  = parlist.sublist("Problem").get("Fluid Density", 8.988e-11);   // kg/mm^3
    Real inRadius = parlist.sublist("Problem").get("Inlet Radius", 0.6875);       // mm
    Real Reynolds = parlist.sublist("Problem").get("Reynolds Number", 50.0);      // dimensionless
    inVelocity_ = Reynolds * dynvisco_ / (density*static_cast<Real>(2)*inRadius); // mm/s
    perm_       = ROL::makePtr<Permeability<Real>>(parlist);
    useNoSlip_  = parlist.sublist("Problem").get("Use No Slip", false);

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
    const int c = fePrs_->gradN()->dimension(0);
    const int f = fePrs_->gradN()->dimension(1);
    const int p = fePrs_->gradN()->dimension(2);
    const int d = fePrs_->gradN()->dimension(3);
    // Initialize residuals.
    res = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> gradP, AgradP, valCtrl, alpha;
    gradP   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    AgradP  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    valCtrl = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alpha   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fePrs_->evaluateGradient(gradP, u_coeff);
    feCtrl_->evaluateValue(valCtrl, z_coeff);
    perm_->compute(alpha, valCtrl, fePrs_->cubPts(), 0);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*AgradP,*alpha,*gradP);
    /*** Evaluate weak form of the residual. ***/
    multiplyByRadius(AgradP,fePrs_->cubPts(),false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *AgradP,              // p
                                                  *fePrs_->gradNdetJ(), // dPhi/dx
                                                  Intrepid::COMP_CPP,
                                                  false);
    // Boundary conditions.
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      // APPLY DIRICHLET CONDITIONS
      for (int i = 0; i < numSideSets; ++i) {
        if (i==0) /* in flow */ {
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            if (numCellsSide > 0) {
              const int numCubPerSide = fePrsBdry_[i][j]->cubPts()->dimension(1);
              ROL::Ptr<Intrepid::FieldContainer<Real>> nRes, nVal;
              nRes = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f);
              nVal = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
              nVal->initialize(-inVelocity_);
              multiplyByRadius(nVal,fePrsBdry_[i][j]->cubPts(),false);
              Intrepid::FunctionSpaceTools::integrate<Real>(*nRes,
                                                            *nVal,
                                                            *fePrsBdry_[i][j]->NdetJ(),
                                                            Intrepid::COMP_CPP,
                                                            false);
              for (int k = 0; k < numCellsSide; ++k) {
                const int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < f; ++l) {
                  (*res)(cidx,l) += (*nRes)(k,l);
                }
              }
            }
          }
        }
        else if (useNoSlip_ && i==1) /* no slip */ {}
        else if (i==2) /* out flow */ {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fpidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                (*res)(cidx,fpidx_[j][l]) = (*u_coeff)(cidx,fpidx_[j][l]) - Patm_;
              }
            }
          }
        }
        else if (i==3) /* no normal */ {} 
      }
    }
  }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c = fePrs_->gradN()->dimension(0);
    const int f = fePrs_->gradN()->dimension(1);
    const int p = fePrs_->gradN()->dimension(2);
    const int d = fePrs_->gradN()->dimension(3);
    // Initialize Jacobians.
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> alpha, alphaPhi, valCtrl;
    valCtrl  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alpha    = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alphaPhi = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p, d);
    feCtrl_->evaluateValue(valCtrl, z_coeff);
    perm_->compute(alpha, valCtrl, fePrs_->cubPts(), 0);
    // Multiply velocity with alpha
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*alphaPhi, *alpha, *fePrs_->gradN());
    /*** Evaluate weak form of the Jacobian. ***/
    multiplyByRadius(alphaPhi,fePrs_->cubPts(),true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *alphaPhi,            // alpha dPhi/dx
                                                  *fePrs_->gradNdetJ(), // dPhi/dx
                                                  Intrepid::COMP_CPP,
                                                  false);
    // APPLY DIRICHLET CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==0) /* in flow */ {}
        else if (useNoSlip_ && i==1) /* no slip */ {}
        else if (i==2) /* out flow */ {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fpidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m = 0; m < f; ++m) {
                  (*jac)(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
                }
                (*jac)(cidx,fpidx_[j][l],fpidx_[j][l]) = static_cast<Real>(1);
              }
            }
          }
        }
        else if (i==3) /* no normal */ {}
      }
    }
  }

  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c  = fePrs_->gradN()->dimension(0);
    const int f  = fePrs_->gradN()->dimension(1);
    const int fc = feCtrl_->gradN()->dimension(1);
    const int p  = fePrs_->gradN()->dimension(2);
    const int d  = fePrs_->gradN()->dimension(3);
    // Initialize Jacobians.
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, fc);
    // Evaluate on FE basis.
    ROL::Ptr<Intrepid::FieldContainer<Real>> valZ, alpha, alphaN, gradP, gradPgradN;
    valZ       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alpha      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alphaN     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc, p);
    gradP      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    gradPgradN = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    fePrs_->evaluateGradient(gradP, u_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    perm_->compute(alpha, valZ, fePrs_->cubPts(), 1);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*alphaN,*alpha,*feCtrl_->N());
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*gradPgradN,*gradP,*fePrs_->gradNdetJ());
    /*** Evaluate weak form of the residual. ***/
    multiplyByRadius(gradPgradN,fePrs_->cubPts(),true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *gradPgradN, // dP/dx dPhi/dx
                                                  *alphaN,     // alpha' Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    // APPLY DIRICHLET CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==0) /* in flow */ {}
        else if (useNoSlip_ && i==1) /* no slip */ {}
        else if (i==2) /* out flow */ {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fpidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < fc; ++m) {
                  (*jac)(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
                }
              }
            }
          }
        }
        else if (i==3) /* no normal */ {}
      }
    }
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
    throw Exception::Zero(">>> Hessian_11 is zero!");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c  = fePrs_->gradN()->dimension(0);
    const int f  = fePrs_->gradN()->dimension(1);
    const int fc = feCtrl_->gradN()->dimension(1);
    const int p  = fePrs_->gradN()->dimension(2);
    const int d  = fePrs_->gradN()->dimension(3);
    // Initialize hessians.
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc, f);
    // Apply Dirichlet conditions
    ROL::Ptr<Intrepid::FieldContainer<Real>> l0_coeff;
    l0_coeff = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    *l0_coeff = *l_coeff;
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==0) /* in flow */ {}
        else if (useNoSlip_ && i==1) /* no slip */ {}
        else if (i==2) /* out flow */ {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fpidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                (*l0_coeff)(cidx,fpidx_[j][l]) = static_cast<Real>(0);
              }
            }
          }
        }
        else if (i==3) /* no normal */ {}
      }
    }
    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> valZ, alpha, gradL, alphaL;
    valZ   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alpha  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    gradL  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    alphaL = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc, p, d);
    fePrs_->evaluateGradient(gradL, l0_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    perm_->compute(alpha, valZ, fePrs_->cubPts(), 1);
    /*** Evaluate weak form of the Hessian. ***/
    // Multiply velocity with alpha
    for (int j = 0; j < c; ++j) {
      for (int k = 0; k < fc; ++k) {
        for (int l = 0; l < p; ++l) {
          for (int m = 0; m < d; ++m) {
            (*alphaL)(j,k,l,m) = (*alpha)(j,l) * (*feCtrl_->N())(j,k,l) * (*gradL)(j,l,m);
          }
        }
      }
    }
    multiplyByRadius(alphaL,fePrs_->cubPts(),true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *alphaL,              // alpha' dL/dx Phi
                                                  *fePrs_->gradNdetJ(), // dPhi/dx
                                                  Intrepid::COMP_CPP,
                                                  false);
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
    const int c  = fePrs_->gradN()->dimension(0);
    const int f  = fePrs_->gradN()->dimension(1);
    const int fc = feCtrl_->gradN()->dimension(1);
    const int p  = fePrs_->gradN()->dimension(2);
    const int d  = fePrs_->gradN()->dimension(3);
    // Initialize hessians.
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, fc);
    // Apply Dirichlet conditions
    ROL::Ptr<Intrepid::FieldContainer<Real>> l0_coeff;
    l0_coeff = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    *l0_coeff = *l_coeff;
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==0) /* in flow */ {}
        else if (useNoSlip_ && i==1) /* no slip */ {}
        else if (i==2) /* out flow */ {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fpidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                (*l0_coeff)(cidx,fpidx_[j][l]) = static_cast<Real>(0);
              }
            }
          }
        }
        else if (i==3) /* no normal */ {}
      }
    }
    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> valZ, alpha, gradL, alphaL;
    valZ   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alpha  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    gradL  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    alphaL = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc, p, d);
    fePrs_->evaluateGradient(gradL, l0_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    perm_->compute(alpha, valZ, fePrs_->cubPts(), 1);
    /*** Evaluate weak form of the Hessian. ***/
    for (int j = 0; j < c; ++j) {
      for (int k = 0; k < fc; ++k) {
        for (int l = 0; l < p; ++l) {
          for (int m = 0; m < d; ++m) {
            (*alphaL)(j,k,l,m) = (*alpha)(j,l) * (*feCtrl_->N())(j,k,l) * (*gradL)(j,l,m);
          }
        }
      }
    }
    multiplyByRadius(alphaL,fePrs_->cubPts(),true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *fePrs_->gradNdetJ(), // dPhi/dx
                                                  *alphaL,              // alpha' dL/dx Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c  = fePrs_->gradN()->dimension(0);
    const int fc = feCtrl_->gradN()->dimension(1);
    const int p  = fePrs_->gradN()->dimension(2);
    const int d  = fePrs_->gradN()->dimension(3);
    // Initialize hessians.
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc, fc);
    // Apply Dirichlet conditions
    ROL::Ptr<Intrepid::FieldContainer<Real>> l0_coeff;
    l0_coeff = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    *l0_coeff = *l_coeff;
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==0) /* in flow */ {}
        else if (useNoSlip_ && i==1) /* no slip */ {}
        else if (i==2) /* out flow */ {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fpidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                (*l0_coeff)(cidx,fpidx_[j][l]) = static_cast<Real>(0);
              }
            }
          }
        }
        else if (i==3) /* no normal */ {}
      }
    }
    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> valZ, alpha, alphaUL, gradL, gradU;
    valZ    = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alpha   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alphaUL = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc, p);
    gradL   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    gradU   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fePrs_->evaluateGradient(gradL, l0_coeff);
    fePrs_->evaluateGradient(gradU, u_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    perm_->compute(alpha, valZ, fePrs_->cubPts(), 2);
    // Multiply velocity with alpha
    Real dot(0);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < fc; ++j) {
        for (int k = 0; k < p; ++k) {
          dot = static_cast<Real>(0);
          for (int l = 0; l < d; ++l) {
            dot += (*gradU)(i,k,l) * (*gradL)(i,k,l);
          }
          (*alphaUL)(i,j,k) = (*alpha)(i,k) * (*feCtrl_->N())(i,j,k) * dot;
        }
      }
    }
    /*** Evaluate weak form of the Hessian. ***/
    // Velocity equation.
    multiplyByRadius(alphaUL,fePrs_->cubPts(),true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *alphaUL,          // alpha'' dU/dx dL/dx Phi
                                                  *feCtrl_->NdetJ(), // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
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
    const int c = fePrs_->gradN()->dimension(0);
    const int f = fePrs_->gradN()->dimension(1);
    const int p = fePrs_->gradN()->dimension(2);
    const int d = fePrs_->gradN()->dimension(3);
    // Initialize Jacobians.
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> valPhi, gradPhi;
    valPhi  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    gradPhi = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p, d);
    *valPhi  = *fePrs_->N();
    *gradPhi = *fePrs_->gradN();
    /*** Evaluate weak form of the Jacobian. ***/
    multiplyByRadius(gradPhi,fePrs_->cubPts(),true);
    multiplyByRadius(valPhi,fePrs_->cubPts(),true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*riesz,
                                                  *gradPhi,             // alpha dPhi/dx
                                                  *fePrs_->gradNdetJ(), // dPhi/dx
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*riesz,
                                                  *valPhi,              // alpha dPhi/dx
                                                  *fePrs_->NdetJ(),     // dPhi/dx
                                                  Intrepid::COMP_CPP,
                                                  true);
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) override {
    // Retrieve dimensions.
    const int c  = feCtrl_->gradN()->dimension(0);
    const int f  = feCtrl_->gradN()->dimension(1);
    const int p  = feCtrl_->gradN()->dimension(2);
    // Initialize Jacobians.
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // Evaluate on FE basis.
    ROL::Ptr<Intrepid::FieldContainer<Real>> valPsi;
    valPsi  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    *valPsi = *feCtrl_->N();
    /*** Evaluate weak form of the residual. ***/
    multiplyByRadius(valPsi,feCtrl_->cubPts(),true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*riesz,
                                                  *valPsi,           // Psi
                                                  *feCtrl_->NdetJ(), // Psi
                                                  Intrepid::COMP_CPP,
                                                  false);
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
    fePrs_  = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrPrs_,cellCub_);
    feCtrl_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrCtrl_,cellCub_,false);
    fpidx_  = fePrs_->getBoundaryDofs();
    // Construct control boundary FE
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      fePrsBdry_.resize(numSideSets);
      //feCtrlBdry_.resize(numSideSets);
      for (int i = 0; i < numSideSets; ++i) {
        int numLocSides = bdryCellNodes[i].size();
        fePrsBdry_[i].resize(numLocSides);
        //feCtrlBdry_[i].resize(numLocSides);
        for (int j = 0; j < numLocSides; ++j) {
          if (bdryCellNodes[i][j] != ROL::nullPtr) {
            fePrsBdry_[i][j] = ROL::makePtr<FE<Real>>(bdryCellNodes[i][j],basisPtrPrs_,bdryCub_,j);
            //feCtrlBdry_[i][j] = ROL::makePtr<FE<Real>>(bdryCellNodes[i][j],basisPtrCtrl_,bdryCub_,j);
          }
        }
      }
    }
  }

  void setFieldPattern(const std::vector<std::vector<int>> &fieldPattern,
                       const std::vector<std::vector<int>> &fieldPattern2) override {
    fieldPattern_ = fieldPattern;
    fieldInfo_    = ROL::makePtr<FieldUtils::FieldInfo>(numFields_,numDofs_,numFieldDofs_,fieldPattern_);
    fieldPatternCtrl_ = fieldPattern2;
    fieldInfoCtrl_    = ROL::makePtr<FieldUtils::FieldInfo>(numFieldsCtrl_,numDofsCtrl_,numFieldDofsCtrl_,fieldPatternCtrl_);
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
    ROL::Ptr<Intrepid::FieldContainer<Real>> gradP, AgradP, valCtrl, alpha;
    gradP   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    AgradP  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    valCtrl = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alpha   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fePrs_->evaluateGradient(gradP, u_coeff);
    feCtrl_->evaluateValue(valCtrl, z_coeff);
    perm_->compute(alpha, valCtrl, fePrs_->cubPts(), 0);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*AgradP,*alpha,*gradP);
    // Print to velocity file
    std::stringstream nameVel;
    nameVel << "velocity" << tag << ".txt";
    std::ofstream fileVel;
    fileVel.open(nameVel.str());
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        fileVel << std::scientific << std::setprecision(15);
        for (int k = 0; k < d; ++k) fileVel << std::setw(25) << (*fePrs_->cubPts())(i,j,k);
        for (int k = 0; k < d; ++k) fileVel << std::setw(25) << -(*AgradP)(i,j,k);
        fileVel << std::endl;
      }
    }
    fileVel.close();
    // Print to permeability file
    std::stringstream namePer;
    namePer << "permeability" << tag << ".txt";
    std::ofstream filePer;
    filePer.open(namePer.str());
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        filePer << std::scientific << std::setprecision(15);
        for (int k = 0; k < d; ++k) filePer << std::setw(25) << (*fePrs_->cubPts())(i,j,k);
        filePer << std::setw(25) << dynvisco_*(*alpha)(i,j);
        filePer << std::endl;
      }
    }
    filePer.close();
  }

  const ROL::Ptr<FE<Real>> getPressureFE(void) const {
    return fePrs_;
  }

  const ROL::Ptr<FE<Real>> getControlFE(void) const {
    return feCtrl_;
  }

  const std::vector<ROL::Ptr<FE<Real>>> getPressureBdryFE(const int sideset = -1) const {
    int side = (sideset < 0 || sideset > 4) ? 0 : sideset;
    return fePrsBdry_[side];
  }

  const std::vector<ROL::Ptr<FE<Real>>> getControlBdryFE(const int sideset = -1) const {
    int side = (sideset < 0 || sideset > 4) ? 0 : sideset;
    return feCtrlBdry_[side];
  }

  const std::vector<std::vector<int>> getBdryCellLocIds(const int sideset = -1) const {
    int side = (sideset < 0 || sideset > 4) ? 0 : sideset;
    return bdryCellLocIds_[side];
  }

  const ROL::Ptr<Permeability<Real>> getPermeability(void) const {
    return perm_;
  }

  const ROL::Ptr<const FieldUtils::FieldInfo> getStateFieldInfo(void) const {
    return fieldInfo_;
  }

  const ROL::Ptr<const FieldUtils::FieldInfo> getControlFieldInfo(void) const {
    return fieldInfoCtrl_;
  }

}; // PDE_Darcy

#endif
