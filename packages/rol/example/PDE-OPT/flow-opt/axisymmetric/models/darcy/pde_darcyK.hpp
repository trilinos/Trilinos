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

#ifndef PDE_DARCYK_HPP
#define PDE_DARCYK_HPP

#include "../../../../TOOLS/pdeK.hpp"
#include "../../../../TOOLS/feK.hpp"
#include "../../../../TOOLS/fieldhelperK.hpp"
#include "permeabilityK.hpp"

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
class PDE_Darcy : public PDE<Real, DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  // Finite element basis information
  basis_ptr basisPtrPrs_, basisPtrCtrl_;
  std::vector<basis_ptr> basisPtrs_, basisPtrsCtrl_;
  // Cell cubature information
  ROL::Ptr<Intrepid2::Cubature<DeviceType,Real,Real>> cellCub_, bdryCub_;
  // Cell node information
  scalar_view volCellNodes_;
  std::vector<std::vector<scalar_view>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<fe_type> fePrs_, feCtrl_;
  std::vector<std::vector<ROL::Ptr<fe_type>>> fePrsBdry_, feCtrlBdry_;
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

  ROL::Ptr<Permeability<Real,DeviceType>> perm_;
  ROL::Ptr<FieldUtils::FieldInfo> fieldInfo_, fieldInfoCtrl_;

  void multiplyByRadius(scalar_view &input, const scalar_view cubPts,
                        bool isField, bool useReciprocal = false) const {
    int c = cubPts.extent_int(0);
    int p = cubPts.extent_int(1);
    scalar_view r("r",c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        if (!useReciprocal) r(i,j) = cubPts(i,j,0);
        else                r(i,j) = static_cast<Real>(1)/cubPts(i,j,0);
      }
    }
    scalar_view in(input);
    if (!isField)
      fst::scalarMultiplyDataData(input,r,in);
    else
      fst::scalarMultiplyDataField(input,r,in);
  }

public:
  PDE_Darcy(ROL::ParameterList &parlist) : Patm_(101.325 /* kg/mm-s^2 */) {
    // Finite element fields.
    int cubDegree     = parlist.sublist("Problem").get("Cubature Degree",4);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2);
    int basisDegPres  = parlist.sublist("Problem").get("Pressure Basis Degree",1);
    int basisDegCtrl  = parlist.sublist("Problem").get("Density Basis Degree",1);
    std::string elemType = parlist.sublist("Problem").get("Element Type","QUAD");
    if (elemType == "TRI") {
      if (basisDegPres == 2)
        basisPtrPrs_ = ROL::makePtr<Intrepid2::Basis_HGRAD_TRI_C2_FEM<DeviceType,Real,Real>>();
      else 
        basisPtrPrs_ = ROL::makePtr<Intrepid2::Basis_HGRAD_TRI_C1_FEM<DeviceType,Real,Real>>();
    }
    else {
      if (basisDegPres == 2)
        basisPtrPrs_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<DeviceType,Real,Real>>();
      else 
        basisPtrPrs_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType,Real,Real>>();
    }
    basisPtrs_.push_back(basisPtrPrs_);       // Pressure component
    shards::CellTopology cellType = basisPtrs_[0]->getBaseCellTopology();        // get the cell type from any basis
    if (basisDegCtrl == 1) {
      if (elemType == "TRI")
        basisPtrCtrl_ = ROL::makePtr<Intrepid2::Basis_HGRAD_TRI_C1_FEM<DeviceType,Real,Real>>();
      else
        basisPtrCtrl_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType,Real,Real>>();
    }
    else {
      basisPtrCtrl_ = ROL::makePtr<Intrepid2::Basis_HVOL_C0_FEM<DeviceType,Real,Real>>(cellType);
    }
    basisPtrsCtrl_.push_back(basisPtrCtrl_);  // Control component
    // Quadrature rules.
    Intrepid2::DefaultCubatureFactory cubFactory;                            // create cubature factory
    cellCub_ = cubFactory.create<DeviceType,Real,Real>(cellType, cubDegree); // create default cubature

    shards::CellTopology bdryCellType = cellType.getCellTopologyData(1, 0);
    bdryCub_ = cubFactory.create<DeviceType,Real,Real>(bdryCellType, bdryCubDegree);

    // Other problem parameters.
    dynvisco_     = parlist.sublist("Problem").get("Dynamic Viscosity", 0.84e-8); // kg/mm-s
    Real density  = parlist.sublist("Problem").get("Fluid Density", 8.988e-11);   // kg/mm^3
    Real inRadius = parlist.sublist("Problem").get("Inlet Radius", 0.6875);       // mm
    Real Reynolds = parlist.sublist("Problem").get("Reynolds Number", 50.0);      // dimensionless
    inVelocity_ = Reynolds * dynvisco_ / (density*static_cast<Real>(2)*inRadius); // mm/s
    perm_       = ROL::makePtr<Permeability<Real,DeviceType>>(parlist);
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

  void residual(scalar_view & res,
                const scalar_view u_coeff,
                const scalar_view z_coeff = scalar_view(),
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c = fePrs_->gradN().extent_int(0);
    const int f = fePrs_->gradN().extent_int(1);
    const int p = fePrs_->gradN().extent_int(2);
    const int d = fePrs_->gradN().extent_int(3);
    // Initialize residuals.
    res = scalar_view("res", c, f);
    // Evaluate/interpolate finite element fields on cells.
    scalar_view gradP("gradP", c, p, d);
    scalar_view AgradP("AgradP", c, p, d);
    scalar_view valCtrl("valCtrl", c, p);
    scalar_view alpha("alpha", c, p);
    fePrs_->evaluateGradient(gradP, u_coeff);
    feCtrl_->evaluateValue(valCtrl, z_coeff);
    perm_->compute(alpha, valCtrl, fePrs_->cubPts(), 0);
    fst::scalarMultiplyDataData(AgradP,alpha,gradP);
    /*** Evaluate weak form of the residual. ***/
    multiplyByRadius(AgradP,fePrs_->cubPts(),false);
    fst::integrate(res,AgradP,fePrs_->gradNdetJ(),false);
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
              const int numCubPerSide = fePrsBdry_[i][j]->cubPts().extent_int(1);
              scalar_view nRes("nRes", numCellsSide, f);
              scalar_view nVal("nVal", numCellsSide, numCubPerSide);
	      Kokkos::deep_copy(nVal, -inVelocity_);
              multiplyByRadius(nVal,fePrsBdry_[i][j]->cubPts(),false);
              fst::integrate(nRes,nVal,fePrsBdry_[i][j]->NdetJ(),false);
              for (int k = 0; k < numCellsSide; ++k) {
                const int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < f; ++l) {
                  res(cidx,l) += nRes(k,l);
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
                res(cidx,fpidx_[j][l]) = u_coeff(cidx,fpidx_[j][l]) - Patm_;
              }
            }
          }
        }
        else if (i==3) /* no normal */ {} 
      }
    }
  }

  void Jacobian_1(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c = fePrs_->gradN().extent_int(0);
    const int f = fePrs_->gradN().extent_int(1);
    const int p = fePrs_->gradN().extent_int(2);
    const int d = fePrs_->gradN().extent_int(3);
    // Initialize Jacobians.
    jac = scalar_view("jac", c, f, f);
    // Evaluate/interpolate finite element fields on cells.
    scalar_view valCtrl("valCtrl", c, p);
    scalar_view alpha("alpha", c, p);
    scalar_view alphaPhi("alphaPhi", c, f, p, d);
    feCtrl_->evaluateValue(valCtrl, z_coeff);
    perm_->compute(alpha, valCtrl, fePrs_->cubPts(), 0);
    // Multiply velocity with alpha
    fst::scalarMultiplyDataField(alphaPhi, alpha, fePrs_->gradN());
    /*** Evaluate weak form of the Jacobian. ***/
    multiplyByRadius(alphaPhi,fePrs_->cubPts(),true);
    fst::integrate(jac,alphaPhi,fePrs_->gradNdetJ(),false);
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
                  jac(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
                }
                jac(cidx,fpidx_[j][l],fpidx_[j][l]) = static_cast<Real>(1);
              }
            }
          }
        }
        else if (i==3) /* no normal */ {}
      }
    }
  }

  void Jacobian_2(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c  = fePrs_->gradN().extent_int(0);
    const int f  = fePrs_->gradN().extent_int(1);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int p  = fePrs_->gradN().extent_int(2);
    const int d  = fePrs_->gradN().extent_int(3);
    // Initialize Jacobians.
    jac = scalar_view("jac", c, f, fc);
    // Evaluate on FE basis.
    scalar_view valZ("valZ", c, p);
    scalar_view alpha("alpha", c, p);
    scalar_view alphaN("alphaN", c, fc, p);
    scalar_view gradP("gradP", c, p, d);
    scalar_view gradPgradN("gradPgradN", c, f, p);
    fePrs_->evaluateGradient(gradP, u_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    perm_->compute(alpha, valZ, fePrs_->cubPts(), 1);
    fst::scalarMultiplyDataField(alphaN,alpha,feCtrl_->N());
    fst::dotMultiplyDataField(gradPgradN,gradP,fePrs_->gradNdetJ());
    /*** Evaluate weak form of the residual. ***/
    multiplyByRadius(gradPgradN,fePrs_->cubPts(),true);
    fst::integrate(jac,gradPgradN,alphaN,false);
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
                  jac(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
                }
              }
            }
          }
        }
        else if (i==3) /* no normal */ {}
      }
    }
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
    throw Exception::Zero(">>> Hessian_11 is zero!");
  }

  void Hessian_12(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c  = fePrs_->gradN().extent_int(0);
    const int f  = fePrs_->gradN().extent_int(1);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int p  = fePrs_->gradN().extent_int(2);
    const int d  = fePrs_->gradN().extent_int(3);
    // Initialize hessians.
    hess = scalar_view("hess", c, fc, f);
    // Apply Dirichlet conditions
    scalar_view l0_coeff("l0_coeff", c, f);
    Kokkos::deep_copy(l0_coeff,l_coeff);
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
                l0_coeff(cidx,fpidx_[j][l]) = static_cast<Real>(0);
              }
            }
          }
        }
        else if (i==3) /* no normal */ {}
      }
    }
    // Evaluate/interpolate finite element fields on cells.
    scalar_view valZ("valZ", c, p);
    scalar_view alpha("alpha", c, p);
    scalar_view gradL("gradL", c, p, d);
    scalar_view alphaL("alphaL", c, fc, p, d);
    fePrs_->evaluateGradient(gradL, l0_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    perm_->compute(alpha, valZ, fePrs_->cubPts(), 1);
    /*** Evaluate weak form of the Hessian. ***/
    // Multiply velocity with alpha
    for (int j = 0; j < c; ++j) {
      for (int k = 0; k < fc; ++k) {
        for (int l = 0; l < p; ++l) {
          for (int m = 0; m < d; ++m) {
            alphaL(j,k,l,m) = alpha(j,l) * (feCtrl_->N())(j,k,l) * gradL(j,l,m);
          }
        }
      }
    }
    multiplyByRadius(alphaL,fePrs_->cubPts(),true);
    fst::integrate(hess,alphaL,fePrs_->gradNdetJ(),false);
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
    const int c  = fePrs_->gradN().extent_int(0);
    const int f  = fePrs_->gradN().extent_int(1);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int p  = fePrs_->gradN().extent_int(2);
    const int d  = fePrs_->gradN().extent_int(3);
    // Initialize hessians.
    hess = scalar_view("hess", c, f, fc);
    // Apply Dirichlet conditions
    scalar_view l0_coeff("l0_coeff", c, f);
    Kokkos::deep_copy(l0_coeff,l_coeff);
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
                l0_coeff(cidx,fpidx_[j][l]) = static_cast<Real>(0);
              }
            }
          }
        }
        else if (i==3) /* no normal */ {}
      }
    }
    // Evaluate/interpolate finite element fields on cells.
    scalar_view valZ("valZ", c, p);
    scalar_view alpha("alpha", c, p);
    scalar_view gradL("gradL", c, p, d);
    scalar_view alphaL("alphaL", c, fc, p, d);
    fePrs_->evaluateGradient(gradL, l0_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    perm_->compute(alpha, valZ, fePrs_->cubPts(), 1);
    /*** Evaluate weak form of the Hessian. ***/
    for (int j = 0; j < c; ++j) {
      for (int k = 0; k < fc; ++k) {
        for (int l = 0; l < p; ++l) {
          for (int m = 0; m < d; ++m) {
            alphaL(j,k,l,m) = alpha(j,l) * (feCtrl_->N())(j,k,l) * gradL(j,l,m);
          }
        }
      }
    }
    multiplyByRadius(alphaL,fePrs_->cubPts(),true);
    fst::integrate(hess,fePrs_->gradNdetJ(),alphaL,false);
  }

  void Hessian_22(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c  = fePrs_->gradN().extent_int(0);
    const int f  = fePrs_->gradN().extent_int(1);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int p  = fePrs_->gradN().extent_int(2);
    const int d  = fePrs_->gradN().extent_int(3);
    // Initialize hessians.
    hess = scalar_view("hess", c, fc, fc);
    // Apply Dirichlet conditions
    scalar_view l0_coeff("l0_coeff", c, f);
    Kokkos::deep_copy(l0_coeff,l_coeff);
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
                l0_coeff(cidx,fpidx_[j][l]) = static_cast<Real>(0);
              }
            }
          }
        }
        else if (i==3) /* no normal */ {}
      }
    }
    // Evaluate/interpolate finite element fields on cells.
    scalar_view valZ("valZ", c, p);
    scalar_view alpha("alpha", c, p);
    scalar_view alphaUL("alphaUL", c, fc, p);
    scalar_view gradL("gradL", c, p, d);
    scalar_view gradU("gradU", c, p, d);
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
            dot += gradU(i,k,l) * gradL(i,k,l);
          }
          alphaUL(i,j,k) = alpha(i,k) * (feCtrl_->N())(i,j,k) * dot;
        }
      }
    }
    /*** Evaluate weak form of the Hessian. ***/
    // Velocity equation.
    multiplyByRadius(alphaUL,fePrs_->cubPts(),true);
    fst::integrate(hess,alphaUL,feCtrl_->NdetJ(),false);
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
    const int c = fePrs_->gradN().extent_int(0);
    const int f = fePrs_->gradN().extent_int(1);
    const int p = fePrs_->gradN().extent_int(2);
    const int d = fePrs_->gradN().extent_int(3);
    // Initialize Jacobians.
    riesz = scalar_view("riesz", c, f, f);
    // Evaluate/interpolate finite element fields on cells.
    scalar_view valPhi("valPhi", c, f, p);
    scalar_view gradPhi("gradPhi", c, f, p, d);
    Kokkos::deep_copy(valPhi, fePrs_->N());
    Kokkos::deep_copy(gradPhi, fePrs_->gradN());
    /*** Evaluate weak form of the Jacobian. ***/
    multiplyByRadius(gradPhi,fePrs_->cubPts(),true);
    multiplyByRadius(valPhi,fePrs_->cubPts(),true);
    fst::integrate(riesz,gradPhi,fePrs_->gradNdetJ(),false);
    fst::integrate(riesz,valPhi,fePrs_->NdetJ(),true);
  }

  void RieszMap_2(scalar_view & riesz) override {
    // Retrieve dimensions.
    const int c = feCtrl_->gradN().extent_int(0);
    const int f = feCtrl_->gradN().extent_int(1);
    const int p = feCtrl_->gradN().extent_int(2);
    // Initialize Jacobians.
    riesz = scalar_view("riesz", c, f, f);
    // Evaluate on FE basis.
    scalar_view valPsi("valPsi", c, f, p);
    Kokkos::deep_copy(valPsi, feCtrl_->N());
    /*** Evaluate weak form of the residual. ***/
    multiplyByRadius(valPsi,feCtrl_->cubPts(),true);
    fst::integrate(riesz,valPsi,feCtrl_->NdetJ(),false);
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
    fePrs_  = ROL::makePtr<fe_type>(volCellNodes_,basisPtrPrs_,cellCub_);
    feCtrl_ = ROL::makePtr<fe_type>(volCellNodes_,basisPtrCtrl_,cellCub_,false);
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
          if (bdryCellNodes[i][j] != scalar_view()) {
            fePrsBdry_[i][j] = ROL::makePtr<fe_type>(bdryCellNodes[i][j],basisPtrPrs_,bdryCub_,j);
            //feCtrlBdry_[i][j] = ROL::makePtr<fe_type>(bdryCellNodes[i][j],basisPtrCtrl_,bdryCub_,j);
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
                 const scalar_view u_coeff,
                 const scalar_view z_coeff = scalar_view(),
                 const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c = fePrs_->gradN().extent_int(0);
    const int p = fePrs_->gradN().extent_int(2);
    const int d = fePrs_->gradN().extent_int(3);
    // Evaluate/interpolate finite element fields on cells.
    scalar_view gradP("gradP", c, p, d);
    scalar_view AgradP("AgradP", c, p, d);
    scalar_view valCtrl("valCtrl", c, p);
    scalar_view alpha("alpha", c, p);
    fePrs_->evaluateGradient(gradP, u_coeff);
    feCtrl_->evaluateValue(valCtrl, z_coeff);
    perm_->compute(alpha, valCtrl, fePrs_->cubPts(), 0);
    fst::scalarMultiplyDataData(AgradP,alpha,gradP);
    // Print to velocity file
    std::stringstream nameVel;
    nameVel << "velocity" << tag << ".txt";
    std::ofstream fileVel;
    fileVel.open(nameVel.str());
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        fileVel << std::scientific << std::setprecision(15);
        for (int k = 0; k < d; ++k) fileVel << std::setw(25) << (fePrs_->cubPts())(i,j,k);
        for (int k = 0; k < d; ++k) fileVel << std::setw(25) << -AgradP(i,j,k);
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
        for (int k = 0; k < d; ++k) filePer << std::setw(25) << (fePrs_->cubPts())(i,j,k);
        filePer << std::setw(25) << dynvisco_*alpha(i,j);
        filePer << std::endl;
      }
    }
    filePer.close();
  }

  const ROL::Ptr<fe_type> getPressureFE(void) const {
    return fePrs_;
  }

  const ROL::Ptr<fe_type> getControlFE(void) const {
    return feCtrl_;
  }

  const std::vector<ROL::Ptr<fe_type>> getPressureBdryFE(const int sideset = -1) const {
    int side = (sideset < 0 || sideset > 4) ? 0 : sideset;
    return fePrsBdry_[side];
  }

  const std::vector<ROL::Ptr<fe_type>> getControlBdryFE(const int sideset = -1) const {
    int side = (sideset < 0 || sideset > 4) ? 0 : sideset;
    return feCtrlBdry_[side];
  }

  const std::vector<std::vector<int>> getBdryCellLocIds(const int sideset = -1) const {
    int side = (sideset < 0 || sideset > 4) ? 0 : sideset;
    return bdryCellLocIds_[side];
  }

  const ROL::Ptr<Permeability<Real,DeviceType>> getPermeability(void) const {
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
