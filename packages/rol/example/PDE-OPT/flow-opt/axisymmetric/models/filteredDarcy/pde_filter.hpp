// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde_darcy.hpp
    \brief Implements the local PDE interface for the filtered Darcy porosity optimization problem.
*/

#ifndef PDE_FILTER_HPP
#define PDE_FILTER_HPP

#include "../../../../TOOLS/pde.hpp"
#include "../../../../TOOLS/fe.hpp"
#include "../../../../TOOLS/fieldhelper.hpp"
#include "../../../../TOOLS/Intrepid_HGRAD_C0_FEM.hpp"
#include "../../../../TOOLS/Intrepid_HGRAD_TRI_C0_FEM.hpp"
#include "../../../../TOOLS/Intrepid_CubatureNodal.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid_HGRAD_TRI_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real>
class PDE_Filter : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>> basisPtrPrs_, basisPtrCtrl_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrs_, basisPtrsCtrl_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real>> cellCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real>> volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real>> feFilter_, feCtrl_;
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
  Real filterRadius2_;

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
  PDE_Filter(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int cubDegree     = parlist.sublist("Problem").get("Filter Cubature Degree",4);
    int basisDegPres  = parlist.sublist("Problem").get("Filter Basis Degree",1);
    int basisDegCtrl  = parlist.sublist("Problem").get("Density Basis Degree",1);
    std::string elemType = parlist.sublist("Problem").get("Element Type","QUAD");
    if (elemType == "TRI") {
      if (basisDegPres == 2)
        basisPtrPrs_ = ROL::makePtr<Intrepid::Basis_HGRAD_TRI_C2_FEM<Real, Intrepid::FieldContainer<Real>>>();
      else
        basisPtrPrs_ = ROL::makePtr<Intrepid::Basis_HGRAD_TRI_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      if (basisDegCtrl == 2)
        basisPtrCtrl_ = ROL::makePtr<Intrepid::Basis_HGRAD_TRI_C2_FEM<Real, Intrepid::FieldContainer<Real>>>();
      else if (basisDegCtrl == 1)
        basisPtrCtrl_ = ROL::makePtr<Intrepid::Basis_HGRAD_TRI_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      else
        basisPtrCtrl_ = ROL::makePtr<Intrepid::Basis_HGRAD_TRI_C0_FEM<Real, Intrepid::FieldContainer<Real>>>();
    }
    else {
      if (basisDegPres == 2)
        basisPtrPrs_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real>>>();
      else
        basisPtrPrs_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      if (basisDegCtrl == 2)
        basisPtrCtrl_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real>>>();
      else if (basisDegCtrl == 1)
        basisPtrCtrl_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      else
        basisPtrCtrl_ = ROL::makePtr<Intrepid::Basis_HGRAD_C0_FEM<Real, Intrepid::FieldContainer<Real>>>();
    }
    basisPtrs_.push_back(basisPtrPrs_);       // Pressure component
    basisPtrsCtrl_.push_back(basisPtrCtrl_);  // Control component
    // Quadrature rules.
    shards::CellTopology cellType = basisPtrs_[0]->getBaseCellTopology();        // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                           // create cubature factory
    if (cubDegree == -1) {  // nodal cubature
      cellCub_ = ROL::makePtr<Intrepid::CubatureNodal<Real, Intrepid::FieldContainer<Real>, Intrepid::FieldContainer<Real>>>(cellType);
    }
    else {                  // default cubature
      cellCub_ = cubFactory.create(cellType, cubDegree);
    }

    // Other problem parameters.
    Real filterRadius = parlist.sublist("Problem").get("Filter Radius", 1e-2);
    filterRadius2_ = filterRadius*filterRadius;

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
    const int c = feFilter_->gradN()->dimension(0);
    const int f = feFilter_->gradN()->dimension(1);
    const int p = feFilter_->gradN()->dimension(2);
    const int d = feFilter_->gradN()->dimension(3);
    // Initialize residuals.
    res = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> valP, gradP, valCtrl;
    valP    = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    gradP   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    valCtrl = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    feFilter_->evaluateValue(valP, u_coeff);
    feFilter_->evaluateGradient(gradP, u_coeff);
    feCtrl_->evaluateValue(valCtrl, z_coeff);
    Intrepid::RealSpaceTools<Real>::scale(*gradP,filterRadius2_);
    Intrepid::RealSpaceTools<Real>::scale(*valCtrl,static_cast<Real>(-1));
    /*** Evaluate weak form of the residual. ***/
    multiplyByRadius(valP,feFilter_->cubPts(),false);
    multiplyByRadius(gradP,feFilter_->cubPts(),false);
    multiplyByRadius(valCtrl,feFilter_->cubPts(),false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *gradP,                  // alpha dp/dx
                                                  *feFilter_->gradNdetJ(), // dPhi/dx
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *valP,                   // p
                                                  *feFilter_->NdetJ(),     // Phi
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *valCtrl,                // z
                                                  *feFilter_->NdetJ(),     // Phi
                                                  Intrepid::COMP_CPP,
                                                  true);
  }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c = feFilter_->gradN()->dimension(0);
    const int f = feFilter_->gradN()->dimension(1);
    const int p = feFilter_->gradN()->dimension(2);
    const int d = feFilter_->gradN()->dimension(3);
    // Initialize Jacobians.
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> valPhi, gradPhi;
    valPhi  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    gradPhi = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p, d);
    *valPhi = *feFilter_->N();
    Intrepid::RealSpaceTools<Real>::scale(*gradPhi, *feFilter_->gradN(), filterRadius2_);
    /*** Evaluate weak form of the Jacobian. ***/
    multiplyByRadius(gradPhi,feFilter_->cubPts(),true);
    multiplyByRadius(valPhi,feFilter_->cubPts(),true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *gradPhi,                // alpha dPhi/dx
                                                  *feFilter_->gradNdetJ(), // dPhi/dx
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *valPhi,                 // alpha dPhi/dx
                                                  *feFilter_->NdetJ(),     // dPhi/dx
                                                  Intrepid::COMP_CPP,
                                                  true);
  }

  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c  = feFilter_->gradN()->dimension(0);
    const int f  = feFilter_->gradN()->dimension(1);
    const int fc = feCtrl_->gradN()->dimension(1);
    const int p  = feFilter_->gradN()->dimension(2);
    // Initialize Jacobians.
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, fc);
    // Evaluate on FE basis.
    ROL::Ptr<Intrepid::FieldContainer<Real>> valPsi;
    valPsi  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc, p);
    Intrepid::RealSpaceTools<Real>::scale(*valPsi,*feCtrl_->N(),static_cast<Real>(-1));
    /*** Evaluate weak form of the residual. ***/
    multiplyByRadius(valPsi,feCtrl_->cubPts(),true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *feFilter_->NdetJ(), // Phi
                                                  *valPsi,             // Psi
                                                  Intrepid::COMP_CPP,
                                                  false);
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
    throw Exception::Zero(">>> Hessian_12 is zero!");
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
    throw Exception::Zero(">>> Hessian_21 is zero!");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> Hessian_22 is zero!");
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
    const int c = feFilter_->gradN()->dimension(0);
    const int f = feFilter_->gradN()->dimension(1);
    const int p = feFilter_->gradN()->dimension(2);
    const int d = feFilter_->gradN()->dimension(3);
    // Initialize Jacobians.
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> valPhi, gradPhi;
    valPhi  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    gradPhi = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p, d);
    *valPhi = *feFilter_->N();
    Intrepid::RealSpaceTools<Real>::scale(*gradPhi, *feFilter_->gradN(), filterRadius2_);
    /*** Evaluate weak form of the Jacobian. ***/
    multiplyByRadius(gradPhi,feFilter_->cubPts(),true);
    multiplyByRadius(valPhi,feFilter_->cubPts(),true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*riesz,
                                                  *gradPhi,                // alpha dPhi/dx
                                                  *feFilter_->gradNdetJ(), // dPhi/dx
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*riesz,
                                                  *valPhi,                 // alpha dPhi/dx
                                                  *feFilter_->NdetJ(),     // dPhi/dx
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
    feFilter_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrPrs_,cellCub_);
    feCtrl_   = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrCtrl_,cellCub_,false);
  }

  void setFieldPattern(const std::vector<std::vector<int>> &fieldPattern,
                       const std::vector<std::vector<int>> &fieldPattern2) override {
    fieldPattern_ = fieldPattern;
    fieldInfo_    = ROL::makePtr<FieldUtils::FieldInfo>(numFields_,numDofs_,numFieldDofs_,fieldPattern_);
    fieldPatternCtrl_ = fieldPattern2;
    fieldInfoCtrl_    = ROL::makePtr<FieldUtils::FieldInfo>(numFieldsCtrl_,numDofsCtrl_,numFieldDofsCtrl_,fieldPatternCtrl_);
  }

  const ROL::Ptr<FE<Real>> getFilterFE() const {
    return feFilter_;
  }

  const ROL::Ptr<FE<Real>> getControlFE() const {
    return feCtrl_;
  }

  const ROL::Ptr<const FieldUtils::FieldInfo> getStateFieldInfo(void) const {
    return fieldInfo_;
  }

  const ROL::Ptr<const FieldUtils::FieldInfo> getControlFieldInfo(void) const {
    return fieldInfoCtrl_;
  }

}; // PDE_Filter

#endif
