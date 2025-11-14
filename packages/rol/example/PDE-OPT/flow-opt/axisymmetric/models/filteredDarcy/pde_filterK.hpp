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

#ifndef PDE_FILTERK_HPP
#define PDE_FILTERK_HPP

#include "../../../../TOOLS/pdeK.hpp"
#include "../../../../TOOLS/feK.hpp"
#include "../../../../TOOLS/fieldhelperK.hpp"
#include "../../../../TOOLS/Intrepid2_CubatureNodal.hpp"

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
class PDE_Filter : public PDE<Real, DeviceType> {
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
  ROL::Ptr<Intrepid2::Cubature<DeviceType,Real,Real>> cellCub_;
  // Cell node information
  scalar_view volCellNodes_;
  std::vector<std::vector<scalar_view>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<fe_type> feFilter_, feCtrl_;
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
  PDE_Filter(ROL::ParameterList &parlist) {
    // Finite element fields.
    int cubDegree     = parlist.sublist("Problem").get("Filter Cubature Degree",4);
    int basisDegPres  = parlist.sublist("Problem").get("Filter Basis Degree",1);
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
    //cellCub_ = cubFactory.create<DeviceType,Real,Real>(cellType, cubDegree); // create default cubature
    if (cubDegree == -1) {  // nodal cubature
      cellCub_ = ROL::makePtr<Intrepid2::CubatureNodal<DeviceType,Real,Real>>(cellType);
    }
    else {                  // default cubature
      cellCub_ = cubFactory.create<DeviceType,Real,Real>(cellType, cubDegree); // create default cubature
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

  void residual(scalar_view & res,
                const scalar_view u_coeff,
                const scalar_view z_coeff = scalar_view(),
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c = feFilter_->gradN().extent_int(0);
    const int f = feFilter_->gradN().extent_int(1);
    const int p = feFilter_->gradN().extent_int(2);
    const int d = feFilter_->gradN().extent_int(3);
    // Initialize residuals.
    res = scalar_view("res", c, f);
    // Evaluate/interpolate finite element fields on cells.
    scalar_view valP("valP", c, p);
    scalar_view gradP("gradP", c, p, d);
    scalar_view valCtrl("valCtrl", c, p);
    feFilter_->evaluateValue(valP, u_coeff);
    feFilter_->evaluateGradient(gradP, u_coeff);
    feCtrl_->evaluateValue(valCtrl, z_coeff);
    rst::scale(gradP,filterRadius2_);
    rst::scale(valCtrl,static_cast<Real>(-1));
    /*** Evaluate weak form of the residual. ***/
    multiplyByRadius(valP,feFilter_->cubPts(),false);
    multiplyByRadius(gradP,feFilter_->cubPts(),false);
    multiplyByRadius(valCtrl,feFilter_->cubPts(),false);
    fst::integrate(res,gradP,feFilter_->gradNdetJ(),false);
    fst::integrate(res,valP,feFilter_->NdetJ(),true);
    fst::integrate(res,valCtrl,feFilter_->NdetJ(),true);
  }

  void Jacobian_1(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c = feFilter_->gradN().extent_int(0);
    const int f = feFilter_->gradN().extent_int(1);
    const int p = feFilter_->gradN().extent_int(2);
    const int d = feFilter_->gradN().extent_int(3);
    // Initialize Jacobians.
    jac = scalar_view("jac", c, f, f);
    // Evaluate/interpolate finite element fields on cells.
    scalar_view valPhi("valPhi", c, f, p);
    scalar_view gradPhi("gradPhi", c, f, p, d);
    Kokkos::deep_copy(valPhi,feFilter_->N());
    rst::scale(gradPhi, feFilter_->gradN(), filterRadius2_);
    /*** Evaluate weak form of the Jacobian. ***/
    multiplyByRadius(gradPhi,feFilter_->cubPts(),true);
    multiplyByRadius(valPhi,feFilter_->cubPts(),true);
    fst::integrate(jac,gradPhi,feFilter_->gradNdetJ(),false);
    fst::integrate(jac,valPhi,feFilter_->NdetJ(),true);
  }

  void Jacobian_2(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c  = feFilter_->gradN().extent_int(0);
    const int f  = feFilter_->gradN().extent_int(1);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int p  = feFilter_->gradN().extent_int(2);
    // Initialize Jacobians.
    jac = scalar_view("jac", c, f, fc);
    // Evaluate on FE basis.
    scalar_view valPsi("valPsi", c, fc, p);
    rst::scale(valPsi,feCtrl_->N(),static_cast<Real>(-1));
    /*** Evaluate weak form of the residual. ***/
    multiplyByRadius(valPsi,feCtrl_->cubPts(),true);
    fst::integrate(jac,feFilter_->NdetJ(),valPsi,false);
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
    throw Exception::Zero(">>> Hessian_12 is zero!");
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
    throw Exception::Zero(">>> Hessian_21 is zero!");
  }

  void Hessian_22(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> Hessian_22 is zero!");
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
    const int c = feFilter_->gradN().extent_int(0);
    const int f = feFilter_->gradN().extent_int(1);
    const int p = feFilter_->gradN().extent_int(2);
    const int d = feFilter_->gradN().extent_int(3);
    // Initialize Jacobians.
    riesz = scalar_view("riesz1", c, f, f);
    // Evaluate/interpolate finite element fields on cells.
    scalar_view valPhi("valPhi", c, f, p);
    scalar_view gradPhi("gradPhi", c, f, p, d);
    Kokkos::deep_copy(valPhi, feFilter_->N());
    rst::scale(gradPhi, feFilter_->gradN(), filterRadius2_);
    /*** Evaluate weak form of the Jacobian. ***/
    multiplyByRadius(gradPhi,feFilter_->cubPts(),true);
    multiplyByRadius(valPhi,feFilter_->cubPts(),true);
    fst::integrate(riesz,gradPhi,feFilter_->gradNdetJ(),false);
    fst::integrate(riesz,valPhi,feFilter_->NdetJ(),true);
  }

  void RieszMap_2(scalar_view & riesz) override {
    // Retrieve dimensions.
    const int c  = feCtrl_->gradN().extent_int(0);
    const int f  = feCtrl_->gradN().extent_int(1);
    const int p  = feCtrl_->gradN().extent_int(2);
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
    feFilter_ = ROL::makePtr<fe_type>(volCellNodes_,basisPtrPrs_,cellCub_);
    feCtrl_   = ROL::makePtr<fe_type>(volCellNodes_,basisPtrCtrl_,cellCub_,false);
  }

  void setFieldPattern(const std::vector<std::vector<int>> &fieldPattern,
                       const std::vector<std::vector<int>> &fieldPattern2) override {
    fieldPattern_ = fieldPattern;
    fieldInfo_    = ROL::makePtr<FieldUtils::FieldInfo>(numFields_,numDofs_,numFieldDofs_,fieldPattern_);
    fieldPatternCtrl_ = fieldPattern2;
    fieldInfoCtrl_    = ROL::makePtr<FieldUtils::FieldInfo>(numFieldsCtrl_,numDofsCtrl_,numFieldDofsCtrl_,fieldPatternCtrl_);
  }

  const ROL::Ptr<fe_type> getFilterFE() const {
    return feFilter_;
  }

  const ROL::Ptr<fe_type> getControlFE() const {
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
