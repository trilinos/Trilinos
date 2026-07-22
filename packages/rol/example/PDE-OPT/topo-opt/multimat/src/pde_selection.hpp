// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde_selection.hpp
    \brief Implements the local selection constraint.
*/

#ifndef PDE_TOPO_OPT_SELECTION_MULTIMAT_HPP
#define PDE_TOPO_OPT_SELECTION_MULTIMAT_HPP

#include "../../../TOOLS/pde.hpp"
#include "../../../TOOLS/fe.hpp"
#include "../../../TOOLS/fieldhelper.hpp"
#include "../../../TOOLS/Intrepid_HGRAD_C0_FEM.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real>
class PDE_MultiMat_Selection : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>> basisPtr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrs_;
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>> basisPtrDens_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrsDens_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real>> cellCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real>> volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real>> fe_, feDens_;
  // Field pattern, offsets, etc.
  std::vector<std::vector<int>> fieldPattern_; // local Field/DOF pattern; set from DOF manager 
  int numFields_;                              // number of fields (equations in the PDE)
  int numDofs_;                                // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                    // for each field, a counting offset
  std::vector<int> numFieldDofs_;              // for each field, number of degrees of freedom
  std::vector<std::vector<int>> fieldPatternDens_; // local Field/DOF pattern; set from DOF manager 
  int numFieldsDens_;                              // number of fields (equations in the PDE)
  int numDofsDens_;                                // total number of degrees of freedom for all (local) fields
  std::vector<int> offsetDens_;                    // for each field, a counting offset
  std::vector<int> numFieldDofsDens_;              // for each field, number of degrees of freedom

  ROL::Ptr<FieldUtils::FieldInfo>             fieldInfo_, fieldInfoDens_;

public:
  PDE_MultiMat_Selection(Teuchos::ParameterList &parlist) {
    // Material properties
    std::vector<Real> ym = ROL::getArrayFromStringParameter<Real>(parlist.sublist("Problem"), "Young's Modulus");
    int T = ym.size();

    // Finite element fields.
    int basisOrderDens = parlist.sublist("Problem").get("Density Basis Order",0);
    int cubDegree      = parlist.sublist("Problem").get("Cubature Degree",4);
    int probDim        = parlist.sublist("Problem").get("Problem Dimension",2);
    if (probDim> 3 || probDim < 2) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> PDE-OPT/binary/elasticity/pde_elasticity.hpp: Problem dimension is not 2 or 3!");
    }
    if (basisOrderDens > 1 || basisOrderDens < 0) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> PDE-OPT/binary/elasticity/pde_elasticity.hpp: Density basis order is not 0 1!");
    }
    if (probDim == 2) {
      if (basisOrderDens == 0) {
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_C0_FEM<Real, Intrepid::FieldContainer<Real>>>();
        basisPtrDens_ = ROL::makePtr<Intrepid::Basis_HGRAD_C0_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
      else if (basisOrderDens == 1) {
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
        basisPtrDens_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
    }
    else if (probDim == 3) {
      if (basisOrderDens == 0) {
        // This does not work in 3D
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_C0_FEM<Real, Intrepid::FieldContainer<Real>>>();
        basisPtrDens_ = ROL::makePtr<Intrepid::Basis_HGRAD_C0_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
      else if (basisOrderDens == 1) {
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_HEX_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
        basisPtrDens_ = ROL::makePtr<Intrepid::Basis_HGRAD_HEX_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
    }
    basisPtrs_.clear(); basisPtrsDens_.clear();
    basisPtrs_.push_back(basisPtr_);  // Displacement component
    for (int i=0; i<T; ++i) {
      basisPtrsDens_.push_back(basisPtrDens_); // Density components
    }

    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology(); // get cell type from basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                // create cubature factory
    cellCub_ = cubFactory.create(cellType, cubDegree);                // create default cubature

    numDofs_ = 0;
    numFields_ = basisPtrs_.size();
    offset_.resize(numFields_);
    numFieldDofs_.resize(numFields_);
    for (int i=0; i<numFields_; ++i) {
      if (i==0) {
        offset_[i]  = 0;
      }
      else {
        offset_[i]  = offset_[i-1] + basisPtrs_[i-1]->getCardinality();
      }
      numFieldDofs_[i] = basisPtrs_[i]->getCardinality();
      numDofs_ += numFieldDofs_[i];
    }

    numDofsDens_ = 0;
    numFieldsDens_ = basisPtrsDens_.size();
    offsetDens_.resize(numFieldsDens_);
    numFieldDofsDens_.resize(numFieldsDens_);
    for (int i=0; i<numFieldsDens_; ++i) {
      if (i==0) {
        offsetDens_[i]  = 0;
      }
      else {
        offsetDens_[i]  = offsetDens_[i-1] + basisPtrsDens_[i-1]->getCardinality();
      }
      numFieldDofsDens_[i] = basisPtrsDens_[i]->getCardinality();
      numDofsDens_ += numFieldDofsDens_[i];
    }
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real>> & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int T = numFieldsDens_;
 
    // Initialize residuals.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> R(1);
    R[0] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f);

    // Split z_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> Z;
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoDens_);

    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> sum, valZ_eval;
    sum       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    sum->initialize(static_cast<Real>(-1));
    for (int t=0; t<T; ++t) {
      valZ_eval->initialize();
      feDens_->evaluateValue(valZ_eval, Z[t]);
      Intrepid::RealSpaceTools<Real>::add(*sum,*valZ_eval);
    }
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[0],
                                                  *sum,          // sum z[t]
                                                  *fe_->NdetJ(), // N
                                                  Intrepid::COMP_CPP,
                                                  false);

    // Combine the residuals.
    FieldUtils::combineFieldCoeff<Real>(res, R, fieldInfo_);
  }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_TopoOpt::Jacobian_1): Jacobian is zero.");
  }


  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c  = fe_->gradN()->dimension(0);
    int f  = fe_->gradN()->dimension(1);
    int fd = feDens_->gradN()->dimension(1);
    int T  = numFieldsDens_;

    // Initialize Jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J(1);
    for (int j=0; j<T; ++j) {
      J[0].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,fd));
    }

    // Evaluate/interpolate finite element fields on cells.
    for (int t=0; t<T; ++t) {
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][t],
                                                    *fe_->NdetJ(),     // N
                                                    *feDens_->N(),     // N
                                                    Intrepid::COMP_CPP,
                                                    false);
    }

    // Combine the jacobians.
    FieldUtils::combineFieldCoeff<Real>(jac, J, fieldInfo_, fieldInfoDens_);
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_TopoOpt::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_TopoOpt::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_TopoOpt::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_TopoOpt::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    //throw Exception::NotImplemented(">>> (PDE_TopoOpt::RieszMap_1): Not implemented.");

    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
 
    // Initialize Riesz maps.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> R(1);
    R[0].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f));

    *(R[0][0]) = *(fe_->stiffMat());

    // Combine the Riesz maps.
    FieldUtils::combineFieldCoeff<Real>(riesz, R, fieldInfo_, fieldInfo_);
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    //throw Exception::NotImplemented(">>> (PDE_TopoOpt::RieszMap_2): Not implemented.");

    // Retrieve dimensions.
    int c = feDens_->gradN()->dimension(0);
    int f = feDens_->gradN()->dimension(1);
    int T = numFieldsDens_;
 
    // Initialize Riesz maps.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> R(T);
    for (int i=0; i<T; ++i) {
      for (int j=0; j<T; ++j) {
        R[i].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f));
      }
    }

    for (int i=0; i<T; ++i) {
      *(R[i][i]) = *(fe_->massMat());
    }

    // Combine the Riesz maps.
    FieldUtils::combineFieldCoeff<Real>(riesz, R, fieldInfoDens_, fieldInfoDens_);
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields() {
    return basisPtrs_;
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields2() {
    return basisPtrsDens_;
  }

  void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real>> &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds) {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fe_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtr_,cellCub_,false);
    feDens_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrDens_,cellCub_,false);
  }

  void setFieldPattern(const std::vector<std::vector<int>> &fieldPattern,
                       const std::vector<std::vector<int>> &fieldPattern2) {
    fieldPattern_     = fieldPattern;
    fieldPatternDens_ = fieldPattern2;
    fieldInfo_     = ROL::makePtr<FieldUtils::FieldInfo>(numFields_, numDofs_, numFieldDofs_, fieldPattern_);
    fieldInfoDens_ = ROL::makePtr<FieldUtils::FieldInfo>(numFieldsDens_, numDofsDens_, numFieldDofsDens_, fieldPatternDens_);
  }

  const ROL::Ptr<FE<Real>> getFE(void) const {
    return fe_;
  }

  const ROL::Ptr<FE<Real>> getDensityFE(void) const {
    return feDens_;
  }

  const ROL::Ptr<FieldUtils::FieldInfo> getFieldInfo(void) const {
    return fieldInfo_;
  }

  const ROL::Ptr<FieldUtils::FieldInfo> getDensityFieldInfo(void) const {
    return fieldInfoDens_;
  }

}; // PDE_TopoOpt

#endif
