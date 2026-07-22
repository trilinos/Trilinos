// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde_filter.hpp
    \brief Implements the local PDE interface for the structural topology
           optimization problem.
*/

#ifndef MULTIMAT_PDE_TOPO_OPT_FILTER_HPP
#define MULTIMAT_PDE_TOPO_OPT_FILTER_HPP

#include "../../../TOOLS/pde.hpp"
#include "../../../TOOLS/fe.hpp"
#include "../../../TOOLS/fieldhelper.hpp"
#include "../../../TOOLS/Intrepid_HGRAD_C0_FEM.hpp"

#include "dirichlet.hpp"
#include "load.hpp"
#include "materialtensor.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"


template <class Real>
class MultiMat_PDE_Filter : public PDE<Real> {
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
  // Finite element definition
  ROL::Ptr<FE<Real>> fe_, feDens_;
  // Problem parameters.
  Real lengthScale_;
  // Number of materials.
  int T_;
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
  MultiMat_PDE_Filter(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder     = parlist.sublist("Problem").get("Filter Basis Order",1);
    int basisOrderDens = parlist.sublist("Problem").get("Density Basis Order",0);
    int cubDegree      = parlist.sublist("Problem").get("Cubature Degree",4);
//    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2);
    int probDim        = parlist.sublist("Problem").get("Problem Dimension",2);
    if (probDim > 3 || probDim < 2) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> PDE-OPT/poisson/pde_poisson.hpp: Problem dimension is not 2 or 3!");
    }
    if (basisOrder > 2 || basisOrder < 1) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> PDE-OPT/poisson/pde_poisson.hpp: Basis order is not 1 or 2!");
    }
    if (probDim == 2) {
      if (basisOrder == 1) {
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >>();
      }
      else if (basisOrder == 2) {
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >>();
      }
      if (basisOrderDens == 1) {
        basisPtrDens_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
      else {
        basisPtrDens_ = ROL::makePtr<Intrepid::Basis_HGRAD_C0_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
    }
    else if (probDim == 3) {
      if (basisOrder == 1) {
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_HEX_C1_FEM<Real, Intrepid::FieldContainer<Real> >>();
      }
      else if (basisOrder == 2) {
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_HEX_C2_FEM<Real, Intrepid::FieldContainer<Real> >>();
      }
      basisPtrDens_ = ROL::makePtr<Intrepid::Basis_HGRAD_HEX_C1_FEM<Real, Intrepid::FieldContainer<Real> >>();
    }
    basisPtrs_.clear(); basisPtrsDens_.clear();
    std::vector<Real> ym = ROL::getArrayFromStringParameter<Real>(parlist.sublist("Problem"), "Young's Modulus");
    T_ = ym.size();
    for (int i = 0; i < T_; ++i) {
      basisPtrs_.push_back(basisPtr_);  // Filtered Density component
      basisPtrsDens_.push_back(basisPtrDens_); // Density component
    }

    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology(); // get cell type from basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                // create cubature factory
    cellCub_ = cubFactory.create(cellType, cubDegree);                // create default cubature

    // Other problem parameters.
    Real filterRadius = parlist.sublist("Problem").get("Filter Radius",  0.1);
    lengthScale_ = std::pow(filterRadius, 2)/static_cast<Real>(12);

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

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real> > & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
 
    // Initialize residuals.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> R(T_);
    for (int i=0; i<d; ++i) {
      R[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f);
    }

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U, Z;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoDens_);

    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval, gradU_eval, valZ_eval;
    valU_eval  =  ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    valZ_eval  =  ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    gradU_eval =  ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    for (int i = 0; i < T_; ++i) {
      valU_eval->initialize(); valZ_eval->initialize(); gradU_eval->initialize();
      fe_->evaluateValue(valU_eval, U[i]);
      feDens_->evaluateValue(valZ_eval, Z[i]);
      fe_->evaluateGradient(gradU_eval, U[i]);

      Intrepid::RealSpaceTools<Real>::scale(*gradU_eval, lengthScale_);
      Intrepid::RealSpaceTools<Real>::scale(*valZ_eval,  static_cast<Real>(-1));

      /*** Evaluate weak form of the residual. ***/
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *gradU_eval,        // R*gradU
                                                    *fe_->gradNdetJ(),  // gradN
                                                    Intrepid::COMP_CPP,
                                                    false);
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *valU_eval,         // U
                                                    *fe_->NdetJ(),      // N
                                                    Intrepid::COMP_CPP,
                                                    true);
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *valZ_eval,         // -Z
                                                    *fe_->NdetJ(),      // N
                                                    Intrepid::COMP_CPP,
                                                    true);
    }

    // Combine the residuals.
    FieldUtils::combineFieldCoeff<Real>(res, R, fieldInfo_);
  }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
 
    // Initialize Jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J(T_);
    for (int i=0; i<T_; ++i) {
      for (int j=0; j<T_; ++j) {
        J[i].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f));
      }
    }

    /*** Evaluate weak form of the Jacobian. ***/
    for (int i = 0; i < T_; ++i) {
      *J[i][i] = *fe_->stiffMat();
      Intrepid::RealSpaceTools<Real>::scale(*J[i][i], lengthScale_);  // ls*gradN1 . gradN2
      Intrepid::RealSpaceTools<Real>::add(*J[i][i],*fe_->massMat());  // + N1 * N2
    }

    // Combine the jacobians.
    FieldUtils::combineFieldCoeff<Real>(jac, J, fieldInfo_, fieldInfo_);
  }


  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c  = fe_->gradN()->dimension(0);
    int f  = fe_->gradN()->dimension(1);
    int fd = feDens_->gradN()->dimension(1);

    // Initialize Jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J(T_);
    for (int i=0; i<T_; ++i) {
      for (int j=0; j<T_; ++j) {
        J[i].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,fd));
      }
    }

    /*** Evaluate weak form of the Jacobian. ***/
    for (int i = 0; i < T_; ++i) {
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][i],
                                                    *fe_->NdetJ(),  // N1
                                                    *feDens_->N(),  // N2
                                                    Intrepid::COMP_CPP,
                                                    false);
      Intrepid::RealSpaceTools<Real>::scale(*J[i][i], static_cast<Real>(-1));  // -N1 * N2
    }

    // Combine the jacobians.
    FieldUtils::combineFieldCoeff<Real>(jac, J, fieldInfo_, fieldInfoDens_);
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Filter::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Filter::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Filter::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Filter::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    //throw Exception::NotImplemented(">>> (PDE_Filter::RieszMap_1): Not implemented.");
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
 
    // Initialize Riesz maps.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> R(T_);
    for (int i=0; i<T_; ++i) {
      for (int j=0; j<T_; ++j) {
        R[i].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f));
      }
    }
    for (int i = 0; i < T_; ++i) {
      *R[i][i] = *fe_->stiffMat();
      Intrepid::RealSpaceTools<Real>::add(*R[i][i],*fe_->massMat());
    }

    // Combine the Riesz maps.
    FieldUtils::combineFieldCoeff<Real>(riesz, R, fieldInfo_, fieldInfo_);
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    //throw Exception::NotImplemented(">>> (PDE_Filter::RieszMap_2): Not implemented.");
    // Retrieve dimensions.
    int c = feDens_->gradN()->dimension(0);
    int f = feDens_->gradN()->dimension(1);
 
    // Initialize Riesz maps.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> R(T_);
    for (int i=0; i<T_; ++i) {
      for (int j=0; j<T_; ++j) {
        R[i].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f));
      }
    }
    for (int i = 0; i < T_; ++i) {
      *R[i][i] = *fe_->massMat();
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
    // Finite element definition.
    fe_     = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtr_,cellCub_);
    feDens_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrDens_,cellCub_,false);
  }

  void setFieldPattern(const std::vector<std::vector<int>> &fieldPattern,
                       const std::vector<std::vector<int>> &fieldPattern2) {
    fieldPattern_     = fieldPattern;
    fieldPatternDens_ = fieldPattern2;
    fieldInfo_     = ROL::makePtr<FieldUtils::FieldInfo>(numFields_, numDofs_, numFieldDofs_, fieldPattern_);
    fieldInfoDens_ = ROL::makePtr<FieldUtils::FieldInfo>(numFieldsDens_, numDofsDens_, numFieldDofsDens_, fieldPatternDens_);
  }

  const ROL::Ptr<FE<Real>> getStateFE(void) const {
    return fe_;
  }

  const ROL::Ptr<FE<Real>> getDensityFE(void) const {
    return feDens_;
  }

}; // PDE_Filter

#endif
