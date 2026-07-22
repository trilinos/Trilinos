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

#ifndef PDE_TOPO_OPT_FILTER_HPP
#define PDE_TOPO_OPT_FILTER_HPP

#include "../../../../TOOLS/pde.hpp"
#include "../../../../TOOLS/fe.hpp"
#include "../../../../TOOLS/Intrepid_HGRAD_C0_FEM.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"


template <class Real>
class PDE_Filter : public PDE<Real> {
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

public:
  PDE_Filter(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder     = parlist.sublist("Problem").get("Filter Basis Order",1);
    int basisOrderDens = parlist.sublist("Problem").get("Density Basis Order",0);
    int cubDegree      = parlist.sublist("Problem").get("Cubature Degree",4);
//    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2);
    int probDim        = parlist.sublist("Problem").get("Problem Dimension",2);
    TEUCHOS_TEST_FOR_EXCEPTION(probDim>3||probDim<2, std::invalid_argument,
      ">>> PDE-OPT/poisson/pde_poisson.hpp: Problem dimension is not 2 or 3!");
    TEUCHOS_TEST_FOR_EXCEPTION(basisOrder>2||basisOrder<1, std::invalid_argument,
      ">>> PDE-OPT/poisson/pde_poisson.hpp: Basis order is not 1 or 2!");
    if (probDim == 2) {
      if (basisOrder == 1)
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      else if (basisOrder == 2)
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real>>>();
      if (basisOrderDens == 1)
        basisPtrDens_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      else
        basisPtrDens_ = ROL::makePtr<Intrepid::Basis_HGRAD_C0_FEM<Real, Intrepid::FieldContainer<Real>>>();
    }
    else if (probDim == 3) {
      if (basisOrder == 1)
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_HEX_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      else if (basisOrder == 2)
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_HEX_C2_FEM<Real, Intrepid::FieldContainer<Real>>>();
      basisPtrDens_ = ROL::makePtr<Intrepid::Basis_HGRAD_HEX_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
    }
    basisPtrs_.clear(); basisPtrsDens_.clear();
    basisPtrs_.push_back(basisPtr_);  // Filtered Density component
    basisPtrsDens_.push_back(basisPtrDens_); // Density component

    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology(); // get cell type from basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                // create cubature factory
    cellCub_ = cubFactory.create(cellType, cubDegree);                // create default cubature

    // Other problem parameters.
    Real filterRadius = parlist.sublist("Problem").get("Filter Radius",  0.1);
    lengthScale_ = std::pow(filterRadius, 2)/static_cast<Real>(12);
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real>> & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    const int d = fe_->gradN()->dimension(3);
 
    // Initialize residuals.
    res = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f);

    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval, gradU_eval, valZ_eval;
    valU_eval  =  ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    valZ_eval  =  ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    gradU_eval =  ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_->evaluateValue(valU_eval, u_coeff);
    feDens_->evaluateValue(valZ_eval, z_coeff);
    fe_->evaluateGradient(gradU_eval, u_coeff);

    Intrepid::RealSpaceTools<Real>::scale(*gradU_eval, lengthScale_);
    Intrepid::RealSpaceTools<Real>::scale(*valZ_eval,  static_cast<Real>(-1));

    /*** Evaluate weak form of the residual. ***/
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *gradU_eval,        // R*gradU
                                                  *fe_->gradNdetJ(),  // gradN
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *valU_eval,         // U
                                                  *fe_->NdetJ(),      // N
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *valZ_eval,         // -Z
                                                  *fe_->NdetJ(),      // N
                                                  Intrepid::COMP_CPP,
                                                  true);
  }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
 
    // Initialize Jacobians.
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f);

    /*** Evaluate weak form of the Jacobian. ***/
    *jac = *(fe_->stiffMat());
    Intrepid::RealSpaceTools<Real>::scale(*jac, lengthScale_);    // ls*gradN1 . gradN2
    Intrepid::RealSpaceTools<Real>::add(*jac,*(fe_->massMat()));  // + N1 * N2
  }


  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    const int c  = fe_->gradN()->dimension(0);
    const int f  = fe_->gradN()->dimension(1);
    const int fd = feDens_->gradN()->dimension(1);
 
    // Initialize Jacobians.
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,fd);

    /*** Evaluate weak form of the Jacobian. ***/
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *fe_->NdetJ(),  // N1
                                                  *feDens_->N(),  // N2
                                                  Intrepid::COMP_CPP,
                                                  false);
      Intrepid::RealSpaceTools<Real>::scale(*jac, static_cast<Real>(-1));  // -N1 * N2
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Filter::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Filter::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Filter::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Filter::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    //throw Exception::NotImplemented(">>> (PDE_Filter::RieszMap_1): Not implemented.");
    // Retrieve dimensions.
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
 
    // Initialize Jacobians.
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f);
    *riesz = *(fe_->stiffMat());
    Intrepid::RealSpaceTools<Real>::add(*riesz,*(fe_->massMat()));
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    //throw Exception::NotImplemented(">>> (PDE_Filter::RieszMap_2): Not implemented.");
    // Retrieve dimensions.
    const int c = feDens_->gradN()->dimension(0);
    const int f = feDens_->gradN()->dimension(1);
 
    // Initialize Jacobians.
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f);
    *riesz = *(feDens_->massMat());
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
                       const std::vector<std::vector<int>> &fieldPattern2) {}

  const ROL::Ptr<FE<Real>> getStateFE(void) const {
    return fe_;
  }

  const ROL::Ptr<FE<Real>> getDensityFE(void) const {
    return feDens_;
  }

}; // PDE_Filter

#endif
