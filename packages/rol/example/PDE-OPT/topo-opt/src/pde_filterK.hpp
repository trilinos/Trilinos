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

#ifndef PDE_TOPO_OPT_FILTERK_HPP
#define PDE_TOPO_OPT_FILTERK_HPP

#include "../../TOOLS/pdeK.hpp"
#include "../../TOOLS/feK.hpp"

#include "Intrepid2_HVOL_C0_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"

#include "ROL_Ptr.hpp"


template <class Real, class DeviceType>
class PDE_Filter : public PDE<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;
private:
  // Finite element basis information
  basis_ptr basisPtr_, basisPtrDens_;
  std::vector<basis_ptr> basisPtrs_, basisPtrsDens_;
  // Cell cubature information
  ROL::Ptr<Intrepid2::Cubature<DeviceType,Real,Real>> cellCub_;
  // Cell node information
  scalar_view volCellNodes_;
  // Finite element definition
  ROL::Ptr<fe_type> fe_, feDens_;
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
        basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType, Real, Real>>();
      else if (basisOrder == 2)
        basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<DeviceType, Real, Real>>();
    }
    else if (probDim == 3) {
      if (basisOrder == 1)
        basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_HEX_C1_FEM<DeviceType, Real, Real>>();
      else if (basisOrder == 2)
        basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_HEX_C2_FEM<DeviceType, Real, Real>>();
    }
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();        // get cell type from basis
    if (probDim == 2) {
      if (basisOrderDens == 1)
        basisPtrDens_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType, Real, Real>>();
      else
        basisPtrDens_ = ROL::makePtr<Intrepid2::Basis_HVOL_C0_FEM<DeviceType, Real, Real>>(cellType);
    }
    else if (probDim == 3) {
      if (basisOrderDens == 1)
        basisPtrDens_ = ROL::makePtr<Intrepid2::Basis_HGRAD_HEX_C1_FEM<DeviceType, Real, Real>>();
      else
        basisPtrDens_ = ROL::makePtr<Intrepid2::Basis_HVOL_C0_FEM<DeviceType, Real, Real>>(cellType);
    }
    basisPtrs_.clear(); basisPtrsDens_.clear();
    basisPtrs_.push_back(basisPtr_);  // Filtered Density component
    basisPtrsDens_.push_back(basisPtrDens_); // Density component

    // Quadrature rules.
    Intrepid2::DefaultCubatureFactory cubFactory;                            // create cubature factory
    cellCub_ = cubFactory.create<DeviceType,Real,Real>(cellType, cubDegree); // create default cubature

    // Other problem parameters.
    Real filterRadius = parlist.sublist("Problem").get("Filter Radius",  0.1);
    lengthScale_ = std::pow(filterRadius, 2)/static_cast<Real>(12);
  }

  void residual(scalar_view & res,
                const scalar_view u_coeff,
                const scalar_view z_coeff = scalar_view(),
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);
    const int p = fe_->gradN().extent_int(2);
    const int d = fe_->gradN().extent_int(3);
 
    // Initialize residuals.
    res = scalar_view("res",c,f);

    // Evaluate/interpolate finite element fields on cells.
    scalar_view valU_eval  = scalar_view("valU_eval", c, p);
    scalar_view valZ_eval  = scalar_view("valZ_eval", c, p);
    scalar_view gradU_eval = scalar_view("gradU_eval", c, p, d);
    fe_->evaluateValue(valU_eval, u_coeff);
    feDens_->evaluateValue(valZ_eval, z_coeff);
    fe_->evaluateGradient(gradU_eval, u_coeff);

    rst::scale(gradU_eval, lengthScale_);
    rst::scale(valZ_eval,  static_cast<Real>(-1));

    /*** Evaluate weak form of the residual. ***/
    fst::integrate(res,gradU_eval,fe_->gradNdetJ(),false);
    fst::integrate(res,valU_eval,fe_->NdetJ(),true);
    fst::integrate(res,valZ_eval,fe_->NdetJ(),true);
  }

  void Jacobian_1(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);
 
    // Initialize Jacobians.
    jac = scalar_view("jac",c,f,f);

    /*** Evaluate weak form of the Jacobian. ***/
    Kokkos::deep_copy(jac,fe_->stiffMat());
    rst::scale(jac, lengthScale_); // ls*gradN1 . gradN2
    rst::add(jac,fe_->massMat());  // + N1 * N2
  }


  void Jacobian_2(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c  = fe_->gradN().extent_int(0);
    const int f  = fe_->gradN().extent_int(1);
    const int fd = feDens_->gradN().extent_int(1);
 
    // Initialize Jacobians.
    jac = scalar_view("jac",c,f,fd);

    /*** Evaluate weak form of the Jacobian. ***/
    fst::integrate(jac,fe_->NdetJ(),feDens_->N(),false);
    rst::scale(jac, static_cast<Real>(-1));
  }

  void Hessian_11(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Filter::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Filter::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Filter::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Filter::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(scalar_view & riesz) override {
    // Retrieve dimensions.
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);
 
    // Initialize Jacobians.
    riesz = scalar_view("riesz1",c,f,f);
    Kokkos::deep_copy(riesz,fe_->stiffMat());
    rst::add(riesz,fe_->massMat());
  }

  void RieszMap_2(scalar_view & riesz) override {
    // Retrieve dimensions.
    const int c = feDens_->gradN().extent_int(0);
    const int f = feDens_->gradN().extent_int(1);
 
    // Initialize Jacobians.
    riesz = scalar_view("riesz2",c,f,f);
    Kokkos::deep_copy(riesz,feDens_->massMat());
  }

  std::vector<basis_ptr> getFields() override {
    return basisPtrs_;
  }

  std::vector<basis_ptr> getFields2() override {
    return basisPtrsDens_;
  }

  void setCellNodes(const scalar_view &volCellNodes,
                    const std::vector<std::vector<scalar_view>> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds) override {
    volCellNodes_ = volCellNodes;
    // Finite element definition.
    fe_     = ROL::makePtr<fe_type>(volCellNodes_,basisPtr_,cellCub_);
    feDens_ = ROL::makePtr<fe_type>(volCellNodes_,basisPtrDens_,cellCub_,false);
  }

  void setFieldPattern(const std::vector<std::vector<int>> &fieldPattern,
                       const std::vector<std::vector<int>> &fieldPattern2) override {}

  const ROL::Ptr<fe_type> getStateFE(void) const {
    return fe_;
  }

  const ROL::Ptr<fe_type> getDensityFE(void) const {
    return feDens_;
  }

}; // PDE_Filter

#endif
