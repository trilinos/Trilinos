// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde.hpp
    \brief Implements the local PDE interface for the obstacle problem.
*/

#ifndef PDE_POISSON_OBSTACLEK_HPP
#define PDE_POISSON_OBSTACLEK_HPP

#include "../TOOLS/pdeK.hpp"
#include "../TOOLS/feK.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real, class DeviceType>
class PDE_Obstacle : public PDE<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  // Finite element basis information
  basis_ptr basisPtr_;
  std::vector<basis_ptr> basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid2::Cubature<DeviceType,Real,Real>> cellCub_;
  // Cell node information
  scalar_view volCellNodes_;
  // Finite element definition
  ROL::Ptr<fe_type> fe_vol_;
  // Load term
  scalar_view load_;
  Real coeff_;
  // Use Riesz map for state variables?
  bool useStateRiesz_;

  Real evaluateLoad(const std::vector<Real> &coord) const {
    return static_cast<Real>(2)*coeff_;
  }

public:
  PDE_Obstacle(ROL::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("Problem").get("Basis Order",1);
    if (basisOrder == 1)
      basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType,Real,Real>>();
    else if (basisOrder == 2)
      basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<DeviceType,Real,Real>>();
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();                  // get the cell type from any basis
    Intrepid2::DefaultCubatureFactory cubFactory;                                      // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 2);              // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create<DeviceType,Real,Real>(cellType, cubDegree);           // create default cubature
    coeff_ = parlist.sublist("Problem").get("Load Magnitude", 1.0);                    // Pointwise load magnitude
    useStateRiesz_ = parlist.sublist("Problem").get("Use State Riesz Map", true);      // use Riesz map for state variables?
  }

  void residual(scalar_view & res,
                const scalar_view u_coeff,
                const scalar_view z_coeff = scalar_view(),
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // GET DIMENSIONS
    int c = fe_vol_->gradN().extent_int(0);
    int f = fe_vol_->gradN().extent_int(1);
    int p = fe_vol_->gradN().extent_int(2);
    int d = fe_vol_->gradN().extent_int(3);
    // INITIALIZE RESIDUAL
    res = scalar_view("res", c, f);
    // COMPUTE STIFFNESS TERM
    scalar_view gradU_eval("gradU_eval", c, p, d);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    fst::integrate(res,gradU_eval,fe_vol_->gradNdetJ(),false);
//    // ADD MASS TERM
//    scalar_view valU_eval("valU_eval", c, p);
//    fe_vol_->evaluateValue(valU_eval, u_coeff);
//    fst::integrate(res,valU_eval,fe_vol_->NdetJ(),true);
    // ADD LOAD TERM
    fst::integrate(res,load_,fe_vol_->NdetJ(),true);
  }

  void Jacobian_1(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // GET DIMENSIONS
    int c = fe_vol_->gradN().extent_int(0);
    int f = fe_vol_->gradN().extent_int(1);
    // INITIALIZE JACOBIAN
    jac = scalar_view("jac", c, f, f);
    // COMPUTE STIFFNESS TERM
    Kokkos::deep_copy(jac,fe_vol_->stiffMat());
//    // ADD MASS TERM
//    rst::add(jac,fe_vol_->massMat());
  }

  void Hessian_11(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Obstacle::Hessian_11): Hessian is zero.");
  }

  void RieszMap_1(scalar_view & riesz) override {
    // Optionally disable Riesz map ...
    if (!useStateRiesz_)
      throw Exception::NotImplemented(">>> (PDE_Obstacle::RieszMap_1): Not implemented.");

    // ...otherwise ...

    int c = fe_vol_->N().extent_int(0);
    int f = fe_vol_->N().extent_int(1);
    // INITIALIZE RIESZ MAP
    riesz = scalar_view("riesz1", c, f, f);
    Kokkos::deep_copy(riesz,fe_vol_->stiffMat());
    rst::add(riesz,fe_vol_->massMat());
  }

  std::vector<basis_ptr> getFields() override {
    return basisPtrs_;
  }

  void setCellNodes(const scalar_view &volCellNodes,
                    const std::vector<std::vector<scalar_view>> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds) override {
    volCellNodes_ = volCellNodes;
    // Finite element definition.
    fe_vol_ = ROL::makePtr<fe_type>(volCellNodes_,basisPtr_,cellCub_);
    // Compute load function
    computeLoad();
  }

  void computeLoad(void) {
    int c = fe_vol_->gradN().extent_int(0);
    int p = fe_vol_->gradN().extent_int(2);
    int d = fe_vol_->gradN().extent_int(3);
    std::vector<Real> coord(d);
    load_ = scalar_view("load_", c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k)
          coord[k] = (fe_vol_->cubPts())(i,j,k);
        load_(i,j) = -evaluateLoad(coord);
      }
    }
  }

  const ROL::Ptr<fe_type> getFE(void) const {
    return fe_vol_;
  }

  const scalar_view getCellNodes(void) const {
    return volCellNodes_;
  }

}; // PDE_Obstacle

#endif
