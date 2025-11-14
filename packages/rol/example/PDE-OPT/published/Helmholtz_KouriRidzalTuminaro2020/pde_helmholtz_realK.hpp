// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde_helmholtz.hpp
    \brief Implements the local PDE interface for the optimal control of
           Helmholtz.
*/

#ifndef PDE_HELMHOLTZ_REALK_HPP
#define PDE_HELMHOLTZ_REALK_HPP

#include "../../TOOLS/pdeK.hpp"
#include "../../TOOLS/feK.hpp"
#include "../../TOOLS/fieldhelperK.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real, class DeviceType>
class PDE_Helmholtz_Real : public PDE<Real,DeviceType> {
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
  std::vector<std::vector<scalar_view>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<fe_type> fe_;

  Real innerAnnulusRadius_;
  Real outerAnnulusRadius_;
  int example_;
  Real waveNumber_;

  scalar_view ctrlWeight_;
  scalar_view ctrlJac_;

  bool insideControlDomain(const std::vector<Real> &x) const {
    bool val = true;
    if (example_==1) {
      const Real eps = std::sqrt(ROL::ROL_EPSILON<Real>());
      Real xnorm(0);
      const int d = x.size();
      for (int i = 0; i < d; ++i) {
        xnorm += x[i]*x[i];
      }
      xnorm = std::sqrt(xnorm);
      val = (xnorm <= outerAnnulusRadius_+eps && xnorm >= innerAnnulusRadius_-eps);
    }
    return val;
  }

  void computeControlWeight(void) {
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);
   
    ctrlWeight_ = scalar_view("ctrlWeight_", c, p);

    const Real zero(0), one(1);
    bool inside(false);
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      inside = false;
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k)
          x[k] = (fe_->cubPts())(i,j,k);
        if ( insideControlDomain(x) ) {
          inside = true;
          break;
        }
      }
      for (int j = 0; j < p; ++j)
        ctrlWeight_(i,j) = (inside ? one : zero);
    }
  }

  void buildControlJacobian(void) {
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);

    // Build force/control term
    scalar_view F("F", c, f, p);
    fst::scalarMultiplyDataField(F, ctrlWeight_, fe_->N());
    ctrlJac_ = scalar_view("ctrlJac", c, f, f);
    fst::integrate(ctrlJac_,F,fe_->NdetJ(),false);
    rst::scale(ctrlJac_,static_cast<Real>(-1));
  }

public:
  PDE_Helmholtz_Real(ROL::ParameterList &parlist)
    : innerAnnulusRadius_(2.5), outerAnnulusRadius_(2.6),
      example_(parlist.sublist("Problem").get("Example",1)),
      waveNumber_(parlist.sublist("Problem").get("Wave Number",10.0)) {
    // Finite element fields.
    basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType,Real,Real>>();
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();        // get the cell type from the basis
    Intrepid2::DefaultCubatureFactory cubFactory;                            // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 2);    // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create<DeviceType,Real,Real>(cellType, cubDegree); // create default cubature
  }

  void residual(scalar_view & res,
                const scalar_view u_coeff,
                const scalar_view z_coeff = scalar_view(),
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    const Real one(1), two(2);
    // Retrieve dimensions.
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);
    // Initialize residuals.
    res = scalar_view("res",c,f);
    // Evaluate/interpolate finite element fields on cells.
    scalar_view valz_eval("valz_eval", c, p);
    scalar_view valu_eval("valu_eval", c, p);
    scalar_view gradu_eval("gradu_eval", c, p, d);
    fe_->evaluateValue(valz_eval, z_coeff);
    fe_->evaluateValue(valu_eval, u_coeff);
    fe_->evaluateGradient(gradu_eval, u_coeff);
    // Integrate PDE term
    fst::integrate(res,valu_eval,fe_->NdetJ(),false);
    rst::scale(res, -std::pow(waveNumber_,two));
    fst::integrate(res,gradu_eval,fe_->gradNdetJ(),true);
    // Build control term
    scalar_view F("F", c, p);
    fst::scalarMultiplyDataData(F, ctrlWeight_, valz_eval);
    rst::scale(F,-one);
    // Volumetric intregration
    fst::integrate(res,F,fe_->NdetJ(),true);
  }

  void Jacobian_1(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    const Real two(2);
    // Retrieve dimensions.
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    // Initialize Jacobian.
    jac = scalar_view("jac",c,f,f);
    // Add PDE terms
    Kokkos::deep_copy(jac,fe_->massMat());
    rst::scale(jac, -std::pow(waveNumber_,two));
    rst::add(jac, fe_->stiffMat());
  }


  void Jacobian_2(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    // Initialize Jacobian.
    jac = scalar_view("jac",c,f,f);
    // Add control term
    Kokkos::deep_copy(jac,ctrlJac_);
    rst::scale(jac, static_cast<Real>(-1));
  }

  void Hessian_11(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Helmholtz::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Helmholtz::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Helmholtz::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Helmholtz::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(scalar_view & riesz) override {
    // Retrieve dimensions.
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    // Initialize Jacobian.
    riesz = scalar_view("riesz1",c,f,f);
    Kokkos::deep_copy(riesz,fe_->stiffMat());
    rst::add(riesz, fe_->massMat());
  }

  void RieszMap_2(scalar_view & riesz) override {
    // Retrieve dimensions.
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    // Initialize Jacobian.
    riesz = scalar_view("riesz2",c,f,f);
    Kokkos::deep_copy(riesz,fe_->massMat());
  }

  std::vector<basis_ptr> getFields() override {
    return basisPtrs_;
  }

  void setCellNodes(const scalar_view &volCellNodes,
                    const std::vector<std::vector<scalar_view>> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds) override {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fe_ = ROL::makePtr<fe_type>(volCellNodes_,basisPtr_,cellCub_);
    // Compute control weight
    computeControlWeight();
    buildControlJacobian();
  }

  const ROL::Ptr<fe_type> getFE(void) const {
    return fe_;
  }

}; // PDE_Helmholtz


#endif
