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

#ifndef PDE_HELMHOLTZK_HPP
#define PDE_HELMHOLTZK_HPP

#include "../TOOLS/pdeK.hpp"
#include "../TOOLS/feK.hpp"
#include "../TOOLS/fieldhelperK.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"

#include "ROL_Ptr.hpp"


template <class Real, class DeviceType>
class PDE_Helmholtz : public PDE<Real, DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type = FE<Real, DeviceType>;
  using fst = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct = Intrepid2::CellTools<DeviceType>;
  using rst = Intrepid2::RealSpaceTools<DeviceType>;

private:
  // Finite element basis information
  basis_ptr basisPtr_;
  std::vector<basis_ptr> basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid2::Cubature<DeviceType, Real, Real>> cellCub_;
  ROL::Ptr<Intrepid2::Cubature<DeviceType, Real, Real>> bdryCub_;
  // Cell node information
  scalar_view volCellNodes_;
  std::vector<std::vector<scalar_view>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<fe_type> fe_;
  std::vector<ROL::Ptr<fe_type>> feBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<scalar_view>> bdryCellDofValues_;
  // Field pattern, offsets, etc.
  std::vector<std::vector<int>> fieldPattern_;   // local Field/DOF pattern; set from DOF manager 
  int numFields_;                                // number of fields (equations in the PDE)
  int numDofs_;                                  // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                      // for each field, a counting offset
  std::vector<int> numFieldDofs_;                // for each field, number of degrees of freedom

  Real waveNumber_;
  Real RoiWaveNumber_;
  Real innerAnnulusRadius_;
  Real outerAnnulusRadius_;
  Real RoiRadius_;

  scalar_view ctrlWeight_;
  scalar_view ctrlJac_;

  ROL::Ptr<FieldHelper<Real, DeviceType>> fieldHelper_;

  void computeRefractiveIndex(scalar_view & kappa) const {
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);

    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = fe_->cubPts()(i, j, k);
        }
        kappa(i, j) = std::pow(evaluateRefractiveIndex(x), 2);
      }
    }
  }

  void computeForce(scalar_view & F, const int component) const {
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);

    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = fe_->cubPts()(i, j, k);
        }
        F(i,j) = evaluateForce(x, component);
      }
    }
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
        for (int k = 0; k < d; ++k) {
          x[k] = fe_->cubPts()(i, j, k);
        }
        if ( insideControlDomain(x) ) {
          inside = true;
          break;
        }
      }
      for (int j = 0; j < p; ++j) {
        ctrlWeight_(i,j) = (inside ? one : zero);
      }
    }
  }

  scalar_view getBoundaryCoeff(const scalar_view cell_coeff, int sideSet, int cell) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][cell];
    const int numCellsSide = bdryCellLocId.size();
    const int f = basisPtr_->getCardinality();

    scalar_view bdry_coeff("bdry_coeff", numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        bdry_coeff(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

  void buildControlJacobian(void) {
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);

    // Build force/control term
    scalar_view F("F", c, f, p);
    fst::scalarMultiplyDataField(F, ctrlWeight_, fe_->N());
    ctrlJac_ = scalar_view("ctrlJac_", c, f, f);
    fst::integrate(ctrlJac_, F, fe_->NdetJ(), false);
  }

public:
  PDE_Helmholtz(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType, Real, Real>>();
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();            // get the cell type from the basis
    Intrepid2::DefaultCubatureFactory cubFactory;                          // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 2);        // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create<DeviceType, Real, Real>(cellType, cubDegree);                           // create default cubature
    int d = cellType.getDimension();

    basisPtrs_.clear();
    for (int i=0; i<2; ++i) {
      basisPtrs_.push_back(basisPtr_);  // Displacement component
    }

    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2); // set cubature degree, e.g., 2
    bdryCub_ = cubFactory.create<DeviceType, Real, Real>(bdryCellType, bdryCubDegree);

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

    waveNumber_         = parlist.sublist("Problem").get("Wave Number",10.0);
    RoiWaveNumber_      = parlist.sublist("Problem").get("ROI Wave Number",10.0);

    Real dist2annulus   = parlist.sublist("Problem").get("Distance to Control Annulus",0.5);
    Real annulusWidth   = parlist.sublist("Problem").get("Control Annulus Width",0.1);
    RoiRadius_          = parlist.sublist("Problem").get("ROI Radius",2.0);
    innerAnnulusRadius_ = RoiRadius_ + dist2annulus;
    outerAnnulusRadius_ = innerAnnulusRadius_ + annulusWidth;
  }

  virtual Real evaluateRefractiveIndex(const std::vector<Real> &x) const {
    Real xnorm(0), val(0);
    const int d = x.size();
    for (int i = 0; i < d; ++i) {
      xnorm += x[i]*x[i];
    }
    xnorm = std::sqrt(xnorm);
    val   = (xnorm <= RoiRadius_) ? RoiWaveNumber_ : waveNumber_;
    return val;
  }

  virtual bool insideControlDomain(const std::vector<Real> &x) const {
    Real xnorm(0);
    const int d = x.size();
    for (int i = 0; i < d; ++i) {
      xnorm += x[i]*x[i];
    }
    xnorm = std::sqrt(xnorm);
    return (xnorm <= outerAnnulusRadius_ && xnorm >= innerAnnulusRadius_);
  }

  virtual Real evaluateForce(const std::vector<Real> &x, const int component) const {
    Real val(0);
    return val;
  }

  void residual(scalar_view & res,
                const scalar_view u_coeff,
                const scalar_view z_coeff = scalar_view(),
                const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);

    // Initialize residuals.
    std::vector<scalar_view> R(2);
    for (int i=0; i<2; ++i) {
      R[i] = scalar_view("R", c, f);
    }

    // Split u_coeff into components.
    std::vector<scalar_view> U;
    std::vector<scalar_view> Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<scalar_view> valZ_eval(2);
    std::vector<scalar_view> valU_eval(2);
    std::vector<scalar_view> gradU_eval(2);
    for (int i=0; i<2; ++i) {
      valZ_eval[i]  = scalar_view("valZ_eval", c, p);
      fe_->evaluateValue(valZ_eval[i], Z[i]);
      valU_eval[i]  = scalar_view("valU_eval", c, p);
      fe_->evaluateValue(valU_eval[i], U[i]);
      gradU_eval[i] = scalar_view("gradU_eval", c, p, d);
      fe_->evaluateGradient(gradU_eval[i], U[i]);
    }

    // Build force/control term
    std::vector<scalar_view> F(2);
    for (int i=0; i<2; ++i) {
      F[i] = scalar_view("F", c, p);
      computeForce(F[i], i);
      scalar_view wZ("wZ", c, p);
      fst::scalarMultiplyDataData(wZ, ctrlWeight_, valZ_eval[i]);
      rst::add(F[i], wZ);
    }

    // Build wave number
    scalar_view kappa("kappa", c, p);
    computeRefractiveIndex(kappa);
    std::vector<scalar_view> kappaU(2);
    for (int i=0; i<2; ++i) {
      kappaU[i] = scalar_view("kappaU", c, p);
      fst::scalarMultiplyDataData(kappaU[i], kappa, valU_eval[i]);
    }

    /*******************************************************************/
    /*** Evaluate weak form of the residual.****************************/
    /*** a(ur,vr) - a(ui,vi) + b(ur,vr) + b(ui,vi) = zr(vr) - zi(vi) ***/
    /*** a(ur,vi) + a(ui,vr) - b(ur,vi) + b(ui,vr) = zr(vi) + zi(vr) ***/
    /*******************************************************************/
    std::vector<scalar_view> stiff(2);
    std::vector<scalar_view> mass(2);
    std::vector<scalar_view> load(2);
    for (int i=0; i<2; ++i) {
      stiff[i] = scalar_view("stiff", c, f);
      mass[i]  = scalar_view("mass", c, f);
      load[i]  = scalar_view("load", c, f);

      fst::integrate(stiff[i],
                     gradU_eval[i],     // grad U
                     fe_->gradNdetJ(),  // grad N
                     false);
      fst::integrate(mass[i],
                     kappaU[i],         // -kappa2 U
                     fe_->NdetJ(),      // N
                     false);
      fst::integrate(load[i],
                     F[i],              // F
                     fe_->NdetJ(),      // N
                     false);
    }

    Kokkos::deep_copy(R[0], stiff[0]);
    rst::subtract(R[0], stiff[1]);
    rst::subtract(R[0], mass[0]);
    rst::add(R[0], mass[1]);
    rst::subtract(R[0], load[0]);
    rst::add(R[0], load[1]);

    Kokkos::deep_copy(R[1], stiff[0]);
    rst::add(R[1], stiff[1]);
    rst::subtract(R[1], mass[0]);
    rst::subtract(R[1], mass[1]);
    rst::subtract(R[1], load[0]);
    rst::subtract(R[1], load[1]);

    // APPLY ROBIN CONTROLS: Sideset 0
    int sideset = 0;
    int numLocalSideIds = bdryCellLocIds_[sideset].size();
    const int numCubPerSide = bdryCub_->getNumPoints();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[sideset][j].size();
      if (numCellsSide) {
        std::vector<scalar_view> robinRes(2);
        for (int i = 0; i < 2; ++i) {
          robinRes[i] = scalar_view("robinRes", numCellsSide, f);
          // Get U coefficients on Robin boundary
          scalar_view u_coeff_bdry = getBoundaryCoeff(U[i], sideset, j);
          // Evaluate U on FE basis
          scalar_view valU_eval_bdry("valU_eval_bdry", numCellsSide, numCubPerSide);
          feBdry_[j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
          // Compute Neumann residual
          fst::integrate(robinRes[i], valU_eval_bdry, feBdry_[j]->NdetJ(), false);
        }
        // Add Neumann control residual to volume residual
        for (int k = 0; k < numCellsSide; ++k) {
          int cidx = bdryCellLocIds_[sideset][j][k];
          for (int l = 0; l < f; ++l) {
            (R[0])(cidx, l) += waveNumber_ * ((robinRes[0])(k, l) + (robinRes[1])(k, l));
            (R[1])(cidx, l) += waveNumber_ * ((robinRes[1])(k, l) - (robinRes[0])(k, l));
          }
        }
      }
    }

    // Combine the residuals.
    fieldHelper_->combineFieldCoeff(res, R);
  }

  void Jacobian_1(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);

    // Initialize Jacobians.
    std::vector<std::vector<scalar_view>> J(2);
    for (int i=0; i<2; ++i) {
      for (int j=0; j<2; ++j) {
        J[i].push_back(scalar_view("J1", c, f, f));
      }
    }

    // Build wave number
    scalar_view kappa("kappa", c, p);
    computeRefractiveIndex(kappa);
    scalar_view kappaN("kappaN", c, f, p);
    fst::scalarMultiplyDataField(kappaN, kappa, fe_->N());

    /*** Evaluate weak form of the Jacobian. ***/
    scalar_view kappaM("kappaM", c, f, f);
    fst::integrate(kappaM, kappaN, fe_->NdetJ(), false);

    Kokkos::deep_copy(J[0][0], fe_->stiffMat());
    rst::subtract(J[0][0], kappaM);

    Kokkos::deep_copy(J[0][1], fe_->stiffMat());
    rst::subtract(J[0][1], kappaM);
    rst::scale(J[0][1], static_cast<Real>(-1));

    Kokkos::deep_copy(J[1][0], fe_->stiffMat());
    rst::subtract(J[1][0], kappaM);

    Kokkos::deep_copy(J[1][1], fe_->stiffMat());
    rst::subtract(J[1][1], kappaM);

    // APPLY ROBIN CONTROL: Sideset 0
    int sideset = 0;
    int numLocalSideIds = bdryCellLocIds_[sideset].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[sideset][j].size();
      if (numCellsSide) {
        // Add Neumann control Jacobian to volume residual
        for (int k = 0; k < numCellsSide; ++k) {
          int cidx = bdryCellLocIds_[sideset][j][k];
          for (int l = 0; l < f; ++l) {
            for (int m = 0; m < f; ++m) {
              (J[0][0])(cidx, l, m) += waveNumber_ * (feBdry_[j]->massMat())(k, l, m);
              (J[0][1])(cidx, l, m) += waveNumber_ * (feBdry_[j]->massMat())(k, l, m);
              (J[1][0])(cidx, l, m) -= waveNumber_ * (feBdry_[j]->massMat())(k, l, m);
              (J[1][1])(cidx, l, m) += waveNumber_ * (feBdry_[j]->massMat())(k, l, m);
            }
          }
        }
      }
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }


  void Jacobian_2(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);

    // Initialize Jacobians.
    std::vector<std::vector<scalar_view>> J(2);
    for (int i=0; i<2; ++i) {
      for (int j=0; j<2; ++j) {
        J[i].push_back(scalar_view("J2", c, f, f));
      }
    }

    /*** Evaluate weak form of the Jacobian. ***/
    Kokkos::deep_copy(J[0][0], ctrlJac_);
    rst::scale(J[0][0], static_cast<Real>(-1));

    Kokkos::deep_copy(J[0][1], ctrlJac_);

    Kokkos::deep_copy(J[1][0], ctrlJac_);
    rst::scale(J[1][0], static_cast<Real>(-1));

    Kokkos::deep_copy(J[1][1], ctrlJac_);
    rst::scale(J[1][1], static_cast<Real>(-1));

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }

  void Hessian_11(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Helmholtz::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Helmholtz::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Helmholtz::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Helmholtz::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(scalar_view & riesz) override {
    //throw Exception::NotImplemented(">>> (PDE_TopoOpt::RieszMap_1): Not implemented.");

    // Retrieve dimensions.
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int d = fe_->gradN().extent_int(3);

    // Initialize Jacobians.
    std::vector<std::vector<scalar_view>> J(d);
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j) {
        J[i].push_back(scalar_view("JRM1", c, f, f));
      }
    }

    for (int i=0; i<d; ++i) {
      Kokkos::deep_copy(J[i][i], fe_->stiffMat());
      rst::add(J[i][i], fe_->massMat());
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  void RieszMap_2(scalar_view & riesz) override {
    //throw Exception::NotImplemented(">>> (PDE_TopoOpt::RieszMap_2): Not implemented.");

    // Retrieve dimensions.
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int d = fe_->gradN().extent_int(3);
 
    // Initialize Jacobians.
    std::vector<std::vector<scalar_view>> J(d);
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j) {
        J[i].push_back(scalar_view("JRM2", c, f, f));
      }
    }

    for (int i=0; i<d; ++i) {
      Kokkos::deep_copy(J[i][i], fe_->massMat());
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  std::vector<basis_ptr> getFields() override {
    return basisPtrs_;
  }

  void setCellNodes(const scalar_view & volCellNodes,
                    const std::vector<std::vector<scalar_view>> & bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> & bdryCellLocIds) override {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fe_ = ROL::makePtr<fe_type>(volCellNodes_, basisPtr_, cellCub_);
    fidx_ = fe_->getBoundaryDofs();
    // Construct boundary FE
    int sideset = 0;
    int numLocSides = bdryCellNodes[sideset].size();
    feBdry_.resize(numLocSides);
    for (int j = 0; j < numLocSides; ++j) {
      if (bdryCellNodes[sideset][j] != scalar_view()) {
        feBdry_[j] = ROL::makePtr<fe_type>(bdryCellNodes[sideset][j], basisPtr_, bdryCub_, j);
      }
    }
    // Compute control weight
    computeControlWeight();
    buildControlJacobian();
  }

  void setFieldPattern(const std::vector<std::vector<int>> & fieldPattern) override {
    fieldPattern_ = fieldPattern;
    fieldHelper_ = ROL::makePtr<FieldHelper<Real, DeviceType>>(numFields_, numDofs_, numFieldDofs_, fieldPattern_);
  }

  const ROL::Ptr<fe_type> getFE(void) const {
    return fe_;
  }

  const std::vector<fe_type> getBdryFE(void) const {
    return feBdry_;
  }

  const std::vector<std::vector<int>> getBdryCellLocIds(const int sideset = 0) const {
    return bdryCellLocIds_[sideset];
  }

  const ROL::Ptr<FieldHelper<Real, DeviceType>> getFieldHelper(void) const {
    return fieldHelper_;
  }

}; // PDE_Helmholtz


#endif
