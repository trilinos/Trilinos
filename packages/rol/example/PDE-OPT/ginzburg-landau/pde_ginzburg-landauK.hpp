// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde_ginzburg-landau.hpp
    \brief Implements the local PDE interface for the optimal control of
           simplified Ginzburg-Landau.
*/

#ifndef PDE_GINZBURGLANDAUK_HPP
#define PDE_GINZBURGLANDAUK_HPP

#include "../TOOLS/pdeK.hpp"
#include "../TOOLS/feK.hpp"
#include "../TOOLS/fieldhelperK.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"

#include "ROL_Ptr.hpp"


template <class Real, class DeviceType>
class PDE_GinzburgLandau : public PDE<Real, DeviceType> {
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
  ROL::Ptr<Intrepid2::Cubature<DeviceType,Real,Real>> cellCub_, bdryCub_;
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
  std::vector<std::vector<int>> fieldPattern_;  // local Field/DOF pattern; set from DOF manager 
  int numFields_;                               // number of fields (equations in the PDE)
  int numDofs_;                                 // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                     // for each field, a counting offset
  std::vector<int> numFieldDofs_;               // for each field, number of degrees of freedom

  Real lambda_;

  ROL::Ptr<FieldHelper<Real,DeviceType>> fieldHelper_;

  void computeMagneticPotential(scalar_view &A) const {
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);
   
    std::vector<Real> x(d), Ax(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (fe_->cubPts())(i,j,k);
        }
        evaluateMagneticPotential(Ax,x);
        for (int k = 0; k < d; ++k) {
          A(i,j,k) = Ax[k];
        }
      }
    } 
  }

  void computeForce(scalar_view &F, int component) const {
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);
   
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (fe_->cubPts())(i,j,k);
        }
        F(i,j) = -evaluateForce(x,component);
      }
    } 
  }

  void computeNeumann(scalar_view &neumann,
                      int locSideId,                              
                      int component) const {
    const int c = feBdry_[locSideId]->gradN().extent_int(0);
    const int p = feBdry_[locSideId]->gradN().extent_int(2);
    const int d = fe_->gradN().extent_int(3);
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (feBdry_[locSideId]->cubPts())(i,j,k);
        }
        neumann(i,j) = evaluateNeumann(x,component)/lambda_;
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

public:
  PDE_GinzburgLandau(ROL::ParameterList &parlist) {
    // Finite element fields.
    basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType,Real,Real>>();
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();        // get the cell type from the basis
    Intrepid2::DefaultCubatureFactory cubFactory;                            // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 2);    // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create<DeviceType,Real,Real>(cellType, cubDegree); // create default cubature
    int d = cellType.getDimension();

    basisPtrs_.clear();
    basisPtrs_.resize(2,basisPtr_); // Displacement component

    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2); // set cubature degree, e.g., 2
    bdryCub_ = cubFactory.create<DeviceType,Real,Real>(bdryCellType, bdryCubDegree);

    numDofs_ = 0;
    numFields_ = basisPtrs_.size();
    offset_.resize(numFields_);
    numFieldDofs_.resize(numFields_);
    for (int i=0; i<numFields_; ++i) {
      if (i==0) offset_[i]  = 0;
      else      offset_[i]  = offset_[i-1] + basisPtrs_[i-1]->getCardinality();
      numFieldDofs_[i] = basisPtrs_[i]->getCardinality();
      numDofs_ += numFieldDofs_[i];
    }

    lambda_ = parlist.sublist("Problem").get("Current Loading",1.0);
  }

  virtual void evaluateMagneticPotential(std::vector<Real> &Ax, const std::vector<Real> &x) const = 0;

  virtual Real evaluateNeumann(const std::vector<Real> &x, const int component) const = 0;

  virtual Real evaluateForce(const std::vector<Real> &x, const int component) const = 0;

  void residual(scalar_view &res,
                const scalar_view u_coeff,
                const scalar_view z_coeff = scalar_view(),
                const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);
    const int p = fe_->gradN().extent_int(2);
    const int d = fe_->gradN().extent_int(3);
 
    // Initialize residuals.
    std::vector<scalar_view> R(2);
    R[0] = scalar_view("res", c,f);
    R[1] = scalar_view("res", c,f);

    // Split u_coeff into components.
    std::vector<scalar_view> U, Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    std::vector<scalar_view> valU_eval(2), gradU_eval(2), F(2);
    scalar_view A("A", c, p, d);
    scalar_view magA("magA", c, p);
    scalar_view magU("magU", c, p);
    scalar_view sqrU1("sqrU1", c, p);
    scalar_view magAU("magAU", c, p);
    scalar_view AgradU("AgradU", c, p);
    scalar_view AU("AU", c, p, d);
    scalar_view magUU("magUU", c, p);
    for (int i=0; i<2; ++i) {
      valU_eval[i]  = scalar_view("valU_eval", c, p);
      gradU_eval[i] = scalar_view("gradU_eval", c, p, d);
      F[i]          = scalar_view("F", c, p);
       // Evaluate/interpolate finite element fields on cells.
      fe_->evaluateValue(valU_eval[i], U[i]);
      fe_->evaluateGradient(gradU_eval[i], U[i]);
      // Build force term
      computeForce(F[i],i);
    }
    // Build magnetic potential
    computeMagneticPotential(A);
    // Compute the magnitude of A
    rst::dot(magA,A,A);
    // Compute magnitude of U
    fst::scalarMultiplyDataData(magU, valU_eval[0],valU_eval[0]);
    fst::scalarMultiplyDataData(sqrU1,valU_eval[1],valU_eval[1]);
    rst::add(magU,sqrU1);
    rst::scale(magU,lambda_);

    /*** Evaluate weak form of the residual. ***/
    int index = 0;
    for (int i = 0; i < 2; ++i) {
      index = (i+1)%2;
      // Symmetric term
      fst::integrate(R[i],gradU_eval[i],fe_->gradNdetJ(),false);
      Kokkos::deep_copy(magAU,static_cast<Real>(0));
      fst::scalarMultiplyDataData(magAU,magA,valU_eval[i]);
      rst::subtract(magAU,valU_eval[i]);
      fst::integrate(R[i],magAU,fe_->NdetJ(),true);
      // Nonlinear term
      Kokkos::deep_copy(magUU,static_cast<Real>(0));
      fst::scalarMultiplyDataData(magUU,magU,valU_eval[i]);
      fst::integrate(R[i],magUU,fe_->NdetJ(),true);
      // Antisymmetric term
      Kokkos::deep_copy(AgradU,static_cast<Real>(0));
      fst::dotMultiplyDataData(AgradU,A,gradU_eval[index]);
      if (index==1) rst::scale(AgradU,static_cast<Real>(-1));
      fst::integrate(R[i],AgradU,fe_->NdetJ(),true);
      Kokkos::deep_copy(AU,static_cast<Real>(0));
      fst::scalarMultiplyDataData(AU,valU_eval[index],A);
      if (index==0) rst::scale(AU,static_cast<Real>(-1));
      fst::integrate(R[i],AU,fe_->gradNdetJ(),true);
      // Force term
      fst::integrate(R[i],F[i],fe_->NdetJ(),true);
    }

    // APPLY ROBIN CONTROLS: Sideset 0
    int sideset = 0;
    int numLocalSideIds = bdryCellLocIds_[sideset].size();
    const int numCubPerSide = bdryCub_->getNumPoints();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[sideset][j].size();
      if (numCellsSide) {
        for (int i = 0; i < 2; ++i) {
          // Get U coefficients on Robin boundary
          scalar_view z_coeff_bdry = getBoundaryCoeff(Z[i], sideset, j);
          // Evaluate U on FE basis
          scalar_view valZ_eval_bdry("valZ_eval_bdry", numCellsSide, numCubPerSide);
          feBdry_[j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
          // Compute uncontrolled Neumann source
          scalar_view neumann("neumann", numCellsSide, numCubPerSide);
          computeNeumann(neumann,j,i);
          // Add uncontrolled Neumann source to control
          rst::add(valZ_eval_bdry, neumann);
          // Compute Neumann residual
          scalar_view neumannRes("neumannRes", numCellsSide, f);
          fst::integrate(neumannRes,valZ_eval_bdry,feBdry_[j]->NdetJ(),false);
          // Add Neumann control residual to volume residual
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[sideset][j][k];
            for (int l = 0; l < f; ++l) { 
              (R[i])(cidx,l) -= lambda_*neumannRes(k,l);
            }
          }
        }
      }
    }

    // Combine the residuals.
    fieldHelper_->combineFieldCoeff(res, R);
  }

  void Jacobian_1(scalar_view &jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);
    const int p = fe_->gradN().extent_int(2);
    const int d = fe_->gradN().extent_int(3);
 
    // Initialize Jacobians.
    std::vector<std::vector<scalar_view>> J(2);
    for (int i=0; i<2; ++i) {
      for (int j=0; j<2; ++j) {
        J[i].push_back(scalar_view("jac",c,f,f));
      }
    }

    // Split u_coeff into components.
    std::vector<scalar_view> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<scalar_view> valU_eval(2), sqrUN(2);
    scalar_view A("A", c, p, d);
    scalar_view magA("magA", c, p);
    scalar_view magAN("magAN", c, f, p);
    scalar_view AgradN("AgradN", c, f, p);
    scalar_view sqrU("sqrU", c, p);
    scalar_view U0U1("U0U1", c, p);
    scalar_view U0U1N("U0U1N", c, f, p);
    scalar_view sqrUN3("sqrUN3", c, f, p);
    for (int i=0; i<2; ++i) {
      valU_eval[i] = scalar_view("valU_eval", c, p);
      fe_->evaluateValue(valU_eval[i], U[i]);
      // Compute jacobian of nonlinearity
      Kokkos::deep_copy(sqrU, static_cast<Real>(0));
      fst::scalarMultiplyDataData(sqrU,valU_eval[i],valU_eval[i]);
      rst::scale(sqrU,lambda_);
      sqrUN[i] = scalar_view("sqrUN", c, f, p);
      fst::scalarMultiplyDataField(sqrUN[i],sqrU,fe_->N());
    }
    // Build magnetic potential
    computeMagneticPotential(A);
    // Compute magnitude of magnetic potential
    rst::dot(magA,A,A);
    // Multiply magnitude of magnetic potential with basis function
    fst::scalarMultiplyDataField(magAN,magA,fe_->N());
    rst::subtract(magAN,fe_->N());
    // Dot magnetic potential with gradient of basis function
    fst::dotMultiplyDataField(AgradN,A,fe_->gradN());
    // Compute jacobian of nonlinearity
    fst::scalarMultiplyDataData(U0U1,valU_eval[0],valU_eval[1]);
    rst::scale(U0U1,static_cast<Real>(2)*lambda_);
    fst::scalarMultiplyDataField(U0U1N,U0U1,fe_->N());

    /*** Evaluate weak form of the Jacobian. ***/
    for (int i=0; i<2; ++i) {
      int index = (i+1)%2;
      // Diagonal linear terms
      Kokkos::deep_copy(J[i][i],fe_->stiffMat());
      fst::integrate(J[i][i],magAN,fe_->NdetJ(),true);
      // Diagonal nonlinear terms
      Kokkos::deep_copy(sqrUN3, static_cast<Real>(0));
      rst::scale(sqrUN3,sqrUN[i],static_cast<Real>(3));
      fst::integrate(J[i][i],sqrUN3,fe_->NdetJ(),true);
      fst::integrate(J[i][i],sqrUN[index],fe_->NdetJ(),true);
      // Off-diagonal nonlinear terms
      fst::integrate(J[i][index],U0U1N,fe_->NdetJ(),false);
      // Off-diagonal linear terms
      fst::integrate(J[i][index],AgradN,fe_->NdetJ(),true);
      rst::scale(AgradN,static_cast<Real>(-1));
      fst::integrate(J[i][index],fe_->NdetJ(),AgradN,true);
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }


  void Jacobian_2(scalar_view &jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);

    // Initialize Jacobians.
    std::vector<std::vector<scalar_view>> J(2);
    for (int i=0; i<2; ++i) {
      for (int j=0; j<2; ++j) {
        J[i].push_back(scalar_view("jac",c,f,f));
      }
    }

    // APPLY ROBIN CONTROL: Sideset 0
    int sideset = 0;
    int numLocalSideIds = bdryCellLocIds_[sideset].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[sideset][j].size();
      if (numCellsSide) {
        for (int i = 0; i < 2; ++i) {
          // Add Neumann control Jacobian to volume residual
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[sideset][j][k];
            for (int l = 0; l < f; ++l) { 
              for (int m = 0; m < f; ++m) { 
                (J[i][i])(cidx,l,m) -= lambda_*(feBdry_[j]->massMat())(k,l,m);
              }
            }
          }
        }
      }
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);

  }

  void Hessian_11(scalar_view &hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);
    const int p = fe_->gradN().extent_int(2);
 
    // Initialize Jacobians.
    std::vector<std::vector<scalar_view>> H(2);
    for (int i=0; i<2; ++i) {
      for (int j=0; j<2; ++j) {
        H[i].push_back(scalar_view("hess",c,f,f));
      }
    }

    // Split u_coeff and l_coeff into components.
    std::vector<scalar_view> U, L;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(L, l_coeff);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<scalar_view> valU_eval(2), valL_eval(2), UL(2);
    scalar_view U0L1_U1L0("U0L1_U1L0",c,p);
    scalar_view U1L0("U1L0",c,p);
    scalar_view U0L1_U1L0_N("U0L1_U1L0_N",c,f,p);
    scalar_view diag("diag",c,p);
    scalar_view diag_N("diag_N",c,f,p);
    for (int i=0; i<2; ++i) {
      valU_eval[i] = scalar_view("valU_eval", c, p);
      fe_->evaluateValue(valU_eval[i], U[i]);
      valL_eval[i] = scalar_view("valL_eval", c, p);
      fe_->evaluateValue(valL_eval[i], L[i]);
      UL[i] = scalar_view("UL", c, p);
      fst::scalarMultiplyDataData(UL[i],valU_eval[i],valL_eval[i]);
    }
    fst::scalarMultiplyDataData(U0L1_U1L0,valU_eval[0],valL_eval[1]);
    fst::scalarMultiplyDataData(U1L0,     valU_eval[1],valL_eval[0]);
    rst::add(U0L1_U1L0,U1L0);
    rst::scale(U0L1_U1L0,static_cast<Real>(2)*lambda_);
    fst::scalarMultiplyDataField(U0L1_U1L0_N,U0L1_U1L0,fe_->N());

    for (int i = 0; i < 2; ++i) {
      int index = (i+1)%2;
      Kokkos::deep_copy(diag, static_cast<Real>(0));
      Kokkos::deep_copy(diag_N, static_cast<Real>(0));
      rst::scale(diag,UL[i],static_cast<Real>(3));
      rst::add(diag,UL[index]);
      rst::scale(diag,static_cast<Real>(2)*lambda_);
      fst::scalarMultiplyDataField(diag_N,diag,fe_->N());
      fst::integrate(H[i][i],diag_N,fe_->NdetJ(),false);
      fst::integrate(H[i][index],U0L1_U1L0_N,fe_->NdetJ(),false);
    }

    // Combine the hessians.
    fieldHelper_->combineFieldCoeff(hess, H);
  }

  void Hessian_12(scalar_view &hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_GinzburgLandau::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(scalar_view &hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_GinzburgLandau::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(scalar_view &hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_GinzburgLandau::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(scalar_view & riesz) override {
    //throw Exception::NotImplemented(">>> (PDE_TopoOpt::RieszMap_1): Not implemented.");

    // Retrieve dimensions.
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);
 
    // Initialize Jacobians.
    std::vector<std::vector<scalar_view>> J(2);
    for (int i=0; i<2; ++i) {
      for (int j=0; j<2; ++j) {
        J[i].push_back(scalar_view("riesz",c,f,f));
      }
    }

    for (int i=0; i<2; ++i) {
      Kokkos::deep_copy(J[i][i],fe_->stiffMat());
      rst::add(J[i][i],fe_->massMat());
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  void RieszMap_2(scalar_view & riesz) override {
    //throw Exception::NotImplemented(">>> (PDE_TopoOpt::RieszMap_2): Not implemented.");

    // Retrieve dimensions.
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);
 
    // Initialize Jacobians.
    std::vector<std::vector<scalar_view>> J(2);
    for (int i=0; i<2; ++i) {
      for (int j=0; j<2; ++j) {
        J[i].push_back(scalar_view("riesz",c,f,f));
      }
    }

    for (int i=0; i<2; ++i) {
      Kokkos::deep_copy(J[i][i],fe_->massMat());
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
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
    fidx_ = fe_->getBoundaryDofs();
    // Construct boundary FE
    int sideset = 0;
    int numLocSides = bdryCellNodes[sideset].size();
    feBdry_.resize(numLocSides);
    for (int j = 0; j < numLocSides; ++j) {
      if (bdryCellNodes[sideset][j] != scalar_view()) {
        feBdry_[j] = ROL::makePtr<fe_type>(bdryCellNodes[sideset][j],basisPtr_,bdryCub_,j);
      }
    }
  }

  void setFieldPattern(const std::vector<std::vector<int>> &fieldPattern) override {
    fieldPattern_ = fieldPattern;
    fieldHelper_ = ROL::makePtr<FieldHelper<Real,DeviceType>>(numFields_, numDofs_, numFieldDofs_, fieldPattern_);
  }

  const ROL::Ptr<fe_type> getFE(void) const {
    return fe_;
  }

  const std::vector<ROL::Ptr<fe_type>> getBdryFE(void) const {
    return feBdry_;
  }

  const std::vector<std::vector<int>> getBdryCellLocIds(const int sideset = 0) const {
    return bdryCellLocIds_[sideset];
  }

  const ROL::Ptr<FieldHelper<Real,DeviceType>> getFieldHelper(void) const {
    return fieldHelper_;
  }

}; // PDE_GinzburgLandau


#endif
