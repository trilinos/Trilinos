// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde_elasticity.hpp
    \brief Implements the local PDE interface for the structural topology
           optimization problem.
*/

#ifndef PDE_TOPO_OPT_ELASTICITY_HPP
#define PDE_TOPO_OPT_ELASTICITY_HPP

#include "../../../TOOLS/pdeK.hpp"
#include "../../../TOOLS/feK.hpp"
#include "../../../TOOLS/fieldhelperK.hpp"

#include "dirichletK.hpp"
#include "tractionK.hpp"
#include "loadK.hpp"
#include "materialtensorK.hpp"

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
class PDE_Elasticity : public PDE<Real, DeviceType> {
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
  ROL::Ptr<Intrepid2::Cubature<DeviceType,Real,Real>> cellCub_, bdryCub_;
  // Cell node information
  scalar_view volCellNodes_;
  std::vector<std::vector<scalar_view>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<fe_type> fe_, feDens_;
  std::vector<std::vector<ROL::Ptr<fe_type>>> feBdry_;
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
  std::vector<std::vector<int>> fieldPatternDens_;   // local Field/DOF pattern; set from DOF manager 
  int numFieldsDens_;                                // number of fields (equations in the PDE)
  int numDofsDens_;                                  // total number of degrees of freedom for all (local) fields
  std::vector<int> offsetDens_;                      // for each field, a counting offset
  std::vector<int> numFieldDofsDens_;                // for each field, number of degrees of freedom

  ROL::Ptr<Load<Real,DeviceType>>            load_; 
  ROL::Ptr<MaterialTensor<Real,DeviceType>>  matTensor_;
  ROL::Ptr<Dirichlet<Real,DeviceType>>       dirichlet_;
  ROL::Ptr<Traction<Real,DeviceType>>        traction_;
  ROL::Ptr<FieldUtils::FieldInfo> fieldInfo_, fieldInfoDens_;

  bool getFields2called_;

public:
  PDE_Elasticity(ROL::ParameterList &parlist) : getFields2called_(false) {
    // Finite element fields.
    int basisOrder     = parlist.sublist("Problem").get("Basis Order",1);
    int basisOrderDens = parlist.sublist("Problem").get("Density Basis Order",0);
    int cubDegree      = parlist.sublist("Problem").get("Cubature Degree",4);
    int bdryCubDegree  = parlist.sublist("Problem").get("Boundary Cubature Degree",2);
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
      if (basisOrder == 1)
        basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType,Real,Real>>();
      else if (basisOrder == 2)
        basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<DeviceType,Real,Real>>();
    }
    else if (probDim == 3) {
      if (basisOrder == 1)
        basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_HEX_C1_FEM<DeviceType,Real,Real>>();
      else if (basisOrder == 2)
        basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_HEX_C2_FEM<DeviceType,Real,Real>>();
    }
    basisPtrs_.clear(); basisPtrs_.resize(probDim,basisPtr_);
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology(); // get cell type from basis
    if (basisOrderDens == 1)
      basisPtrDens_ = ROL::makePtr<Intrepid2::Basis_HVOL_C0_FEM<DeviceType,Real,Real>>(cellType);
    else {
      if (probDim == 2)
        basisPtrDens_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType,Real,Real>>();
      else if (probDim == 3)
        basisPtrDens_ = ROL::makePtr<Intrepid2::Basis_HGRAD_HEX_C1_FEM<DeviceType,Real,Real>>();
    }
    basisPtrsDens_.clear(); basisPtrsDens_.push_back(basisPtrDens_); // Density component

    // Quadrature rules.
    Intrepid2::DefaultCubatureFactory cubFactory;                            // create cubature factory
    cellCub_ = cubFactory.create<DeviceType,Real,Real>(cellType, cubDegree); // create default cubature

    shards::CellTopology bdryCellType = cellType.getCellTopologyData(probDim-1, 0);
    bdryCub_ = cubFactory.create<DeviceType,Real,Real>(bdryCellType, bdryCubDegree);

    matTensor_ = ROL::makePtr<MaterialTensor<Real,DeviceType>>(parlist.sublist("Problem"));
    std::string example = parlist.sublist("Problem").get("Example","Default");
    load_    = ROL::makePtr<Load<Real,DeviceType>>(parlist.sublist("Problem"),example);
    traction_= ROL::makePtr<Traction<Real,DeviceType>>(parlist.sublist("Problem"),example);
    dirichlet_ = ROL::makePtr<Dirichlet<Real,DeviceType>>(parlist.sublist("Problem"),example);

    numDofs_ = 0;
    numFields_ = basisPtrs_.size();
    offset_.resize(numFields_);
    numFieldDofs_.resize(numFields_);
    for (int i=0; i<numFields_; ++i) {
      if (i==0)
        offset_[i]  = 0;
      else
        offset_[i]  = offset_[i-1] + basisPtrs_[i-1]->getCardinality();
      numFieldDofs_[i] = basisPtrs_[i]->getCardinality();
      numDofs_ += numFieldDofs_[i];
    }

    numDofsDens_ = 0;
    numFieldsDens_ = basisPtrsDens_.size();
    offsetDens_.resize(numFieldsDens_);
    numFieldDofsDens_.resize(numFieldsDens_);
    for (int i=0; i<numFieldsDens_; ++i) {
      if (i==0)
        offsetDens_[i]  = 0;
      else
        offsetDens_[i]  = offsetDens_[i-1] + basisPtrsDens_[i-1]->getCardinality();
      numFieldDofsDens_[i] = basisPtrsDens_[i]->getCardinality();
      numDofsDens_ += numFieldDofsDens_[i];
    }
  }

  void residual(scalar_view &res,
                const scalar_view u_coeff,
                const scalar_view z_coeff = scalar_view(),
                const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);
    int matd = matTensor_->getMatrixDim();
 
    // Initialize residuals.
    std::vector<scalar_view> R(d);
    for (int i=0; i<d; ++i) R[i] = scalar_view("res",c,f);

    // Split u_coeff into components.
    std::vector<scalar_view> U, Z;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoDens_);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<scalar_view> gradDisp_eval(d);
    scalar_view UMat;
    scalar_view rho("rho", c, p);
    scalar_view valZ_eval("valZ_eval", c, p);
    scalar_view rhoUMat("rhoUMat", c, p, matd);
    for (int i=0; i<d; ++i) {
      gradDisp_eval[i] =  scalar_view("gradDisp_eval", c, p, d);
      fe_->evaluateGradient(gradDisp_eval[i], U[i]);
    }
    feDens_->evaluateValue(valZ_eval, Z[0]);

    // EVALUATE MATERIAL TENSOR
    matTensor_->computeUmat(UMat, gradDisp_eval);
    matTensor_->computeDensity(rho, valZ_eval);
    fst::scalarMultiplyDataData(rhoUMat, rho, UMat);
    for (int i=0; i<d; ++i)
      fst::integrate(R[i],rhoUMat,matTensor_->CBdetJ(i),false);

    // EVALUATE LOAD
    if (!load_->isNull()) {
      std::vector<scalar_view> load(d);
      for (int i=0; i<d; ++i)
        load[i] = scalar_view("load", c, p);
      load_->compute(load, fe_, PDE<Real,DeviceType>::getParameter(), static_cast<Real>(-1));
      for (int i=0; i<d; ++i)
        fst::integrate(R[i],load[i],fe_->NdetJ(),true);
    }

    // APPLY TRACTION CONDITIONS
    if (!traction_->isNull()) {
      traction_->apply(R, feBdry_, PDE<Real,DeviceType>::getParameter(), static_cast<Real>(-1));
    }

    // APPLY DIRICHLET CONDITIONS
    dirichlet_->applyResidual(R,U);

    // Combine the residuals.
    FieldUtils::combineFieldCoeff<Real>(res, R, fieldInfo_);
  }

  void Jacobian_1(scalar_view &jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);
    int matd = matTensor_->getMatrixDim();
 
    // Initialize Jacobians.
    std::vector<std::vector<scalar_view>> J(d);
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j)
        J[i].push_back(scalar_view("jac1",c,f,f));
    }

    // Split z_coeff into components.
    std::vector<scalar_view> Z;
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoDens_);

    std::vector<scalar_view> rhoBMat(d);
    scalar_view rho("rho", c, p);
    scalar_view valZ_eval("valZ_eval", c, p);

    feDens_->evaluateValue(valZ_eval, Z[0]);
    matTensor_->computeDensity(rho, valZ_eval);
    for (int i=0; i<d; ++i) {
      rhoBMat[i] =  scalar_view("rhoBMat", c, f, p, matd);
      fst::scalarMultiplyDataField(rhoBMat[i], rho, matTensor_->B(i));
    }

    /*** Evaluate weak form of the Jacobian. ***/
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j)
        fst::integrate(J[i][j],rhoBMat[i],matTensor_->CBdetJ(j),false);
    }

    // APPLY DIRICHLET CONDITIONS
    dirichlet_->applyJacobian1(J);

    // Combine the jacobians.
    FieldUtils::combineFieldCoeff<Real>(jac, J, fieldInfo_, fieldInfo_);
  }

  void Jacobian_2(scalar_view &jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    int c  = fe_->gradN().extent_int(0);
    int f  = fe_->gradN().extent_int(1);
    int fd = feDens_->gradN().extent_int(1);
    int p  = fe_->gradN().extent_int(2);
    int d  = fe_->gradN().extent_int(3);
    int matd = matTensor_->getMatrixDim();

    // Initialize Jacobians.
    std::vector<std::vector<scalar_view>> J(d);
    for (int i=0; i<d; ++i)
      J[i].push_back(scalar_view("jac1",c,f,fd));

    // Split u_coeff into components.
    std::vector<scalar_view> U, Z;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoDens_);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<scalar_view> gradDisp_eval(d), CBrhoUMat(d);
    scalar_view UMat;
    scalar_view rho("rho", c, p);
    scalar_view valZ_eval("valZ_eval", c, p);
    scalar_view rhoUMat("rhoUMat", c, p, matd);
    for (int i=0; i<d; ++i) {
      CBrhoUMat[i]     = scalar_view("CBrhoUMat", c, f, p);
      gradDisp_eval[i] = scalar_view("gradDisp_eval", c, p, d);
      fe_->evaluateGradient(gradDisp_eval[i], U[i]);
    }
    feDens_->evaluateValue(valZ_eval, Z[0]);
    matTensor_->computeUmat(UMat, gradDisp_eval);
    matTensor_->computeDensity(rho, valZ_eval, 1);  // first derivative
    fst::scalarMultiplyDataData(rhoUMat, rho, UMat);

    for (int i=0; i<d; ++i)
      fst::dotMultiplyDataField(CBrhoUMat[i], rhoUMat, matTensor_->CBdetJ(i));

    /*** Evaluate weak form of the residual. ***/
    for (int i=0; i<d; ++i)
      fst::integrate(J[i][0],CBrhoUMat[i],feDens_->N(),false);

    // APPLY DIRICHLET CONDITIONS
    dirichlet_->applyJacobian2(J);

    // Combine the jacobians.
    FieldUtils::combineFieldCoeff<Real>(jac, J, fieldInfo_, fieldInfoDens_);
  }

  void Hessian_11(scalar_view &hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_TopoOpt::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(scalar_view &hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    int c  = fe_->gradN().extent_int(0);
    int f  = fe_->gradN().extent_int(1);
    int fd = feDens_->gradN().extent_int(1);
    int p  = fe_->gradN().extent_int(2);
    int d  = fe_->gradN().extent_int(3);
    int matd = matTensor_->getMatrixDim();

    // Initialize Hessians.
    std::vector<std::vector<scalar_view>> J(1);
    for (int i=0; i<d; ++i)
      J[0].push_back(scalar_view("hes12",c,fd,f));

    // Split u_coeff into components.
    std::vector<scalar_view> L, Z;
    FieldUtils::splitFieldCoeff<Real>(L, l_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoDens_);

    // Apply Dirichlet conditions to the multipliers.
    dirichlet_->applyMultiplier(L);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<scalar_view> gradDisp_eval(d), CBrhoLMat(d);
    scalar_view LMat;
    scalar_view rho("rho", c, p);
    scalar_view valZ_eval("valZ_eval", c, p);
    scalar_view rhoLMat("rhoLMat", c, p, matd);
    for (int i=0; i<d; ++i) {
      CBrhoLMat[i]     = scalar_view("CBrhoLMat", c, f, p);
      gradDisp_eval[i] = scalar_view("gradDisp_eval", c, p, d);
      fe_->evaluateGradient(gradDisp_eval[i], L[i]);
    }
    feDens_->evaluateValue(valZ_eval, Z[0]);
    matTensor_->computeUmat(LMat, gradDisp_eval);
    matTensor_->computeDensity(rho, valZ_eval, 1);  // first derivative
    fst::scalarMultiplyDataData(rhoLMat, rho, LMat);

    for (int i=0; i<d; ++i)
      fst::dotMultiplyDataField(CBrhoLMat[i], rhoLMat, matTensor_->CBdetJ(i));

    /*** Evaluate weak form of the residual. ***/
    for (int i=0; i<d; ++i)
      fst::integrate(J[0][i],feDens_->N(),CBrhoLMat[i],false);

    // Combine the Hessians.
    FieldUtils::combineFieldCoeff<Real>(hess, J, fieldInfoDens_, fieldInfo_);

  }

  void Hessian_21(scalar_view &hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    int c  = fe_->gradN().extent_int(0);
    int f  = fe_->gradN().extent_int(1);
    int fd = feDens_->gradN().extent_int(1);
    int p  = fe_->gradN().extent_int(2);
    int d  = fe_->gradN().extent_int(3);
    int matd = matTensor_->getMatrixDim();

    // Initialize Hessians.
    std::vector<std::vector<scalar_view>> J(d);
    for (int i=0; i<d; ++i)
      J[i].push_back(scalar_view("hes21",c,f,fd));

    // Split u_coeff into components.
    std::vector<scalar_view> L, Z;
    FieldUtils::splitFieldCoeff<Real>(L, l_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoDens_);

    // Apply Dirichlet conditions to the multipliers.
    dirichlet_->applyMultiplier(L);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<scalar_view> gradDisp_eval(d), CBrhoLMat(d);
    scalar_view LMat;
    scalar_view rho("rho", c, p);
    scalar_view valZ_eval("valZ_eval", c, p);
    scalar_view rhoLMat("rhoLMat", c, p, matd);
    for (int i=0; i<d; ++i) {
      CBrhoLMat[i]     = scalar_view("CBrhoLMat", c, f, p);
      gradDisp_eval[i] = scalar_view("gradDisp_eval", c, p, d);
      fe_->evaluateGradient(gradDisp_eval[i], L[i]);
    }
    feDens_->evaluateValue(valZ_eval, Z[0]);
    matTensor_->computeUmat(LMat, gradDisp_eval);
    matTensor_->computeDensity(rho, valZ_eval, 1);  // first derivative
    fst::scalarMultiplyDataData(rhoLMat, rho, LMat);

    for (int i=0; i<d; ++i)
      fst::dotMultiplyDataField(CBrhoLMat[i], rhoLMat, matTensor_->CBdetJ(i));

    /*** Evaluate weak form of the residual. ***/
    for (int i=0; i<d; ++i)
      fst::integrate(J[i][0],CBrhoLMat[i],feDens_->N(),false);

    // Combine the Hessians.
    FieldUtils::combineFieldCoeff(hess, J, fieldInfo_, fieldInfoDens_);
  }

  void Hessian_22(scalar_view &hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    int c  = fe_->gradN().extent_int(0);
    int fd = feDens_->gradN().extent_int(1);
    int p  = fe_->gradN().extent_int(2);
    int d  = fe_->gradN().extent_int(3);
    int matd = matTensor_->getMatrixDim();

    // Initialize Hessians.
    std::vector<std::vector<scalar_view>> J(1);
    J[0].push_back(scalar_view("hes22",c,fd,fd));

    // Split u_coeff into components.
    std::vector<scalar_view> U, L, Z;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(L, l_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoDens_);

    // Apply Dirichlet conditions to the multipliers.
    dirichlet_->applyMultiplier(L);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<scalar_view> gradDispU_eval(d), gradDispL_eval(d);
    scalar_view UMat, LMat;
    scalar_view rho("rho", c, p);
    scalar_view valZ_eval("valZ_eval", c, p);
    scalar_view rhoLMat("rhoLMat", c, p, matd);
    scalar_view CUMat("CUMat", c, p, matd);
    scalar_view CUrhoLMat("CUrhoLMat", c, p);
    scalar_view NCUrhoLMat("NCUrhoLMat", c, fd, p);
    for (int i=0; i<d; ++i) {
      gradDispU_eval[i] = scalar_view("gradDispU_eval", c, p, d);
      fe_->evaluateGradient(gradDispU_eval[i], U[i]);
      gradDispL_eval[i] = scalar_view("gradDispL_eval", c, p, d);
      fe_->evaluateGradient(gradDispL_eval[i], L[i]);
    }
    feDens_->evaluateValue(valZ_eval, Z[0]);
    matTensor_->computeUmat(UMat, gradDispU_eval);
    matTensor_->computeUmat(LMat, gradDispL_eval);
    matTensor_->computeDensity(rho, valZ_eval, 2);  // second derivative
    fst::scalarMultiplyDataData(rhoLMat, rho, LMat);
    matTensor_->applyTensor(CUMat, UMat);

    fst::dotMultiplyDataData(CUrhoLMat, rhoLMat, CUMat);
    fst::scalarMultiplyDataField(NCUrhoLMat, CUrhoLMat, feDens_->N());

    /*** Evaluate weak form of the residual. ***/
    fst::integrate(J[0][0],NCUrhoLMat,feDens_->NdetJ(),false);

    // Combine the Hessians.
    FieldUtils::combineFieldCoeff(hess, J, fieldInfoDens_, fieldInfoDens_);
  }

  void RieszMap_1(scalar_view & riesz) override {
    // Retrieve dimensions.
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int d = fe_->gradN().extent_int(3);
 
    // Initialize Jacobians.
    std::vector<std::vector<scalar_view>> J(d);
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j)
        J[i].push_back(scalar_view("riesz1",c,f,f));
    }

    for (int i=0; i<d; ++i) {
      Kokkos::deep_copy(J[i][i],fe_->stiffMat());
      rst::add(J[i][i],fe_->massMat());
    }

    // Combine the jacobians.
    FieldUtils::combineFieldCoeff<Real>(riesz, J, fieldInfo_, fieldInfo_);
  }

  void RieszMap_2(scalar_view & riesz) override {
    // Retrieve dimensions.
    int c = feDens_->gradN().extent_int(0);
    int f = feDens_->gradN().extent_int(1);

    riesz = scalar_view("riesz2",c,f,f);
    Kokkos::deep_copy(riesz,feDens_->massMat());
  }

  // This must be called before getFields2
  void setDensityFields(const std::vector<basis_ptr> &basisPtrs) {
    if (getFields2called_) {
      TEUCHOS_TEST_FOR_EXCEPTION(getFields2called_, std::invalid_argument,
        ">>> PDE-OPT/topo-opt/elasticity/src/pde_elasticity.hpp: Must call before getFields2!");
    }
    else {
      basisPtrDens_  = basisPtrs[0];
      basisPtrsDens_ = basisPtrs;

      numDofsDens_ = 0;
      numFieldsDens_ = basisPtrsDens_.size();
      offsetDens_.resize(numFieldsDens_);
      numFieldDofsDens_.resize(numFieldsDens_);
      for (int i=0; i<numFieldsDens_; ++i) {
        if (i==0)
          offsetDens_[i]  = 0;
        else
          offsetDens_[i]  = offsetDens_[i-1] + basisPtrsDens_[i-1]->getCardinality();
        numFieldDofsDens_[i] = basisPtrsDens_[i]->getCardinality();
        numDofsDens_ += numFieldDofsDens_[i];
      }
    }
  }

  std::vector<basis_ptr> getFields() override {
    return basisPtrs_;
  }

  std::vector<basis_ptr> getFields2() override {
    getFields2called_ = true;
    return basisPtrsDens_;
  }

  void setCellNodes(const scalar_view &volCellNodes,
                    const std::vector<std::vector<scalar_view>> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds) override {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fe_ = ROL::makePtr<fe_type>(volCellNodes_,basisPtr_,cellCub_);
    feDens_ = ROL::makePtr<fe_type>(volCellNodes_,basisPtrDens_,cellCub_,false);
    fidx_ = fe_->getBoundaryDofs();
    if (!traction_->isNull()) {
      traction_->setCellNodes(bdryCellNodes_,bdryCellLocIds_);
    }
    dirichlet_->setCellNodes(bdryCellNodes_,bdryCellLocIds_,fidx_);
    matTensor_->setFE(fe_);
    // Construct boundary FE
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      feBdry_.resize(numSideSets);
      for (int i = 0; i < numSideSets; ++i) {
        int numLocSides = bdryCellNodes[i].size();
        feBdry_[i].resize(numLocSides);
        for (int j = 0; j < numLocSides; ++j) {
          if (bdryCellNodes[i][j] != scalar_view()) {
            feBdry_[i][j] = ROL::makePtr<fe_type>(bdryCellNodes[i][j],basisPtr_,bdryCub_,j);
          }
        }
      }
    }
  }

  void setFieldPattern(const std::vector<std::vector<int>> &fieldPattern,
                       const std::vector<std::vector<int>> &fieldPattern2) override {
    fieldPattern_     = fieldPattern;
    fieldPatternDens_ = fieldPattern2;
    fieldInfo_     = ROL::makePtr<FieldUtils::FieldInfo>(numFields_, numDofs_, numFieldDofs_, fieldPattern_);
    fieldInfoDens_ = ROL::makePtr<FieldUtils::FieldInfo>(numFieldsDens_, numDofsDens_, numFieldDofsDens_, fieldPatternDens_);
  }

  const ROL::Ptr<fe_type> getStateFE(void) const {
    return fe_;
  }

  const ROL::Ptr<fe_type> getDensityFE(void) const {
    return feDens_;
  }

  const std::vector<std::vector<ROL::Ptr<fe_type>>> getBdryFE(void) const {
    return feBdry_;
  }

  const std::vector<std::vector<std::vector<int>>> getBdryCellLocIds(void) const {
    return bdryCellLocIds_;
  }

  const ROL::Ptr<FieldUtils::FieldInfo> getStateFieldInfo(void) const {
    return fieldInfo_;
  }

  const ROL::Ptr<FieldUtils::FieldInfo> getDensityFieldInfo(void) const {
    return fieldInfoDens_;
  }

  const ROL::Ptr<Load<Real,DeviceType>> getLoad(void) const {
    return load_;
  }

  const ROL::Ptr<Traction<Real,DeviceType>> getTraction(void) const {
    return traction_;
  }

  const ROL::Ptr<MaterialTensor<Real,DeviceType>> getMaterialTensor(void) const {
    return matTensor_;
  }

}; // PDE_TopoOpt

#endif
