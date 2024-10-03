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

#ifndef PDE_TOPO_OPT_MULTIMAT_HPP
#define PDE_TOPO_OPT_MULTIMAT_HPP

#include "../../../TOOLS/pde.hpp"
#include "../../../TOOLS/fe.hpp"
#include "../../../TOOLS/fieldhelper.hpp"
#include "../../../TOOLS/Intrepid_HGRAD_C0_FEM.hpp"

#include "dirichlet.hpp"
#include "traction.hpp"
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
class PDE_MultiMat_Elasticity : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>> basisPtr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrs_;
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>> basisPtrDens_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrsDens_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real>> cellCub_;
  ROL::Ptr<Intrepid::Cubature<Real>> bdryCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real>> volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real>> fe_, feDens_;
  std::vector<std::vector<ROL::Ptr<FE<Real>>>> feBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellDofValues_;
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

  ROL::Ptr<MultiMat_Load<Real>>                        load_; 
  std::vector<ROL::Ptr<MultiMat_MaterialTensor<Real>>> matTensor_;
  ROL::Ptr<MultiMat_Dirichlet<Real>>                   dirichlet_;
  ROL::Ptr<MultiMat_Traction<Real>>                    traction_;
  ROL::Ptr<FieldUtils::FieldInfo>             fieldInfo_, fieldInfoDens_;

  bool getFields2called_;

public:
  PDE_MultiMat_Elasticity(Teuchos::ParameterList &parlist) : getFields2called_(false) {
    // Material properties
    ROL::ParameterList list;
    list.set("Use Plain Stress", parlist.sublist("Problem").get("Use Plain Stress", true));
    list.set("SIMP Power", parlist.sublist("Problem").get("SIMP Power", 1.0));
    std::vector<Real> ym = ROL::getArrayFromStringParameter<Real>(parlist.sublist("Problem"), "Young's Modulus");
    std::vector<Real> pr = ROL::getArrayFromStringParameter<Real>(parlist.sublist("Problem"), "Poisson Ratio");
    int T = ym.size();
    matTensor_.resize(T);
    for (int t=0; t<T; ++t) {
      list.set("Young's Modulus",ym[t]);
      list.set("Poisson Ratio",pr[t]);
      matTensor_[t] = ROL::makePtr<MultiMat_MaterialTensor<Real>>(list);
    }
    std::string example = parlist.sublist("Problem").get("Example","Default");
    load_      = ROL::makePtr<MultiMat_Load<Real>>(parlist.sublist("Problem"),example);
    traction_  = ROL::makePtr<MultiMat_Traction<Real>>(parlist.sublist("Problem"),example);
    dirichlet_ = ROL::makePtr<MultiMat_Dirichlet<Real>>(parlist.sublist("Problem"),example);

    // Finite element fields.
    int basisOrder     = parlist.sublist("Problem").get("Basis Order",1);
    int basisOrderDens = parlist.sublist("Problem").get("Density Basis Order",0);
    int cubDegree      = parlist.sublist("Problem").get("Cubature Degree",4);
    int bdryCubDegree  = parlist.sublist("Problem").get("Boundary Cubature Degree",2);
    int probDim        = parlist.sublist("Problem").get("Problem Dimension",2);
    if (probDim> 3 || probDim < 2) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> PDE-OPT/binary/elasticity/pde_elasticity.hpp: Problem dimension is not 2 or 3!");
    }
    if (basisOrder> 2 || basisOrder < 1) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> PDE-OPT/binary/elasticity/pde_elasticity.hpp: Basis order is not 1 or 2!");
    }
    if (basisOrderDens > 1 || basisOrderDens < 0) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> PDE-OPT/binary/elasticity/pde_elasticity.hpp: Density basis order is not 0 1!");
    }
    if (probDim == 2) {
      if (basisOrder == 1) {
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
      else if (basisOrder == 2) {
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
      if (basisOrderDens == 0) {
        basisPtrDens_ = ROL::makePtr<Intrepid::Basis_HGRAD_C0_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
      else if (basisOrderDens == 1) {
        basisPtrDens_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
    }
    else if (probDim == 3) {
      if (basisOrder == 1) {
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_HEX_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
      else if (basisOrder == 2) {
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_HEX_C2_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
      if (basisOrderDens == 0) {
        // This does not work in 3D
        basisPtrDens_ = ROL::makePtr<Intrepid::Basis_HGRAD_C0_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
      else if (basisOrderDens == 1) {
        basisPtrDens_ = ROL::makePtr<Intrepid::Basis_HGRAD_HEX_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
    }
    basisPtrs_.clear(); basisPtrsDens_.clear();
    for (int i=0; i<probDim; ++i) {
      basisPtrs_.push_back(basisPtr_);  // Displacement component
    }
    for (int i=0; i<T; ++i) {
      basisPtrsDens_.push_back(basisPtrDens_); // Density components
    }

    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology(); // get cell type from basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                // create cubature factory
    cellCub_ = cubFactory.create(cellType, cubDegree);                // create default cubature

    shards::CellTopology bdryCellType = cellType.getCellTopologyData(probDim-1, 0);
    bdryCub_ = cubFactory.create(bdryCellType, bdryCubDegree);

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
    int d = fe_->gradN()->dimension(3);
    int matd = matTensor_[0]->getMatrixDim();
    int T    = matTensor_.size(); 
 
    // Initialize residuals.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> R(d);
    for (int i=0; i<d; ++i) {
      R[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f);
    }

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U, Z;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoDens_);

    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> UMat, rho, valZ_eval, rhoUMat;
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> gradDisp_eval(d);
    rho       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    rhoUMat   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
    for (int i=0; i<d; ++i) {
      gradDisp_eval[i] =  ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      fe_->evaluateGradient(gradDisp_eval[i], U[i]);
    }
    matTensor_[0]->computeUmat(UMat, gradDisp_eval);
    for (int t=0; t<T; ++t) {
      valZ_eval->initialize(); rho->initialize(); rhoUMat->initialize();
      feDens_->evaluateValue(valZ_eval, Z[t]);
      matTensor_[t]->computeDensity(rho, valZ_eval);
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*rhoUMat, *rho, *UMat);
      for (int i=0; i<d; ++i) {
        Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                      *rhoUMat,                  // C B U
                                                      *matTensor_[t]->CBdetJ(i), // B'
                                                      Intrepid::COMP_CPP,
                                                      true);
      }
    }

    // EVALUATE LOAD
    if (!load_->isNull()) {
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> load(d);
      for (int i=0; i<d; ++i) {
        load[i] =  ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      }
      load_->compute(load, fe_, PDE<Real>::getParameter(), static_cast<Real>(-1));
      for (int i=0; i<d; ++i) {
        Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                      *load[i],           // F
                                                      *fe_->NdetJ(),      // N
                                                      Intrepid::COMP_CPP,
                                                      true);
      }
    }

    // APPLY TRACTION CONDITIONS
    if (!traction_->isNull()) {
      traction_->apply(R, feBdry_, PDE<Real>::getParameter(), static_cast<Real>(-1));
    }

    // APPLY DIRICHLET CONDITIONS
    dirichlet_->applyResidual(R, U);

    // Combine the residuals.
    FieldUtils::combineFieldCoeff<Real>(res, R, fieldInfo_);
  }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    int matd = matTensor_[0]->getMatrixDim();
    int T    = matTensor_.size();
 
    // Initialize Jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J(d);
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j) {
        J[i].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f));
      }
    }

    // Split z_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> Z;
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoDens_);

    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> rho, valZ_eval, rhoBMat;
    rho       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    rhoBMat   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p, matd);
    for (int t=0; t<T; ++t) {
      valZ_eval->initialize(); rho->initialize();
      feDens_->evaluateValue(valZ_eval, Z[t]);
      matTensor_[t]->computeDensity(rho, valZ_eval);

      /*** Evaluate weak form of the Jacobian. ***/
      for (int i=0; i<d; ++i) {
        rhoBMat->initialize();
        Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*rhoBMat, *rho, *matTensor_[t]->B(i));
        for (int j=0; j<d; ++j) {
          Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][j],
                                                        *rhoBMat,                  // rho B
                                                        *matTensor_[t]->CBdetJ(j), // B' C
                                                        Intrepid::COMP_CPP,
                                                        true);
        }
      }
    }

    // APPLY DIRICHLET CONDITIONS
    dirichlet_->applyJacobian1(J);

    // Combine the jacobians.
    FieldUtils::combineFieldCoeff<Real>(jac, J, fieldInfo_, fieldInfo_);
  }


  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c  = fe_->gradN()->dimension(0);
    int f  = fe_->gradN()->dimension(1);
    int fd = feDens_->gradN()->dimension(1);
    int p  = fe_->gradN()->dimension(2);
    int d  = fe_->gradN()->dimension(3);
    int matd = matTensor_[0]->getMatrixDim();
    int T    = matTensor_.size();

    // Initialize Jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J(d);
    for (int i=0; i<d; ++i) {
      for (int j=0; j<T; ++j) {
        J[i].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,fd));
      }
    }

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U, Z;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoDens_);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> gradDisp_eval(d);
    ROL::Ptr<Intrepid::FieldContainer<Real>> rho, valZ_eval, UMat, rhoUMat, CBrhoUMat;
    rho       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    rhoUMat   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
    CBrhoUMat = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    for (int i=0; i<d; ++i) {
      gradDisp_eval[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      fe_->evaluateGradient(gradDisp_eval[i], U[i]);
    }
    matTensor_[0]->computeUmat(UMat, gradDisp_eval);
    for (int t=0; t<T; ++t) {
      valZ_eval->initialize(); rho->initialize(); rhoUMat->initialize();
      feDens_->evaluateValue(valZ_eval, Z[t]);
      matTensor_[t]->computeDensity(rho, valZ_eval, 1);  // first derivative
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*rhoUMat, *rho, *UMat);

      for (int i=0; i<d; ++i) {
        CBrhoUMat->initialize();
        Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*CBrhoUMat, *rhoUMat, *matTensor_[t]->CBdetJ(i));
        /*** Evaluate weak form of the residual. ***/
        Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][t],
                                                      *CBrhoUMat,        // B' C drho B U
                                                      *feDens_->N(),     // N
                                                      Intrepid::COMP_CPP,
                                                      false);
      }
    }

    // APPLY DIRICHLET CONDITIONS
    dirichlet_->applyJacobian2(J);

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
    // Retrieve dimensions.
    int c  = fe_->gradN()->dimension(0);
    int f  = fe_->gradN()->dimension(1);
    int fd = feDens_->gradN()->dimension(1);
    int p  = fe_->gradN()->dimension(2);
    int d  = fe_->gradN()->dimension(3);
    int matd = matTensor_[0]->getMatrixDim();
    int T    = matTensor_.size();

    // Initialize Hessians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J(T);
    for (int t=0; t<T; ++t) {
      for (int i=0; i<d; ++i) {
        J[t].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,fd,f));
      }
    }

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> L, Z;
    FieldUtils::splitFieldCoeff<Real>(L, l_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoDens_);

    // Apply Dirichlet conditions to the multipliers.
    dirichlet_->applyMultiplier(L);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> gradDisp_eval(d);
    ROL::Ptr<Intrepid::FieldContainer<Real>> rho, valZ_eval, LMat, rhoLMat, CBrhoLMat;
    rho       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    rhoLMat   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
    CBrhoLMat = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    for (int i=0; i<d; ++i) {
      gradDisp_eval[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      fe_->evaluateGradient(gradDisp_eval[i], L[i]);
    }
    matTensor_[0]->computeUmat(LMat, gradDisp_eval);
    for (int t=0; t<T; ++t) {
      valZ_eval->initialize(); rho->initialize(); rhoLMat->initialize();
      feDens_->evaluateValue(valZ_eval, Z[t]);
      matTensor_[t]->computeDensity(rho, valZ_eval, 1);  // first derivative
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*rhoLMat, *rho, *LMat);

      for (int i=0; i<d; ++i) {
        CBrhoLMat->initialize();
        Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*CBrhoLMat, *rhoLMat, *matTensor_[t]->CBdetJ(i));
        /*** Evaluate weak form of the residual. ***/
        Intrepid::FunctionSpaceTools::integrate<Real>(*J[t][i],
                                                      *feDens_->N(), // N
                                                      *CBrhoLMat,    // B' C drho B U
                                                      Intrepid::COMP_CPP,
                                                      false);
      }
    }

    // Combine the Hessians.
    FieldUtils::combineFieldCoeff<Real>(hess, J, fieldInfoDens_, fieldInfo_);
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c  = fe_->gradN()->dimension(0);
    int f  = fe_->gradN()->dimension(1);
    int fd = feDens_->gradN()->dimension(1);
    int p  = fe_->gradN()->dimension(2);
    int d  = fe_->gradN()->dimension(3);
    int matd = matTensor_[0]->getMatrixDim();
    int T    = matTensor_.size();

    // Initialize Hessians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J(d);
    for (int i=0; i<d; ++i) {
      for (int t=0; t<T; ++t) {
        J[i].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,fd));
      }
    }

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> L, Z;
    FieldUtils::splitFieldCoeff<Real>(L, l_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoDens_);

    // Apply Dirichlet conditions to the multipliers.
    dirichlet_->applyMultiplier(L);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> gradDisp_eval(d);
    ROL::Ptr<Intrepid::FieldContainer<Real>> rho, valZ_eval, LMat, rhoLMat, CBrhoLMat;
    rho       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    rhoLMat   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
    CBrhoLMat = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    for (int i=0; i<d; ++i) {
      gradDisp_eval[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      fe_->evaluateGradient(gradDisp_eval[i], L[i]);
    }
    matTensor_[0]->computeUmat(LMat, gradDisp_eval);
    for (int t=0; t<T; ++t) {
      valZ_eval->initialize(); rho->initialize(); rhoLMat->initialize();
      feDens_->evaluateValue(valZ_eval, Z[t]);
      matTensor_[t]->computeDensity(rho, valZ_eval, 1);  // first derivative
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*rhoLMat, *rho, *LMat);

      for (int i=0; i<d; ++i) {
        CBrhoLMat->initialize();
        Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*CBrhoLMat, *rhoLMat, *matTensor_[t]->CBdetJ(i));
        /*** Evaluate weak form of the residual. ***/
        Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][t],
                                                      *CBrhoLMat,    // B' C drho B U
                                                      *feDens_->N(), // N
                                                      Intrepid::COMP_CPP,
                                                      false);
      }
    }

    // Combine the Hessians.
    FieldUtils::combineFieldCoeff(hess, J, fieldInfo_, fieldInfoDens_);
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c  = fe_->gradN()->dimension(0);
    int fd = feDens_->gradN()->dimension(1);
    int p  = fe_->gradN()->dimension(2);
    int d  = fe_->gradN()->dimension(3);
    int matd = matTensor_[0]->getMatrixDim();
    int T    = matTensor_.size();

    // Initialize Hessians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J(T);
    for (int t=0; t<T; ++t) {
      for (int s=0; s<T; ++s) {
        J[t].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,fd,fd));
      }
    }

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U, L, Z;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(L, l_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoDens_);

    // Apply Dirichlet conditions to the multipliers.
    dirichlet_->applyMultiplier(L);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> gradDispU_eval(d), gradDispL_eval(d);
    ROL::Ptr<Intrepid::FieldContainer<Real>> rho, valZ_eval, UMat, LMat, rhoLMat, CUMat, CUrhoLMat, NCUrhoLMat;
    rho        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    valZ_eval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    rhoLMat    = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
    CUMat      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
    CUrhoLMat  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    NCUrhoLMat = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fd, p);
    for (int i=0; i<d; ++i) {
      gradDispU_eval[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      fe_->evaluateGradient(gradDispU_eval[i], U[i]);
      gradDispL_eval[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      fe_->evaluateGradient(gradDispL_eval[i], L[i]);
    }
    matTensor_[0]->computeUmat(UMat, gradDispU_eval);
    matTensor_[0]->computeUmat(LMat, gradDispL_eval);
    for (int t=0; t<T; ++t) {
      valZ_eval->initialize(); rho->initialize(); rhoLMat->initialize();
      CUMat->initialize(); CUrhoLMat->initialize(); NCUrhoLMat->initialize();
      feDens_->evaluateValue(valZ_eval, Z[t]);
      matTensor_[t]->computeDensity(rho, valZ_eval, 2);  // second derivative
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*rhoLMat, *rho, *LMat);
      matTensor_[t]->applyTensor(CUMat, UMat);

      Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(*CUrhoLMat, *rhoLMat, *CUMat);
      Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*NCUrhoLMat, *CUrhoLMat, *feDens_->N());

      /*** Evaluate weak form of the residual. ***/
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[t][t],
                                                    *NCUrhoLMat,        // N L' B' C ddrho B U
                                                    *feDens_->NdetJ(),  // N
                                                    Intrepid::COMP_CPP,
                                                    false);
    }

    // Combine the Hessians.
    FieldUtils::combineFieldCoeff(hess, J, fieldInfoDens_, fieldInfoDens_);
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    //throw Exception::NotImplemented(">>> (PDE_TopoOpt::RieszMap_1): Not implemented.");

    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int d = fe_->gradN()->dimension(3);
 
    // Initialize Riesz maps.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> R(d);
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j) {
        R[i].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f));
      }
    }

    for (int i=0; i<d; ++i) {
      *(R[i][i]) = *(fe_->stiffMat());
      Intrepid::RealSpaceTools<Real>::add(*(R[i][i]),*(fe_->massMat()));
    }

    // Combine the Riesz maps.
    FieldUtils::combineFieldCoeff<Real>(riesz, R, fieldInfo_, fieldInfo_);
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    //throw Exception::NotImplemented(">>> (PDE_TopoOpt::RieszMap_2): Not implemented.");

    // Retrieve dimensions.
    int c = feDens_->gradN()->dimension(0);
    int f = feDens_->gradN()->dimension(1);
    int T = matTensor_.size();
 
    // Initialize Riesz maps.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> R(numFieldsDens_);
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

  // This must be called before getFields2
  void setDensityFields(const std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> &basisPtrs) {
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
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields() {
    return basisPtrs_;
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields2() {
    getFields2called_ = true;
    return basisPtrsDens_;
  }

  void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real>> &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds) {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fe_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtr_,cellCub_);
    feDens_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrDens_,cellCub_,false);
    fidx_ = fe_->getBoundaryDofs();
    if (!traction_->isNull()) {
      traction_->setCellNodes(bdryCellNodes_,bdryCellLocIds_);
    }
    dirichlet_->setCellNodes(bdryCellNodes_,bdryCellLocIds_,fidx_);
    int T = matTensor_.size();
    for (int t=0; t<T; ++t) {
      matTensor_[t]->setFE(fe_);
    }
    // Construct boundary FE
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets> 0) {
      feBdry_.resize(numSideSets);
      for (int i = 0; i < numSideSets; ++i) {
        int numLocSides = bdryCellNodes[i].size();
        feBdry_[i].resize(numLocSides);
        for (int j = 0; j < numLocSides; ++j) {
          if (bdryCellNodes[i][j] != ROL::nullPtr) {
            feBdry_[i][j] = ROL::makePtr<FE<Real>>(bdryCellNodes[i][j],basisPtr_,bdryCub_,j);
          }
        }
      }
    }
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

  const std::vector<std::vector<ROL::Ptr<FE<Real>>>> getBdryFE(void) const {
    return feBdry_;
  }

  const std::vector<std::vector<std::vector<int>>> getBdryCellLocIds(void) const {
    return bdryCellLocIds_;
  }

  const ROL::Ptr<FieldUtils::FieldInfo> getFieldInfo(void) const {
    return fieldInfo_;
  }

  const ROL::Ptr<FieldUtils::FieldInfo> getDensityFieldInfo(void) const {
    return fieldInfoDens_;
  }

  const ROL::Ptr<MultiMat_Load<Real>> getLoad(void) const {
    return load_;
  }

  const ROL::Ptr<MultiMat_Traction<Real>> getTraction(void) const {
    return traction_;
  }

  const std::vector<ROL::Ptr<MultiMat_MaterialTensor<Real>>> getMaterialTensor(void) const {
    return matTensor_;
  }

}; // PDE_TopoOpt

#endif
