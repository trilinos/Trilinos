// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde.hpp
    \brief Implements the local PDE interface for the Stefan_Boltzmann control problem.
*/

#ifndef BINARY_PDE_STEFAN_BOLTZMANN_HPP
#define BINARY_PDE_STEFAN_BOLTZMANN_HPP

#include "../../TOOLS/pde.hpp"
#include "../../TOOLS/fe.hpp"
#include "../../TOOLS/Intrepid_HGRAD_C0_FEM.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"
#include "Teuchos_LAPACK.hpp"

template <class Real>
class BinaryStefanBoltzmannPDE : public PDE<Real> {
private:

  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>> basisPtr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrs_;
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>> basisPtr2_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrs2_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real>> cellCub_;
  ROL::Ptr<Intrepid::Cubature<Real>> bdryCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real>> volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real>> fe_vol_, fe_ctrl_;
  std::vector<std::vector<ROL::Ptr<FE<Real>>>> fe_bdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fidx_;

  Real kappa_, sigma_, u0_, zpow_;
  int probDim_;
  bool useParam_;
  int nx_, ny_;
  Real dx_, dy_;
  std::vector<std::multimap<int,int>> ctrlPts_;
  std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> ctrlJac_;

public:

  BinaryStefanBoltzmannPDE(Teuchos::ParameterList &parlist) {
    kappa_ = parlist.sublist("Problem").get("Thermal Conductivity",    16.0);     // Watts/(meter Kelvin)
    sigma_ = parlist.sublist("Problem").get("Radiation Constant",      1.92e-10); // W/(meter Kelvin^3) in 2D and W/(meter^2 Kelvin^4) in 3D
    u0_    = parlist.sublist("Problem").get("Surrounding Temperature", 293.0);    // Kelvin
    zpow_  = parlist.sublist("Problem").get("Control Power",           2500.0);   // Watts
    // Finite element fields.
    int basisOrder = parlist.sublist("Problem").get("Order of FE Discretization", 1);
    int cubDegree  = parlist.sublist("Problem").get("Cubature Degree",            4);
    probDim_       = parlist.sublist("Problem").get("Problem Dimension",          2);
    if (probDim_ == 2) {
      basisPtr2_ = ROL::makePtr<Intrepid::Basis_HGRAD_C0_FEM<Real, Intrepid::FieldContainer<Real>>>();
      if (basisOrder == 1)
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      else if (basisOrder == 2)
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real>>>();
    }
    else if (probDim_ == 3) {
      basisPtr2_ = ROL::makePtr<Intrepid::Basis_HGRAD_C0_FEM<Real, Intrepid::FieldContainer<Real>>>();
      if (basisOrder == 1)
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_HEX_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      else if (basisOrder == 2)
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_HEX_C2_FEM<Real, Intrepid::FieldContainer<Real>>>();
    }
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    basisPtrs2_.clear(); basisPtrs2_.push_back(basisPtr2_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();        // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                       // create cubature factory
    cellCub_ = cubFactory.create(cellType, cubDegree);                       // create default cubature

    int d = cellType.getDimension();
    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2); // set cubature degree, e.g., 2
    bdryCub_ = cubFactory.create(bdryCellType, bdryCubDegree);

    // Nondimensionalize
    bool nondim = parlist.sublist("Problem").get("Nondimensionalize", true);
    if (nondim) {
      Real T0 = parlist.sublist("Problem").get("Reference Temperature",1000.0);
      u0_    /= T0;
      sigma_ *= std::pow(T0,probDim_)/kappa_;
      zpow_  /= T0*kappa_;
      kappa_  = static_cast<Real>(1);
    }

    // Parametric Control Information
    useParam_ = parlist.sublist("Problem").get("Use Parametric Control",   true);
    nx_       = parlist.sublist("Problem").get("Number X Control Patches", 4);
    ny_       = parlist.sublist("Problem").get("Number Y Control Patches", 4);
    dx_       = parlist.sublist("Problem").get("Width X Control Patches",  0.02);
    dy_       = parlist.sublist("Problem").get("Width Y Control Patches",  0.02);
  }
  
  void residual(ROL::Ptr<Intrepid::FieldContainer<Real>>             & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real>>              & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN()->dimension(0);
    const int f = fe_vol_->gradN()->dimension(1);
    const int p = fe_vol_->gradN()->dimension(2);
    const int d = fe_vol_->gradN()->dimension(3);
    // INITIALIZE RESIDUAL
    res = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // EVALUATE STATE ON FE BASIS
    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval, gradU_eval, valZ_eval;
    valU_eval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    gradU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    valZ_eval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    // INTEGRATE (kappa * grad(U)) . grad(N)
    Intrepid::RealSpaceTools<Real>::scale(*gradU_eval,kappa_);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *gradU_eval,
                                                  *fe_vol_->gradNdetJ(),
                                                  Intrepid::COMP_CPP, false);
    // INTEGRATE zpow * z
    if (useParam_) {
      applyParametricControl(valZ_eval,z_param);
    }
    else {
      fe_ctrl_->evaluateValue(valZ_eval, z_coeff);
      Intrepid::RealSpaceTools<Real>::scale(*valZ_eval,-zpow_);
    }
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *valZ_eval,
                                                  *fe_vol_->NdetJ(),
                                                  Intrepid::COMP_CPP, true);

    // APPLY STEFAN-BOLTZMANN CONDITIONS
    const int numLocalSideIds = bdryCellLocIds_[0].size();
    const int numCubPerSide = bdryCub_->getNumPoints();
    for (int j = 0; j < numLocalSideIds; ++j) {
      const int numCellsSide = bdryCellLocIds_[0][j].size();
      if (numCellsSide) {
        ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff_bdry, valU_eval_bdry, sb_valU, sbRes;
        valU_eval_bdry = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        sb_valU        = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        sbRes          = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f);
        // Get U coefficients on Stefan-Boltzmann boundary
        u_coeff_bdry = getBoundaryCoeff(*u_coeff, 0, j);
        // Evaluate U on FE basis
        fe_bdry_[0][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
        // Compute Stefan-Boltzmann residual
        computeStefanBoltzmann(sb_valU,valU_eval_bdry,j,0);
        Intrepid::FunctionSpaceTools::integrate<Real>(*sbRes,
                                                      *sb_valU,
                                                      *fe_bdry_[0][j]->NdetJ(),
                                                      Intrepid::COMP_CPP, false);
        // Add Stefan-Boltzmann residual to volume residual
        for (int k = 0; k < numCellsSide; ++k) {
          int cidx = bdryCellLocIds_[0][j][k];
          for (int l = 0; l < f; ++l) {
            (*res)(cidx,l) += (*sbRes)(k,l);
          }
        }
      }
    }
  }
  
  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real>>             & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN()->dimension(0);
    const int f = fe_vol_->gradN()->dimension(1);
    // INITILAIZE JACOBIAN
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // INTEGRATE (kappa * grad(N)) . grad(N)
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *fe_vol_->gradN(),
                                                  *fe_vol_->gradNdetJ(),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*jac,kappa_);

    // APPLY STEFAN-BOLTZMANN CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    const int numCubPerSide = bdryCub_->getNumPoints();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      if (numCellsSide) {
        ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff_bdry, valU_eval_bdry, sb_derivU, sb_derivU_N, sbJac;
        valU_eval_bdry = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        sb_derivU      = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        sb_derivU_N    = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f, numCubPerSide);
        sbJac          = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f, f);
        // Get U coefficients on Stefan-Boltzmann boundary
        u_coeff_bdry = getBoundaryCoeff(*u_coeff, 0, j);
        // Evaluate U on FE basis
        fe_bdry_[0][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
        // Compute Stefan-Boltzmann residual
        computeStefanBoltzmann(sb_derivU,valU_eval_bdry,j,1);
        Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*sb_derivU_N,
                                                                    *sb_derivU,
                                                                    *fe_bdry_[0][j]->N());
        Intrepid::FunctionSpaceTools::integrate<Real>(*sbJac,
                                                      *sb_derivU_N,
                                                      *fe_bdry_[0][j]->NdetJ(),
                                                      Intrepid::COMP_CPP, false);
        // Add Stefan-Boltzmann residual to volume residual
        for (int k = 0; k < numCellsSide; ++k) {
          int cidx = bdryCellLocIds_[0][j][k];
          for (int l = 0; l < f; ++l) {
            for (int m = 0; m < f; ++m) {
              (*jac)(cidx,l,m) += (*sbJac)(k,l,m);
            }
          }
        }
      }
    }
  }
  
  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real>>             & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              & z_param = ROL::nullPtr) {
    if (!useParam_) {
      // GET DIMENSIONS
      const int c  = fe_vol_->N()->dimension(0);
      const int f1 = fe_vol_->N()->dimension(1);
      const int f2 = fe_ctrl_->N()->dimension(1);
      // INITILAIZE JACOBIAN
      jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f1, f2);
      // INTEGRATE control Jacobian
      Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                    *fe_vol_->N(),
                                                    *fe_ctrl_->NdetJ(),
                                                    Intrepid::COMP_CPP, false);
      Intrepid::RealSpaceTools<Real>::scale(*jac,-zpow_);
    }
    else {
      throw Exception::Zero(">>> (BinaryStefanBoltzmannPDE::Jacobian_2): Jacobian is zero.");
    }
  }

  void Jacobian_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>>  & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>>  & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>               & z_param = ROL::nullPtr) {
    if (useParam_) {
      // GET DIMENSIONS
      const int c = fe_vol_->gradN()->dimension(0);
      const int f = fe_vol_->gradN()->dimension(1);
      // ADD CONTROL TERM TO RESIDUAL
      for (int i = 0; i < ny_; ++i) {
        for (int j = 0; j < nx_; ++j) {
          int ind = j + i*nx_;
          jac[ind] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
          Intrepid::FunctionSpaceTools::integrate<Real>(*jac[ind],
                                                        *ctrlJac_[ind],
                                                        *fe_vol_->NdetJ(),
                                                        Intrepid::COMP_CPP, false);
        }
      }
    }
    else {
      throw Exception::Zero(">>> (BinaryStefanBoltzmannPDE::Jacobian_3): Jacobian is zero.");
    }
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real>>             & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN()->dimension(0);
    const int f = fe_vol_->gradN()->dimension(1);
    // INITILAIZE JACOBIAN
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // APPLY STEFAN-BOLTZMANN CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    const int numCubPerSide = bdryCub_->getNumPoints();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      if (numCellsSide) {
        ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff_bdry, l_coeff_bdry, valU_eval_bdry, valL_eval_bdry;
        ROL::Ptr<Intrepid::FieldContainer<Real>> sb_derivU, sb_derivU_L, sb_derivU_L_N, sbHess;
        valU_eval_bdry = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        valL_eval_bdry = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        sb_derivU      = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        sb_derivU_L    = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        sb_derivU_L_N  = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f, numCubPerSide);
        sbHess         = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f, f);
        // Get U coefficients on Stefan-Boltzmann boundary
        u_coeff_bdry = getBoundaryCoeff(*u_coeff, 0, j);
        l_coeff_bdry = getBoundaryCoeff(*l_coeff, 0, j);
        // Evaluate U on FE basis
        fe_bdry_[0][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
        fe_bdry_[0][j]->evaluateValue(valL_eval_bdry, l_coeff_bdry);
        // Compute Stefan-Boltzmann residual
        computeStefanBoltzmann(sb_derivU,valU_eval_bdry,j,2);
        Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*sb_derivU_L,
                                                                   *sb_derivU,
                                                                   *valL_eval_bdry);
        Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*sb_derivU_L_N,
                                                                    *sb_derivU_L,
                                                                    *fe_bdry_[0][j]->N());
        Intrepid::FunctionSpaceTools::integrate<Real>(*sbHess,
                                                      *sb_derivU_L_N,
                                                      *fe_bdry_[0][j]->NdetJ(),
                                                      Intrepid::COMP_CPP, false);
        // Add Stefan-Boltzmann residual to volume residual
        for (int k = 0; k < numCellsSide; ++k) {
          int cidx = bdryCellLocIds_[0][j][k];
          for (int l = 0; l < f; ++l) { 
            for (int m = 0; m < f; ++m) { 
              (*hess)(cidx,l,m) += (*sbHess)(k,l,m);
            }
          }
        }
      }
    }
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real>>             & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (BinaryStefanBoltzmannPDE::Hessian_12): Hessian is zero.");
  }

  void Hessian_13(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>>  & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>>  & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>>  & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>               & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (BinaryStefanBoltzmannPDE::Hessian_13): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real>>             & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (BinaryStefanBoltzmannPDE::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real>>             & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (BinaryStefanBoltzmannPDE::Hessian_22): Hessian is zero.");
  }

  void Hessian_23(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>>  & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>>  & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>>  & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>               & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (BinaryStefanBoltzmannPDE::Hessian_23): Hessian is zero.");
  }

  void Hessian_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>>  & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>>  & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>>  & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>               & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (BinaryStefanBoltzmannPDE::Hessian_31): Hessian is zero.");
  }

  void Hessian_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>>  & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>>  & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>>  & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>               & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (BinaryStefanBoltzmannPDE::Hessian_32): Hessian is zero.");
  }

  void Hessian_33(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>>               & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>>               & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>>               & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>                            & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (BinaryStefanBoltzmannPDE::Hessian_33): Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    // GET DIMENSIONS
    int c = fe_vol_->N()->dimension(0);
    int f = fe_vol_->N()->dimension(1);
    // INITIALIZE RIESZ
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    *riesz = *fe_vol_->stiffMat();
    Intrepid::RealSpaceTools<Real>::add(*riesz,*(fe_vol_->massMat()));
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    // GET DIMENSIONS
    int c = fe_vol_->N()->dimension(0);
    int f = fe_ctrl_->N()->dimension(1);
    // INITIALIZE RIESZ
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    *riesz = *fe_ctrl_->massMat();
  }
 
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields(void) {
    return basisPtrs_;
  }
 
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields2(void) {
    if (useParam_) return basisPtrs_;
    return basisPtrs2_;
  }

  void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real>> &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> &bdryCellNodes, 
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds ) {
    volCellNodes_   = volCellNodes;
    bdryCellNodes_  = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition
    fe_vol_  = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtr_,cellCub_);
    //fe_ctrl_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtr2_,cellCub_,false);
    // Set local boundary DOFs
    fidx_ = fe_vol_->getBoundaryDofs();
    // Construct boundary FEs
    const int numSidesets = bdryCellNodes.size();
    fe_bdry_.resize(numSidesets);
    for(int i = 0; i < numSidesets; ++i) {
      int numLocSides = bdryCellNodes[i].size();
      fe_bdry_[i].resize(numLocSides);
      for (int j = 0; j < numLocSides; ++j) {
        if (bdryCellNodes_[i][j] != ROL::nullPtr) {
          fe_bdry_[i][j] = ROL::makePtr<FE<Real>>(bdryCellNodes_[i][j],basisPtr_,bdryCub_,j);
        }
      }
    }
    // Compute control cells
    if (useParam_) {
      const int c = fe_vol_->N()->dimension(0);
      const int p = fe_vol_->N()->dimension(2);
      std::vector<Real> x(probDim_);
      ctrlPts_.clear(); ctrlPts_.resize(nx_*ny_);
      ctrlJac_.clear(); ctrlJac_.resize(nx_*ny_);
      for (int i = 0; i < nx_*ny_; ++i) {
        ctrlPts_[i].clear();
        ctrlJac_[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
        ctrlJac_[i]->initialize();
      }
      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          for (int k = 0; k < probDim_; ++k) {
            x[k] = (*fe_vol_->cubPts())(i,j,k);
          }
          bool found = false;
          for (int l = 0; l < ny_ && !found; ++l) {
            Real cy = static_cast<Real>(0.5)*static_cast<Real>(2*l+1)/static_cast<Real>(ny_);
            if (std::abs(x[1]-cy)<=dy_) {
              for (int m = 0; m < nx_ && !found; ++m) {
                Real cx = static_cast<Real>(0.5)*static_cast<Real>(2*m+1)/static_cast<Real>(nx_);
                if (std::abs(x[0]-cx)<=dx_) {
                  ctrlPts_[m+l*nx_].insert({i,j});
                  (*ctrlJac_[m+l*nx_])(i,j) = -zpow_;
                  found = true;
                }
              }
            }
          }
        }
      }
    }
  }

  const ROL::Ptr<FE<Real>> getVolFE(void) const {
    return fe_vol_;
  }

  const ROL::Ptr<FE<Real>> getControlFE(void) const {
    return fe_ctrl_;
  }

  const std::vector<ROL::Ptr<FE<Real>>> getBdryFE(void) const {
    return fe_bdry_[0];
  }

  const std::vector<std::vector<int>> getBdryCellLocIds(void) const {
    return bdryCellLocIds_[0];
  }
 
private:

  Real evaluateStefanBoltzmann(Real u, int locSideId, int deriv = 0) const {
    if ( deriv == 1 ) {
      return sigma_*static_cast<Real>(probDim_+1)*std::pow(u,probDim_);
    }
    if ( deriv == 2 ) {
      return sigma_*static_cast<Real>(probDim_+1)*static_cast<Real>(probDim_)*std::pow(u,probDim_-1);
    }
    return sigma_*(std::pow(std::abs(u),probDim_)*u-std::pow(u0_,probDim_+1));
  }

  /***************************************************************************/
  /************** COMPUTE PDE COEFFICIENTS AT DOFS ***************************/
  /***************************************************************************/
  void computeStefanBoltzmann(ROL::Ptr<Intrepid::FieldContainer<Real>> &sb,
                              const ROL::Ptr<Intrepid::FieldContainer<Real>> &u,
                              int locSideId,
                              int deriv = 0) const {
    const int c = u->dimension(0);
    const int p = u->dimension(1);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*sb)(i,j) = evaluateStefanBoltzmann((*u)(i,j),locSideId,deriv);
      }
    }
  }

  void applyParametricControl(ROL::Ptr<Intrepid::FieldContainer<Real>> &Bz,
                              const ROL::Ptr<const std::vector<Real>>  &z) const {
    Bz->initialize();
    for (int i = 0; i < ny_; ++i) {
      for (int j = 0; j < nx_; ++j) {
        int ind = j + i*nx_;
        for (auto it = ctrlPts_[ind].begin(); it != ctrlPts_[ind].end(); ++it) {
          (*Bz)(it->first,it->second) = -zpow_*(*z)[ind];
        }
      }
    }
  }

  /***************************************************************************/
  /************** EXTRACT COEFFICIENTS ON BOUNDARY ***************************/
  /***************************************************************************/
  ROL::Ptr<Intrepid::FieldContainer<Real>> getBoundaryCoeff(
      const Intrepid::FieldContainer<Real> & cell_coeff,
      int sideSet, int cell) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][cell];
    const int numCellsSide = bdryCellLocId.size();
    const int f = basisPtr_->getCardinality();
    
    ROL::Ptr<Intrepid::FieldContainer<Real > > bdry_coeff = 
      ROL::makePtr<Intrepid::FieldContainer<Real >>(numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        (*bdry_coeff)(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }
}; // PDE_stefan_boltzmann

#endif
