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

#ifndef PDEOPT_PDE_STOCH_STEFAN_BOLTZMANNK_HPP
#define PDEOPT_PDE_STOCH_STEFAN_BOLTZMANNK_HPP

#include "../TOOLS/pdeK.hpp"
#include "../TOOLS/feK.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"

#include "ROL_Ptr.hpp"
#include "Teuchos_LAPACK.hpp"

template <class Real, class DeviceType>
class StochasticStefanBoltzmannPDE : public PDE<Real,DeviceType> {
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
  ROL::Ptr<fe_type> fe_vol_;
  std::vector<std::vector<ROL::Ptr<fe_type>>> fe_bdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fidx_;

  const Real scale_;
  Real xmid_;
  Real engTemp_;
  Real airTemp_;
  Real H2OTemp_;
  Real advMag_;
  Real SBscale_;
  Real nonLin_;

public:

  StochasticStefanBoltzmannPDE(ROL::ParameterList &parlist) : scale_(1) {
    xmid_ = parlist.sublist("Geometry").get<Real>("Step height");
    engTemp_ = parlist.sublist("Problem").get("Engine: Ambient Temperature",450.0);
    airTemp_ = parlist.sublist("Problem").get("Air: Ambient Temperature",293.0);
    H2OTemp_ = parlist.sublist("Problem").get("Water: Ambient Temperature",303.0);
    advMag_  = parlist.sublist("Problem").get("Advection Magnitude",6.3);
    bool useSB = parlist.sublist("Problem").get("Use Stefan-Boltzmann",true);
    SBscale_ = (useSB) ? static_cast<Real>(1) : static_cast<Real>(0);
    bool useND = parlist.sublist("Problem").get("Use Nonlinear Conductivities",true);
    nonLin_ = (useND) ? static_cast<Real>(1) : static_cast<Real>(0);
    // Finite element fields.
    int basisOrder = parlist.sublist("Problem").get("Basis Order",1);
    if (basisOrder == 1)
      basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType,Real,Real>>();
    else if (basisOrder == 2)
      basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<DeviceType,Real,Real>>();
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();        // get the cell type from any basis
    Intrepid2::DefaultCubatureFactory cubFactory;                       // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Volume Cubature Degree",2); // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create<DeviceType,Real,Real>(cellType, cubDegree);                       // create default cubature

    int d = cellType.getDimension();
    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2); // set cubature degree, e.g., 2
    bdryCub_ = cubFactory.create<DeviceType,Real,Real>(bdryCellType, bdryCubDegree);
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
    // EVALUATE STATE ON FE BASIS
    scalar_view U_eval("U_eval", c, p);
    fe_vol_->evaluateValue(U_eval, u_coeff);
    scalar_view gradU_eval("gradU_eval", c, p, d);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    // COMPUTE CONSTANT PDE COEFFICIENTS
    scalar_view V("V", c, p, d);
    scalar_view rhs("rhs", c, p);
    computeCoefficients(V,rhs,z_param);
    // COMPUTE DIFFUSIVITY
    scalar_view kappa("kappa", c, p);
    computeDiffusivity(kappa,U_eval,0);
    // MULTIPLY kappa * grad(U)
    scalar_view kappa_gradU("kappa_gradU", c, p, d);
    fst::tensorMultiplyDataData(kappa_gradU,kappa,gradU_eval);
    // INTEGRATE (kappa * grad(U)) . grad(N)
    fst::integrate(res,kappa_gradU,fe_vol_->gradNdetJ(),false);
    // MULTIPLY V . grad(U)
    scalar_view V_gradU("V_gradU", c, p);
    fst::dotMultiplyDataData(V_gradU,V,gradU_eval);
    // INTEGRATE (V . grad(U)) * N
    fst::integrate(res,V_gradU,fe_vol_->NdetJ(),true);
    // ADD RHS TO RESIDUAL
    fst::integrate(res,rhs,fe_vol_->NdetJ(),true);

    // APPLY NEUMANN CONDITIONS: Sideset 3, 6
    // ---> Nothing to do
    int numLocalSideIds(0);
    // APPLY STEFAN-BOLTZMANN CONDITIONS: Sideset 1, 2, 4, 5
    std::vector<int> sidesets = {1, 2, 4, 5};
    for (int i = 0; i < 4; ++i) {
      numLocalSideIds = bdryCellLocIds_[sidesets[i]].size();
      const int numCubPerSide = bdryCub_->getNumPoints();
      for (int j = 0; j < numLocalSideIds; ++j) {
        int numCellsSide = bdryCellLocIds_[sidesets[i]][j].size();
        if (numCellsSide) {
          // Get U coefficients on Stefan-Boltzmann boundary
          scalar_view u_coeff_bdry = getBoundaryCoeff(u_coeff, sidesets[i], j);
          // Evaluate U on FE basis
          scalar_view valU_eval_bdry("valU_eval_bdry", numCellsSide, numCubPerSide);
          fe_bdry_[sidesets[i]][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
          // Compute Stefan-Boltzmann residual
          scalar_view sb_valU("sb_valU", numCellsSide, numCubPerSide);
          computeStefanBoltzmann(sb_valU,valU_eval_bdry,sidesets[i],j,0);
          scalar_view sbRes("sbRes", numCellsSide, f);
          fst::integrate(sbRes,sb_valU,fe_bdry_[sidesets[i]][j]->NdetJ(),false);
          // Add Stefan-Boltzmann residual to volume residual
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[sidesets[i]][j][k];
            for (int l = 0; l < f; ++l) { 
              res(cidx,l) += sbRes(k,l);
            }
          }
        }
      }
    }
    // APPLY ROBIN CONTROLS: Sideset 0
    int sideset = 0;
    numLocalSideIds = bdryCellLocIds_[sideset].size();
    const int numCubPerSide = bdryCub_->getNumPoints();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[sideset][j].size();
      if (numCellsSide) {
        // Get U coefficients on Robin boundary
        scalar_view u_coeff_bdry = getBoundaryCoeff(u_coeff, sideset, j);
        // Evaluate U on FE basis
        scalar_view valU_eval_bdry("valU_eval_bdry", numCellsSide, numCubPerSide);
        fe_bdry_[sideset][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
        // Compute Stefan-Boltzmann residual
        scalar_view robinVal("robinVal", numCellsSide, numCubPerSide);
        scalar_view valZ_eval_bdry;
        if (z_coeff != scalar_view()) {
          // Get Z coefficients on Robin boundary
          scalar_view z_coeff_bdry = getBoundaryCoeff(z_coeff, sideset, j);
          // Evaluate Z on FE basis
          valZ_eval_bdry = scalar_view("valZ_eval_bdry", numCellsSide, numCubPerSide);
          fe_bdry_[sideset][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
        }
        else {
          valZ_eval_bdry = scalar_view();
        }
        computeRobin(robinVal,valU_eval_bdry,valZ_eval_bdry,sideset,j,0);
        scalar_view robinRes("robinRes", numCellsSide, f);
        fst::integrate(robinRes,robinVal,fe_bdry_[sideset][j]->NdetJ(),false);
        // Add Stefan-Boltzmann residual to volume residual
        for (int k = 0; k < numCellsSide; ++k) {
          int cidx = bdryCellLocIds_[sideset][j][k];
          for (int l = 0; l < f; ++l) { 
            res(cidx,l) += robinRes(k,l);
          }
        }
      }
    }
//    // APPLY DIRICHLET CONTROLS: Sideset 0
//    int sideset = 0;
//    numLocalSideIds = bdryCellLocIds_[sideset].size();
//    for (int j = 0; j < numLocalSideIds; ++j) {
//      int numCellsSide = bdryCellLocIds_[sideset][j].size();
//      int numBdryDofs = fidx_[j].size();
//      for (int k = 0; k < numCellsSide; ++k) {
//        int cidx = bdryCellLocIds_[sideset][j][k];
//        for (int l = 0; l < numBdryDofs; ++l) {
//          (*res)(cidx,fidx_[j][l])
//            = (*u_coeff)(cidx,fidx_[j][l]) - (*z_coeff)(cidx,fidx_[j][l]);
//        }
//      }
//    }
  }
  
  void Jacobian_1(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // GET DIMENSIONS
    int c = fe_vol_->gradN().extent_int(0);
    int f = fe_vol_->gradN().extent_int(1);
    int p = fe_vol_->gradN().extent_int(2);
    int d = fe_vol_->gradN().extent_int(3);
    // INITILAIZE JACOBIAN
    jac = scalar_view("jac", c, f, f);
    // EVALUATE STATE ON FE BASIS
    scalar_view U_eval("U_eval", c, p);
    fe_vol_->evaluateValue(U_eval, u_coeff);
    scalar_view gradU_eval("gradU_eval", c, p, d);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    // COMPUTE CONSTNAT PDE COEFFICIENTS
    scalar_view V("V", c, p, d);
    scalar_view rhs("rhs", c, p);
    computeCoefficients(V,rhs,z_param);
    // COMPUTE DIFFUSIVITY
    scalar_view kappa("kappa", c, p);
    computeDiffusivity(kappa,U_eval,0);
    scalar_view d_kappa("d_kappa", c, p);
    computeDiffusivity(d_kappa,U_eval,1);
    // MULTIPLY kappa * grad(N)
    scalar_view kappa_gradN("kappa_gradN", c, f, p, d);
    fst::tensorMultiplyDataField(kappa_gradN,kappa,fe_vol_->gradN());
    // INTEGRATE (kappa * grad(N)) . grad(N)
    fst::integrate(jac,kappa_gradN,fe_vol_->gradNdetJ(),false);
    // MULTIPLY d_kappa * grad(U)
    scalar_view d_kappa_gradU("kappa_gradU", c, p, d);
    fst::tensorMultiplyDataData(d_kappa_gradU,d_kappa,gradU_eval);
    // MULTIPLY (d_kappa * grad(U)) . grad(N)
    scalar_view d_kappa_gradU_gradN("d_kappa_gradU_gradN", c, f, p);
    fst::dotMultiplyDataField(d_kappa_gradU_gradN,d_kappa_gradU,fe_vol_->gradNdetJ());
    // INTEGRATE (d_kappa * grad(U)) . grad(N) * N
    fst::integrate(jac,d_kappa_gradU_gradN,fe_vol_->N(),true);
    // MULTIPLY V . grad(N)
    scalar_view V_gradN("V_gradN", c, f, p);
    fst::dotMultiplyDataField(V_gradN,V,fe_vol_->gradN());
    // INTEGRATE (V . grad(U)) * N
    fst::integrate(jac,fe_vol_->NdetJ(),V_gradN,true);

    // APPLY NEUMANN CONDITIONS: Sideset 3, 6
    // ---> Nothing to do
    int numLocalSideIds(0);
    // APPLY STEFAN-BOLTZMANN CONDITIONS: Sideset 1, 2, 4, 5
    std::vector<int> sidesets = {1, 2, 4, 5};
    for (int i = 0; i < 4; ++i) {
      numLocalSideIds = bdryCellLocIds_[sidesets[i]].size();
      const int numCubPerSide = bdryCub_->getNumPoints();
      for (int j = 0; j < numLocalSideIds; ++j) {
        int numCellsSide = bdryCellLocIds_[sidesets[i]][j].size();
        if (numCellsSide) {
          // Get U coefficients on Stefan-Boltzmann boundary
          scalar_view u_coeff_bdry = getBoundaryCoeff(u_coeff, sidesets[i], j);
          // Evaluate U on FE basis
          scalar_view valU_eval_bdry("valU_eval_bdry", numCellsSide, numCubPerSide);
          fe_bdry_[sidesets[i]][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
          // Compute Stefan-Boltzmann residual
          scalar_view sb_derivU("sb_derivU", numCellsSide, numCubPerSide);
          computeStefanBoltzmann(sb_derivU,valU_eval_bdry,sidesets[i],j,1);
          scalar_view sb_derivU_N("sb_derivU_N", numCellsSide, f, numCubPerSide);
          fst::scalarMultiplyDataField(sb_derivU_N,sb_derivU,fe_bdry_[sidesets[i]][j]->N());
          scalar_view sbJac("sbJac", numCellsSide, f, f);
          fst::integrate(sbJac,sb_derivU_N,fe_bdry_[sidesets[i]][j]->NdetJ(),false);
          // Add Stefan-Boltzmann residual to volume residual
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[sidesets[i]][j][k];
            for (int l = 0; l < f; ++l) { 
              for (int m = 0; m < f; ++m) { 
                jac(cidx,l,m) += sbJac(k,l,m);
              }
            }
          }
        }
      }
    }
    // APPLY ROBIN CONTROL: Sideset 0
    int sideset = 0;
    numLocalSideIds = bdryCellLocIds_[sideset].size();
    const int numCubPerSide = bdryCub_->getNumPoints();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[sideset][j].size();
      if (numCellsSide) {
        // Get U coefficients on Robin boundary
        scalar_view u_coeff_bdry = getBoundaryCoeff(u_coeff, sideset, j);
        // Evaluate U on FE basis
        scalar_view valU_eval_bdry("valU_eval_bdry", numCellsSide, numCubPerSide);
        fe_bdry_[sideset][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
        // Compute Stefan-Boltzmann residual
        scalar_view robinVal("robinVal", numCellsSide, numCubPerSide);
        scalar_view valZ_eval_bdry;
        if (z_coeff != scalar_view()) {
          // Get Z coefficients on Robin boundary
          scalar_view z_coeff_bdry = getBoundaryCoeff(z_coeff, sideset, j);
          // Evaluate Z on FE basis
          valZ_eval_bdry = scalar_view("valZ_eval_bdry", numCellsSide, numCubPerSide);
          fe_bdry_[sideset][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
        }
        else {
          valZ_eval_bdry = scalar_view();
        }
        computeRobin(robinVal,valU_eval_bdry,valZ_eval_bdry,sideset,j,1,1);
        scalar_view robinVal_N("robinVal_N", numCellsSide, f, numCubPerSide);
        fst::scalarMultiplyDataField(robinVal_N,robinVal,fe_bdry_[sideset][j]->N());
        scalar_view robinJac("robinJac", numCellsSide, f, f);
        fst::integrate(robinJac,robinVal_N,fe_bdry_[sideset][j]->NdetJ(),false);
        // Add Stefan-Boltzmann residual to volume residual
        for (int k = 0; k < numCellsSide; ++k) {
          int cidx = bdryCellLocIds_[sideset][j][k];
          for (int l = 0; l < f; ++l) { 
            for (int m = 0; m < f; ++m) { 
              jac(cidx,l,m) += robinJac(k,l,m);
            }
          }
        }
      }
    }
//    // APPLY DIRICHLET CONDITIONS: Sideset 0
//    int sideset = 0;
//    numLocalSideIds = bdryCellLocIds_[sideset].size();
//    for (int j = 0; j < numLocalSideIds; ++j) {
//      int numCellsSide = bdryCellLocIds_[sideset][j].size();
//      int numBdryDofs = fidx_[j].size();
//      for (int k = 0; k < numCellsSide; ++k) {
//        int cidx = bdryCellLocIds_[sideset][j][k];
//        for (int l = 0; l < numBdryDofs; ++l) {
//          for (int m = 0; m < f; ++m) {
//            (*jac)(cidx,fidx_[j][l],m) = static_cast<Real>(0);
//          }
//          (*jac)(cidx,fidx_[j][l],fidx_[j][l]) = static_cast<Real>(1);
//        }
//      }
//    }
  }
  
  void Jacobian_2(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if (z_coeff != scalar_view()) {
      // GET DIMENSIONS
      int c = fe_vol_->gradN().extent_int(0);
      int f = fe_vol_->gradN().extent_int(1);
      // INITILAIZE JACOBIAN
      jac = scalar_view("jac", c, f, f);
      // APPLY ROBIN CONTROL: Sideset 0
      int sideset = 0;
      int numLocalSideIds = bdryCellLocIds_[sideset].size();
      const int numCubPerSide = bdryCub_->getNumPoints();
      for (int j = 0; j < numLocalSideIds; ++j) {
        int numCellsSide = bdryCellLocIds_[sideset][j].size();
        if (numCellsSide) {
          // Get U coefficients on Robin boundary
          scalar_view u_coeff_bdry = getBoundaryCoeff(u_coeff, sideset, j);
          // Get Z coefficients on Robin boundary
          scalar_view z_coeff_bdry = getBoundaryCoeff(z_coeff, sideset, j);
          // Evaluate U on FE basis
          scalar_view valU_eval_bdry("valU_eval_bdry", numCellsSide, numCubPerSide);
          fe_bdry_[sideset][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
          // Evaluate Z on FE basis
          scalar_view valZ_eval_bdry("valZ_eval_bdry", numCellsSide, numCubPerSide);
          fe_bdry_[sideset][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
          // Compute Stefan-Boltzmann residual
          scalar_view robinVal("robinVal", numCellsSide, numCubPerSide);
          computeRobin(robinVal,valU_eval_bdry,valZ_eval_bdry,sideset,j,1,2);
          scalar_view robinVal_N("robinVal_N", numCellsSide, f, numCubPerSide);
          fst::scalarMultiplyDataField(robinVal_N,robinVal,fe_bdry_[sideset][j]->N());
          scalar_view robinJac("robinJac", numCellsSide, f, f);
          fst::integrate(robinJac,robinVal_N,fe_bdry_[sideset][j]->NdetJ(),false);
          // Add Stefan-Boltzmann residual to volume residual
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[sideset][j][k];
            for (int l = 0; l < f; ++l) { 
              for (int m = 0; m < f; ++m) { 
                jac(cidx,l,m) += robinJac(k,l,m);
              }
            }
          }
        }
      }
    }
    else {
      throw Exception::Zero(">>> (StochasticStefanBoltzmannPDE::Jacobian_2): Jacobian_2 is zero.");
    }
//    // APPLY DIRICHLET CONTROLS: Sideset 0
//    int sideset = 0;
//    int numLocalSideIds = bdryCellLocIds_[sideset].size();
//    for (int j = 0; j < numLocalSideIds; ++j) {
//      int numCellsSide = bdryCellLocIds_[sideset][j].size();
//      int numBdryDofs = fidx_[j].size();
//      for (int k = 0; k < numCellsSide; ++k) {
//        int cidx = bdryCellLocIds_[sideset][j][k];
//        for (int l = 0; l < numBdryDofs; ++l) {
//          (*jac)(cidx,fidx_[j][l],fidx_[j][l]) = static_cast<Real>(-1);
//        }
//      }
//    }
  }

  void Jacobian_3(std::vector<scalar_view> & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if (z_param != ROL::nullPtr) {
      // GET DIMENSIONS
      int c = fe_vol_->gradN().extent_int(0);
      int f = fe_vol_->gradN().extent_int(1);
      int p = fe_vol_->gradN().extent_int(2);
      int d = fe_vol_->gradN().extent_int(3);
      // INITIALIZE RESIDUAL
      jac.resize(z_param->size());
      jac[0] = scalar_view("jac", c, f);
      // EVALUATE STATE ON FE BASIS
      scalar_view U_eval("U_eval", c, p);
      fe_vol_->evaluateValue(U_eval, u_coeff);
      scalar_view gradU_eval("gradU_eval", c, p, d);
      fe_vol_->evaluateGradient(gradU_eval, u_coeff);
      ROL::Ptr<std::vector<Real>> one = ROL::makePtr<std::vector<Real>>(z_param->size(), 1); 
      // COMPUTE CONSTANT PDE COEFFICIENTS
      scalar_view V("V", c, p, d);
      scalar_view rhs("rhs", c, p);
      computeCoefficients(V,rhs,one);
      // MULTIPLY V . grad(U)
      scalar_view V_gradU("V_gradU", c, p);
      fst::dotMultiplyDataData(V_gradU,V,gradU_eval);
      // INTEGRATE (V . grad(U)) * N
      fst::integrate(jac[0],V_gradU,fe_vol_->NdetJ(),true);
    }
    else{
      throw Exception::Zero(">>> (StochasticStefanBoltzmannPDE::Jacobian_3): Jacobian_3 is zero.");
    }
  }

  void Hessian_11(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // GET DIMENSIONS
    int c = fe_vol_->gradN().extent_int(0);
    int f = fe_vol_->gradN().extent_int(1);
    int p = fe_vol_->gradN().extent_int(2);
    int d = fe_vol_->gradN().extent_int(3);
    // INITILAIZE JACOBIAN
    hess = scalar_view("hess", c, f, f);
    scalar_view l_coeff_dbc("l_coeff_dbc", c, f);
    Kokkos::deep_copy(l_coeff_dbc, l_coeff);
//    // APPLY DIRICHLET CONDITIONS TO LAGRANGE MULTIPLIERS: Sideset 0
//    int sideset = 0;
//    int numLocalSideIds = bdryCellLocIds_[sideset].size();
//    for (int j = 0; j < numLocalSideIds; ++j) {
//      int numCellsSide = bdryCellLocIds_[sideset][j].size();
//      int numBdryDofs = fidx_[j].size();
//      for (int k = 0; k < numCellsSide; ++k) {
//        int cidx = bdryCellLocIds_[sideset][j][k];
//        for (int l = 0; l < numBdryDofs; ++l) {
//          (*l_coeff_dbc)(cidx,fidx_[j][l]) = static_cast<Real>(0);
//        }
//      }
//    }
    // EVALUATE STATE ON FE BASIS
    scalar_view U_eval("U_eval", c, p);
    fe_vol_->evaluateValue(U_eval, u_coeff);
    scalar_view gradU_eval("gradU_eval", c, p, d);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    scalar_view gradL_eval("gradL_eval", c, p, d);
    fe_vol_->evaluateGradient(gradL_eval, l_coeff_dbc);
    // COMPUTE DIFFUSIVITY
    scalar_view d1_kappa("d1_kappa", c, p);
    computeDiffusivity(d1_kappa,U_eval,1);
    scalar_view d2_kappa("d2_kappa", c, p);
    computeDiffusivity(d2_kappa,U_eval,2);
    // MULTIPLY d1_kappa * grad(L)
    scalar_view d1_kappa_gradL("d1_kappa_gradL", c, p, d);
    fst::tensorMultiplyDataData(d1_kappa_gradL,d1_kappa,gradL_eval);
    // MULTIPLY (d_1kappa * grad(L)) . grad(N)det(J)
    scalar_view d1_kappa_gradL_gradNdetJ("d1_kappa_gradL_gradNdetJ", c, f, p);
    fst::dotMultiplyDataField(d1_kappa_gradL_gradNdetJ,d1_kappa_gradL,fe_vol_->gradNdetJ());
    // INTEGRATE (d1_kappa * grad(L)) . grad(N) * N
    fst::integrate(hess,fe_vol_->N(),d1_kappa_gradL_gradNdetJ,false);
    // MULTIPLY d1_kappa * grad(L) . grad(N)
    scalar_view d1_kappa_gradL_gradN("d1_kappa_gradL_gradN", c, f, p);
    fst::dotMultiplyDataField(d1_kappa_gradL_gradN,d1_kappa_gradL,fe_vol_->gradN());
    // INTEGRATE (d1_kappa * grad(L)) . grad(N) * N
    fst::integrate(hess,d1_kappa_gradL_gradN,fe_vol_->NdetJ(),true);
    // MULTIPLY grad(U) . grad(L)
    scalar_view gradU_gradL("gradU_gradL", c, p);
    fst::dotMultiplyDataData(gradU_gradL,gradU_eval,gradL_eval);
    // MULTIPLY d2_kappa * grad(U) . grad(L)
    scalar_view d2_kappa_gradU_gradL("d2_kappa_gradU_gradL", c, p);
    fst::scalarMultiplyDataData(d2_kappa_gradU_gradL,d2_kappa,gradU_gradL);
    // MULTIPLY d2_kappa * grad(U) . grad(L) * N
    scalar_view d2_kappa_gradU_gradL_N("d2_kappa_gradU_gradL_N", c, f, p);
    fst::scalarMultiplyDataField(d2_kappa_gradU_gradL_N,d2_kappa_gradU_gradL,fe_vol_->N());
    // INTEGRATE (d2_kappa * grad(U) . grad(L) ) * N * N
    fst::integrate(hess,d2_kappa_gradU_gradL_N,fe_vol_->NdetJ(),true);
    // APPLY NEUMANN CONDITIONS: Sideset 3, 6
    // ---> Nothing to do
    // APPLY STEFAN-BOLTZMANN CONDITIONS: Sideset 1, 2, 4, 5
    std::vector<int> sidesets = {1, 2, 4, 5};
    for (int i = 0; i < 4; ++i) {
      int numLocalSideIds = bdryCellLocIds_[sidesets[i]].size();
      const int numCubPerSide = bdryCub_->getNumPoints();
      for (int j = 0; j < numLocalSideIds; ++j) {
        int numCellsSide = bdryCellLocIds_[sidesets[i]][j].size();
        if (numCellsSide) {
          // Get U coefficients on Stefan-Boltzmann boundary
          scalar_view u_coeff_bdry = getBoundaryCoeff(u_coeff, sidesets[i], j);
          scalar_view l_coeff_bdry = getBoundaryCoeff(l_coeff_dbc, sidesets[i], j);
          // Evaluate U on FE basis
          scalar_view valU_eval_bdry("valU_eval_bdry", numCellsSide, numCubPerSide);
          fe_bdry_[sidesets[i]][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
          scalar_view valL_eval_bdry("valL_eval_bdry", numCellsSide, numCubPerSide);
          fe_bdry_[sidesets[i]][j]->evaluateValue(valL_eval_bdry, l_coeff_bdry);
          // Compute Stefan-Boltzmann residual
          scalar_view sb_derivU("sb_derivU", numCellsSide, numCubPerSide);
          computeStefanBoltzmann(sb_derivU,valU_eval_bdry,sidesets[i],j,2);
          scalar_view sb_derivU_L("sb_derivU_L", numCellsSide, numCubPerSide);
          fst::scalarMultiplyDataData(sb_derivU_L,sb_derivU,valL_eval_bdry);
          scalar_view sb_derivU_L_N("sb_derivU_L_N", numCellsSide, f, numCubPerSide);
          fst::scalarMultiplyDataField(sb_derivU_L_N,sb_derivU_L,fe_bdry_[sidesets[i]][j]->N());
          scalar_view sbHess("sbHess", numCellsSide, f, f);
          fst::integrate(sbHess,sb_derivU_L_N,fe_bdry_[sidesets[i]][j]->NdetJ(),false);
          // Add Stefan-Boltzmann residual to volume residual
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[sidesets[i]][j][k];
            for (int l = 0; l < f; ++l) { 
              for (int m = 0; m < f; ++m) { 
                hess(cidx,l,m) += sbHess(k,l,m);
              }
            }
          }
        }
      }
    }
    // APPLY ROBIN CONTROL: Sideset 0
    // --> Nothing to do
  }

  void Hessian_12(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (StochasticStefanBoltzmannPDE::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (StochasticStefanBoltzmannPDE::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (StochasticStefanBoltzmannPDE::Hessian_22): Hessian is zero.");
  }

  void Hessian_13(std::vector<scalar_view> & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if ( z_param != ROL::nullPtr ) {
      // GET DIMENSIONS
      int c = fe_vol_->gradN().extent_int(0);
      int f = fe_vol_->gradN().extent_int(1);
      int p = fe_vol_->gradN().extent_int(2);
      int d = fe_vol_->gradN().extent_int(3);
      // INITILAIZE HESSIAN
      hess.resize(z_param->size());
      hess[0] = scalar_view("hess", c, f);
      // EVALUATE STATE ON FE BASIS
      scalar_view L_eval("L_eval", c, p);
      fe_vol_->evaluateValue(L_eval, l_coeff);
      // COMPUTE CONSTANT PDE COEFFICIENTS
      scalar_view V("V", c, p, d);
      scalar_view rhs("rhs", c, p);
      ROL::Ptr<std::vector<Real>> one = ROL::makePtr<std::vector<Real>>(z_param->size(), 1); 
      computeCoefficients(V,rhs,one);
      // MULTIPLY V . grad(N)
      scalar_view V_gradN("V_gradN", c, f, p);
      fst::dotMultiplyDataField(V_gradN,V,fe_vol_->gradNdetJ());
      // INTEGRATE (V . grad(U)) * N
      fst::integrate(hess[0],L_eval,V_gradN,false);
    }
    else {
      throw Exception::Zero(">>> (StochasticStefanBoltzmannPDE::Hessian_13): Hessian_13 is zero.");
    }
  }

  void Hessian_23(std::vector<scalar_view> & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (StochasticStefanBoltzmannPDE::Hessian_23): Hessian_23 is zero.");
  }

  void Hessian_31(std::vector<scalar_view> & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if ( z_param != ROL::nullPtr ) {
      // GET DIMENSIONS
      int c = fe_vol_->gradN().extent_int(0);
      int f = fe_vol_->gradN().extent_int(1);
      int p = fe_vol_->gradN().extent_int(2);
      int d = fe_vol_->gradN().extent_int(3);
      // INITILAIZE HESSIAN
      hess.resize(z_param->size());
      hess[0] = scalar_view("hess", c, f);
      // EVALUATE STATE ON FE BASIS
      scalar_view L_eval("L_eval", c, p);
      fe_vol_->evaluateValue(L_eval, l_coeff);
      // COMPUTE CONSTANT PDE COEFFICIENTS
      scalar_view V("V", c, p, d);
      scalar_view rhs("rhs", c, p);
      ROL::Ptr<std::vector<Real>> one = ROL::makePtr<std::vector<Real>>(z_param->size(), 1); 
      computeCoefficients(V,rhs,one);
      // MULTIPLY V . grad(N)
      scalar_view V_gradN("V_gradN", c, f, p);
      fst::dotMultiplyDataField(V_gradN,V,fe_vol_->gradNdetJ());
      // INTEGRATE (V . grad(U)) * N
      fst::integrate(hess[0],L_eval,V_gradN,false);
    }
    else {
      throw Exception::Zero(">>> (StochasticStefanBoltzmannPDE::Hessian_31): Hessian_31 is zero.");
    }
  }

  void Hessian_32(std::vector<scalar_view> & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (StochasticStefanBoltzmannPDE::Hessian_32): Hessian_32 is zero.");
  }

  void Hessian_33(std::vector<std::vector<scalar_view>> & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (StochasticStefanBoltzmannPDE::Hessian_33): Hessian_33 is zero.");
  }

  void RieszMap_1(scalar_view & riesz) override {
    // GET DIMENSIONS
    int c = fe_vol_->N().extent_int(0);
    int f = fe_vol_->N().extent_int(1);
    // INITIALIZE RIESZ
    riesz = scalar_view("riesz", c, f, f);
    Kokkos::deep_copy(riesz, fe_vol_->stiffMat());
    rst::add(riesz,fe_vol_->massMat());
  }

  void RieszMap_2(scalar_view & riesz) override {
    throw Exception::NotImplemented(">>> (StochasticStefanBoltzmannPDE::RieszMap2): Not implemented.");
    // GET DIMENSIONS
    int c = fe_vol_->N().extent_int(0);
    int f = fe_vol_->N().extent_int(1);
    // INITIALIZE RIESZ
    riesz = scalar_view("riesz", c, f, f);
    Kokkos::deep_copy(riesz,fe_vol_->massMat());
  }
 
  std::vector<basis_ptr> getFields(void) override {
    return basisPtrs_;
  }

  void setCellNodes(const scalar_view &volCellNodes,
                    const std::vector<std::vector<scalar_view>> &bdryCellNodes, 
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds ) override {
    volCellNodes_   = volCellNodes;
    bdryCellNodes_  = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition
    fe_vol_ = ROL::makePtr<fe_type>(volCellNodes_,basisPtr_,cellCub_);
    // Set local boundary DOFs
    fidx_ = fe_vol_->getBoundaryDofs();
    // Construct boundary FEs
    const int numSidesets = bdryCellNodes.size();
    fe_bdry_.resize(numSidesets);
    for(int i = 0; i < numSidesets; ++i) {
      int numLocSides = bdryCellNodes[i].size();
      fe_bdry_[i].resize(numLocSides);
      for (int j = 0; j < numLocSides; ++j) {
        if (bdryCellNodes_[i][j] != scalar_view()) {
          fe_bdry_[i][j] = ROL::makePtr<fe_type>(bdryCellNodes_[i][j],basisPtr_,bdryCub_,j);
        }
      }
    }
  }

  const ROL::Ptr<fe_type> getVolFE(void) const {
    return fe_vol_;
  }

  const std::vector<ROL::Ptr<fe_type>> getBdryFE(const int sideset) const {
    return fe_bdry_[sideset];
  }

  const std::vector<std::vector<int>> getBdryCellLocIds(const int sideset) const {
    return bdryCellLocIds_[sideset];
  }
 
private:
    
  /***************************************************************************/
  /************** EVALUATE PDE COEFFICIENTS AT DOF COORDINATES ***************/
  /***************************************************************************/
  Real evaluateDiffusivity(Real u, const std::vector<Real> & x, int deriv = 0) const {
    const std::vector<Real> param = PDE<Real,DeviceType>::getParameter();
    const int size_w = 11, size_a = 9;
    if ( x[1] < xmid_ ) {
      // Water Thermal Conductivity: 0.5818 at 280K and 0.6797 at 370K
      std::vector<Real> param_w(size_w,0);
      for (int i = 0; i < size_w; ++i) {
        Real pi = (static_cast<int>(param.size()) > i ? param[i] : static_cast<Real>(0.5));
        param_w[i] = 0.01*pi;
      }
      std::vector<Real> c(3,0);
      getWaterCoeff(c,param_w);
      Real val = c[0] + nonLin_*(c[1] * u + c[2] * u * u);
      const Real min = 0.1;
      if ( deriv == 1 ) {
        return (val < min ? static_cast<Real>(0) : nonLin_*(c[1] + static_cast<Real>(2)*c[2]*u));
      }
      if ( deriv == 2 ) {
        return (val < min ? static_cast<Real>(0) : nonLin_*(static_cast<Real>(2)*c[2]));
      }
      return (val < min ? min : c[0] + nonLin_*(c[1] * u + c[2] * u * u));
    }
    else {
      // Aluminum Thermal Conductivity: 236 at 273K and 240 at 400K
      std::vector<Real> param_a(size_a,0);
      for (int i = 0; i < size_a; ++i) {
        Real pi = (static_cast<int>(param.size()) > size_w+i ? param[size_w+i] : static_cast<Real>(0.5));
        param_a[i] = pi;
      }
      std::vector<Real> c(5,0);
      getAluminumCoeff(c,param_a);
      Real u2 = u*u, u3 = u2*u, u4 = u3*u;
      Real val = c[0] + nonLin_*(c[1]*u + c[2]*u2 + c[3]*u3 + c[4]*u4);
      const Real min = 100.0;
      if ( deriv == 1 ) {
        return (val < min ? static_cast<Real>(0) : 
               nonLin_*(c[1] + static_cast<Real>(2)*c[2]*u
                             + static_cast<Real>(3)*c[3]*u2
                             + static_cast<Real>(4)*c[4]*u3));
      }
      if ( deriv == 2 ) {
        return (val < min ? static_cast<Real>(0) : 
               nonLin_*(static_cast<Real>(2)*c[2] +  static_cast<Real>(6)*c[3]*u
                                                  + static_cast<Real>(12)*c[4]*u2));
      }
      return (val < min ? min : c[0] + nonLin_*(c[1]*u + c[2]*u2 + c[3]*u3 + c[4]*u4));
    }
  }

  void evaluateVelocity(std::vector<Real> &adv, const std::vector<Real> &x, const std::vector<Real> &z_param) const {
    if ( x[1] < xmid_ ) {
      const std::vector<Real> param = PDE<Real,DeviceType>::getParameter();
      Real p20 = (param.size() > 20 ? param[20] : static_cast<Real>(0.5));
      const Real min = static_cast<Real>(0.1)*xmid_;
      const Real max = static_cast<Real>(0.9)*xmid_;
      const Real x1  = static_cast<Real>(0.5)*((max-min)*p20 + (max+min));
      const Real mag = ((x[1] <  x1) ? x[1]/x1 : (xmid_-x[1])/(xmid_-x1));
      adv[0] = -z_param[0]*mag;
      adv[1] = static_cast<Real>(0);
    }
    else {
      adv[0] = static_cast<Real>(0);
      adv[1] = static_cast<Real>(0);
    }
  }

  Real evaluateRHS(const std::vector<Real> &x) const {
    return static_cast<Real>(0);
  }

  Real evaluateStefanBoltzmann(Real u, const std::vector<Real> &x,
                               int sideset, int locSideId,
                               int deriv = 0) const {
    const std::vector<Real> param = PDE<Real,DeviceType>::getParameter();
    // sig is the Stefan-Boltzmann constant
    // c1 is sig times the emissivity [0.5,1.5]
    // c2 is the ambient temperature away from the aluminum
    // c3 is the thermal convectivity of air (5), oil (40), and water (440)
    Real c1(0), c2(0), c3(0), sig(5.67e-8);
    if ( sideset == 2 ) {
      Real p21 = (param.size()>21 ? param[21] : static_cast<Real>(0.5));
      Real p22 = (param.size()>22 ? param[22] : static_cast<Real>(0.5));
      Real p23 = (param.size()>23 ? param[23] : static_cast<Real>(0.5));
      c1 = SBscale_ * sig * (static_cast<Real>(0.09) + static_cast<Real>(5.e-3) * p21);
      c2 = airTemp_             + static_cast<Real>(0.02*airTemp_) * p22;
      c3 = static_cast<Real>(5) + static_cast<Real>(0.5)           * p23;
    }
    else if ( sideset == 4 || sideset == 5 ) {
      Real p24 = (param.size()>24 ? param[24] : static_cast<Real>(0.5));
      Real p25 = (param.size()>25 ? param[25] : static_cast<Real>(0.5));
      Real p26 = (param.size()>26 ? param[26] : static_cast<Real>(0.5));
      c1 = SBscale_ * sig * (static_cast<Real>(0.09) + static_cast<Real>(5.e-3) * p24);
      c2 = engTemp_ + static_cast<Real>(0.2)*engTemp_ * p25;
      c3 = static_cast<Real>(40) + static_cast<Real>(2) * p26;
    }
    else if ( sideset == 1 ) {
      Real p27 = (param.size()>27 ? param[27] : static_cast<Real>(0.5));
      Real p28 = (param.size()>28 ? param[28] : static_cast<Real>(0.5));
      Real p29 = (param.size()>29 ? param[29] : static_cast<Real>(0.5));
      c1 = SBscale_ * sig * (static_cast<Real>(0.09) + static_cast<Real>(5.e-3) * p27);
      c2 = H2OTemp_ + static_cast<Real>(0.05)*H2OTemp_ * (p28 + static_cast<Real>(1));
      c3 = static_cast<Real>(440) + static_cast<Real>(20) * p29;
    }
    if ( deriv == 1 ) {
      return c1 * static_cast<Real>(4) * std::pow(u,3) + c3;
    }
    if ( deriv == 2 ) {
      return c1 * static_cast<Real>(4) * static_cast<Real>(3) * std::pow(u,2);
    }
    return c1 * (std::pow(u,4) - std::pow(c2,4)) + c3 * (u - c2);
  }

  Real evaluateRobin(Real u, Real z, const std::vector<Real> &x,
                     int sideset, int locSideId,
                     int deriv = 0, int component = 1) const {
    const std::vector<Real> param = PDE<Real,DeviceType>::getParameter();
    Real p30 = (param.size()>30 ? param[30] : static_cast<Real>(0.5));
    // c is the thermal convectivity of water (440)
    Real c = static_cast<Real>(440) + static_cast<Real>(20) * p30;
    if ( deriv == 1 ) {
      return (component==1) ? c : -c;
    }
    if ( deriv > 1 ) {
      return static_cast<Real>(0);
    }
    return c * (u - z);
  }

  Real evaluateRobin(Real u, std::vector<Real> &x,
                     int sideset, int locSideId,
                     int deriv = 0, int component = 1) const {
    const std::vector<Real> param = PDE<Real,DeviceType>::getParameter();
    Real p31 = (param.size()>31 ? param[31] : static_cast<Real>(0.5));
    // c is the thermal convectivity of water (440)
    Real c = static_cast<Real>(440) + static_cast<Real>(20) * p31;
    if ( deriv == 1 ) {
      return (component==1) ? c : static_cast<Real>(0);
    }
    if ( deriv > 1 ) {
      return static_cast<Real>(0);
    }
    return c * (u - static_cast<Real>(293));
  }
 
  /***************************************************************************/
  /************** COMPUTE PDE COEFFICIENTS AT DOFS ***************************/
  /***************************************************************************/
  void computeDiffusivity(scalar_view &kappa, const scalar_view u, int deriv = 0 ) const {
    // GET DIMENSIONS
    int c = fe_vol_->gradN().extent_int(0);
    int p = fe_vol_->gradN().extent_int(2);
    int d = fe_vol_->gradN().extent_int(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k)
          pt[k] = (fe_vol_->cubPts())(i,j,k);
        // Compute diffusivity
        kappa(i,j) = scale_*evaluateDiffusivity(u(i,j),pt,deriv);
      }
    }
  }

  void computeCoefficients(scalar_view &V, scalar_view &rhs,
       const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) const {
    // GET DIMENSIONS
    int c = fe_vol_->gradN().extent_int(0);
    int p = fe_vol_->gradN().extent_int(2);
    int d = fe_vol_->gradN().extent_int(3);
    std::vector<Real> pt(d), adv(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k)
          pt[k] = (fe_vol_->cubPts())(i,j,k);
        // Compute advection velocity field V
        if (z_param != ROL::nullPtr) { 
          evaluateVelocity(adv,pt,*z_param);
        }
        else {
          std::vector<Real> param = {advMag_};
          evaluateVelocity(adv,pt,param);
        }
        for (int k = 0; k < d; ++k) {
          V(i,j,k) = scale_*adv[k];
        }
        // Compute forcing term f
        rhs(i,j) = -scale_*evaluateRHS(pt);
      }
    }
  }

  void computeStefanBoltzmann(scalar_view &sb,
                              const scalar_view u,
                              int sideset,
                              int locSideId,
                              int deriv = 0) const {
    const int c = u.extent_int(0);
    const int p = u.extent_int(1);
    const int d = fe_vol_->gradN().extent_int(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (fe_bdry_[sideset][locSideId]->cubPts())(i,j,k);
        }
        sb(i,j) = evaluateStefanBoltzmann(u(i,j),pt,sideset,locSideId,deriv);
      }
    }
  }

  void computeRobin(scalar_view &robin,
                    const scalar_view u,
                    const scalar_view z,
                    int sideset,
                    int locSideId,
                    int deriv = 0,
                    int component = 1) const {
    const int c = u.extent_int(0);
    const int p = u.extent_int(1);
    const int d = fe_vol_->gradN().extent_int(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (fe_bdry_[sideset][locSideId]->cubPts())(i,j,k);
        }
        if (z != scalar_view()) {
          robin(i,j) = evaluateRobin(u(i,j),z(i,j),pt,sideset,locSideId,deriv,component);
        }
        else {
          robin(i,j) = evaluateRobin(u(i,j),pt,sideset,locSideId,deriv,component);
        }
      }
    }
  }

  /***************************************************************************/
  /************** EXTRACT COEFFICIENTS ON BOUNDARY ***************************/
  /***************************************************************************/
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

  /***************************************************************************/
  /************** COMPUTE LEAST SQUARES COEFFICIENTS *************************/
  /***************************************************************************/
  void getAluminumCoeff(std::vector<Real> &c, const std::vector<Real> &param) const {
    const std::vector<Real> Ta = {273, 300, 350, 400, 500, 600, 700, 800, 900};
    const std::vector<Real> Ka = {236, 237, 240, 240, 237, 232, 226, 220, 213};
    Teuchos::LAPACK<int,Real> lapack;
    const char trans = 'N';
    const int m = Ta.size();
    const int n = 5;
    const int nrhs = 1;
    const int lda = m;
    const int ldb = m;
    std::vector<Real> A(m*n,1);
    std::vector<Real> b(m,1);
    for (int i = 0; i < m; ++i) {
      b[i] = Ka[i] + param[i];
      for (int j = 0; j < n; ++j) {
        A[j*m + i] = std::pow(Ta[i],j);
      }
    }
    const int lwork = n + m;
    std::vector<Real> work(lwork,0);
    int info;
    lapack.GELS(trans,m,n,nrhs,&A[0],lda,&b[0],ldb,&work[0],lwork,&info);
    c.clear(); c.resize(n,0);
    for (int i = 0; i < n; ++i) {
      c[i] = b[i];
    }
  }

  void getWaterCoeff(std::vector<Real> &c, const std::vector<Real> &param) const {
    const std::vector<Real> Tw = {270, 280, 290, 300, 310, 320, 330, 340, 350, 370, 400};
    const std::vector<Real> Kw = {0.5551, 0.5818, 0.5918, 0.6084, 0.6233, 0.6367, 0.6485, 0.6587, 0.6673, 0.6797, 0.6864};
    Teuchos::LAPACK<int,Real> lapack;
    const char trans = 'N';
    const int m = Tw.size();
    const int n = 3;
    const int nrhs = 1;
    const int lda = m;
    const int ldb = m;
    std::vector<Real> A(m*n,1);
    std::vector<Real> b(m,1);
    for (int i = 0; i < m; ++i) {
      b[i] = Kw[i] + param[i];
      for (int j = 0; j < n; ++j) {
        A[j*m + i] = std::pow(Tw[i],j);
      }
    }
    const int lwork = n + m;
    std::vector<Real> work(lwork,0);
    int info;
    lapack.GELS(trans,m,n,nrhs,&A[0],lda,&b[0],ldb,&work[0],lwork,&info);
    c.clear(); c.resize(n,0);
    for (int i = 0; i < n; ++i) {
      c[i] = b[i];
    }
  }

}; // PDE_stefan_boltzmann

#endif
