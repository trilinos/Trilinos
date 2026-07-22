// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde.hpp
    \brief Implements the local PDE interface for the Poisson-Boltzmann control problem.
*/

#ifndef PDE_POISSON_BOLTZMANN_EX02K_HPP
#define PDE_POISSON_BOLTZMANN_EX02K_HPP

#include "../TOOLS/pdeK.hpp"
#include "../TOOLS/feK.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real, class DeviceType>
class PDE_Poisson_Boltzmann_ex02 : public PDE<Real, DeviceType> {
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
  ROL::Ptr<Intrepid2::Cubature<DeviceType, Real, Real> > cellCub_;
  ROL::Ptr<Intrepid2::Cubature<DeviceType, Real, Real> > bdryCub_;
  // Cell node information
  scalar_view volCellNodes_;
  std::vector<std::vector<scalar_view>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<fe_type> fe_vol_;
  std::vector<std::vector<ROL::Ptr<fe_type>>> fe_bdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fidx_;

  Real wx_, wy_, dx_;
  Real uScale_, vScale_;
  bool useRobin_;
  Real robinCoeff_;

public:
  PDE_Poisson_Boltzmann_ex02(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("Problem").get("Order of FE discretization",1);
    if (basisOrder == 1) {
      basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType, Real, Real>>();
    }
    else if (basisOrder == 2) {
      basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<DeviceType, Real, Real>>();
    }
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();                  // get the cell type from any basis
    Intrepid2::DefaultCubatureFactory cubFactory;                                      // create cubature factory
    int cubDegree = parlist.sublist("PDE Poisson Boltzmann").get("Cubature Degree",2); // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create<DeviceType, Real, Real>(cellType, cubDegree);         // create default cubature

    int d = cellType.getDimension();
    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2); // set cubature degree, e.g., 2
    bdryCub_ = cubFactory.create<DeviceType, Real, Real>(bdryCellType, bdryCubDegree);

    wx_      = parlist.sublist("Geometry").get("Width",0.6);
    wy_      = parlist.sublist("Geometry").get("Height",0.2);
    dx_      = wx_/static_cast<Real>(6);
    uScale_  = parlist.sublist("Problem").get("Electron Density",1.0);
    vScale_  = parlist.sublist("Problem").get("Hole Density",1.0);

    useRobin_ = parlist.sublist("Problem").get("Use Robin Conditions",false);
    robinCoeff_ = parlist.sublist("Problem").get("Robin Coefficient",1e2);
  }

  void residual(scalar_view & res,
                const scalar_view u_coeff,
                const scalar_view z_coeff = scalar_view(),
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN().extent_int(0);
    const int f = fe_vol_->gradN().extent_int(1);
    const int p = fe_vol_->gradN().extent_int(2);
    const int d = fe_vol_->gradN().extent_int(3);
    // INITIALIZE RESIDUAL
    res = scalar_view("", c, f);
    // COMPUTE PDE COEFFICIENTS
    scalar_view lambda2("lambda2", c, p), delta2("delta2", c, p), scale("scale", c, p);
    computeCoefficients(lambda2,delta2,scale);
    // ADD STIFFNESS TERM TO RESIDUAL
    scalar_view gradU_eval("gradU_eval", c, p, d);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    scalar_view lambda2_gradU_eval("lambda2_gradU_eval", c, p, d);
    fst::tensorMultiplyDataData(lambda2_gradU_eval, lambda2, gradU_eval);
    fst::integrate(res, lambda2_gradU_eval, fe_vol_->gradNdetJ(), false);
    // ADD NONLINEAR TERM TO RESIDUAL
    scalar_view valU_eval("valU_eval", c, p);
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    scalar_view phi_valU_eval("phi_valU_eval", c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        phi_valU_eval(i,j) = delta2(i,j)*(uScale_*std::exp((valU_eval)(i,j)) - vScale_*std::exp(-(valU_eval)(i,j)));
      }
    }
    fst::integrate(res, phi_valU_eval, fe_vol_->NdetJ(), true);
    // ADD CONTROL TERM TO RESIDUAL
    scalar_view valZ_eval("valZ_eval", c, p);
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    scalar_view valZ_scal("valZ_scal", c, p);
    fst::scalarMultiplyDataData(valZ_scal, scale, valZ_eval);
    fst::integrate(res, valZ_scal, fe_vol_->NdetJ(), true);
    // APPLY DIRICHLET CONDITIONS
    std::vector<std::vector<scalar_view>> bdryCellDofValues;
    computeDirichlet(bdryCellDofValues);
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      if ( i==1 || i==2 || i==3 ) {
        int numLocalSideIds = bdryCellLocIds_[i].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          int numCellsSide = bdryCellLocIds_[i][j].size();
          if (!useRobin_) {
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                res(cidx,fidx_[j][l])
                  = u_coeff(cidx,fidx_[j][l]) - (bdryCellDofValues[i][j])(k,l);
              }
            }
          }
          else {
            const int numCubPerSide = bdryCub_->getNumPoints();
            if (numCellsSide) {
              // Get U and Z coefficients on Robin boundary
              scalar_view u_coeff_bdry, z_coeff_bdry;
              u_coeff_bdry = getBoundaryCoeff(u_coeff, i, j);
              z_coeff_bdry = getBoundaryCoeff(z_coeff, i, j);
              // Evaluate U and Z on FE basis
              scalar_view valU_eval_bdry, valZ_eval_bdry;
              valU_eval_bdry = scalar_view("valU_eval_bdry", numCellsSide, numCubPerSide);
              valZ_eval_bdry = scalar_view("valZ_eval_bdry", numCellsSide, numCubPerSide);
              fe_bdry_[i][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
              fe_bdry_[i][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
              // Compute Robin residual
              scalar_view robinRes("robinRes", numCellsSide, f);
              scalar_view robinVal("robinVal", numCellsSide, numCubPerSide);
              computeRobin(robinVal,valU_eval_bdry,valZ_eval_bdry,i,j,0);
              fst::integrate(robinRes, robinVal, (fe_bdry_[i][j])->NdetJ(), false);
              // Add Stefan-Boltzmann residual to volume residual
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < f; ++l) {
                  res(cidx,l) += robinRes(k,l);
                }
              }
            }
          }
        }
      }
    }
  }

  void Jacobian_1(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN().extent_int(0);
    const int f = fe_vol_->gradN().extent_int(1);
    const int p = fe_vol_->gradN().extent_int(2);
    const int d = fe_vol_->gradN().extent_int(3);
    // INITIALIZE JACOBIAN
    jac = scalar_view("jac1", c, f, f);
    // COMPUTE PDE COEFFICIENTS
    scalar_view lambda2("lambda2", c, p), delta2("delta2", c, p), scale("scale", c, p);
    computeCoefficients(lambda2, delta2, scale);
    // ADD STIFFNESS TERM TO JACOBIAN
    scalar_view lambda2_gradN_eval("lambda2_gradN_eval", c, f, p, d);
    fst::tensorMultiplyDataField(lambda2_gradN_eval, lambda2, fe_vol_->gradN());
    // COMPUTE STIFFNESS TERM
    fst::integrate(jac, lambda2_gradN_eval, fe_vol_->gradNdetJ(), false);
    // ADD NONLINEAR TERM
    scalar_view valU_eval("valU_eval", c, p);
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    scalar_view dphi_valU_eval("dphi_valU_eval", c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        dphi_valU_eval(i,j) = delta2(i,j)*(uScale_*std::exp((valU_eval)(i,j)) + vScale_*std::exp(-(valU_eval)(i,j)));
      }
    }
    scalar_view NexpU("NexpU", c, f, p);
    fst::scalarMultiplyDataField(NexpU, dphi_valU_eval, fe_vol_->N());
    fst::integrate(jac, NexpU, fe_vol_->NdetJ(), true);
    // APPLY DIRICHLET CONDITIONS
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      if ( i==1 || i==2 || i==3 ) {
        int numLocalSideIds = bdryCellLocIds_[i].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          int numCellsSide = bdryCellLocIds_[i][j].size();
          if (!useRobin_) {
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
                for (int m = 0; m < f; ++m) {
                  jac(cidx, fidx_[j][l], m) = static_cast<Real>(0);
                }
                jac(cidx, fidx_[j][l], fidx_[j][l]) = static_cast<Real>(1);
              }
            }
          }
          else {
            const int numCubPerSide = bdryCub_->getNumPoints();
            if (numCellsSide) {
              // Get U and Z coefficients on Robin boundary
              scalar_view u_coeff_bdry, z_coeff_bdry;
              u_coeff_bdry = getBoundaryCoeff(u_coeff, i, j);
              z_coeff_bdry = getBoundaryCoeff(z_coeff, i, j);
              // Evaluate U and Z on FE basis
              scalar_view valU_eval_bdry, valZ_eval_bdry;
              valU_eval_bdry = scalar_view("valU_eval_bdry", numCellsSide, numCubPerSide);
              valZ_eval_bdry = scalar_view("valZ_eval_bdry", numCellsSide, numCubPerSide);
              fe_bdry_[i][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
              fe_bdry_[i][j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
              // Compute Robin residual
              scalar_view robinVal_N("robinVal_N", numCellsSide, f, numCubPerSide);
              scalar_view robinJac("robinJac", numCellsSide, f, f);
              scalar_view robinVal("robinVal", numCellsSide, numCubPerSide);
              computeRobin(robinVal, valU_eval_bdry, valZ_eval_bdry, i, j, 1, 1);
              fst::scalarMultiplyDataField(robinVal_N, robinVal, (fe_bdry_[i][j])->N());
              fst::integrate(robinJac, robinVal_N, (fe_bdry_[i][j])->NdetJ(), false);
              // Add Stefan-Boltzmann residual to volume residual
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < f; ++l) {
                  for (int m = 0; m < f; ++m) {
                    jac(cidx, l, m) += robinJac(k, l, m);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  void Jacobian_2(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN().extent_int(0);
    const int f = fe_vol_->gradN().extent_int(1);
    const int p = fe_vol_->gradN().extent_int(2);
    // INITIALIZE JACOBIAN
    jac = scalar_view("jac2", c, f, f);
    // COMPUTE PDE COEFFICIENTS
    scalar_view lambda2("lambda2", c, p), delta2("delta2", c, p), scale("scale", c, p);
    computeCoefficients(lambda2,delta2,scale);
    // ADD CONTROL TERM
    scalar_view valN_scal("valN_scal", c, f, p);
    fst::scalarMultiplyDataField(valN_scal, scale, fe_vol_->N());
    fst::integrate(jac, valN_scal, fe_vol_->NdetJ(), false);
    // APPLY DIRICHLET CONDITIONS
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      if ( i==1 || i==2 || i==3 ) {
        int numLocalSideIds = bdryCellLocIds_[i].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          int numCellsSide = bdryCellLocIds_[i][j].size();
          if (!useRobin_) {
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
                for (int m = 0; m < f; ++m) {
                  jac(cidx, fidx_[j][l], m) = static_cast<Real>(0);
                }
              }
            }
          }
          else {
            const int numCubPerSide = bdryCub_->getNumPoints();
            if (numCellsSide) {
              // Get U and Z coefficients on Robin boundary
              scalar_view u_coeff_bdry, z_coeff_bdry;
              u_coeff_bdry = getBoundaryCoeff(u_coeff, i, j);
              z_coeff_bdry = getBoundaryCoeff(z_coeff, i, j);
              // Evaluate U and Z on FE basis
              scalar_view valU_eval_bdry, valZ_eval_bdry;
              valU_eval_bdry = scalar_view("valU_eval_bdry", numCellsSide, numCubPerSide);
              valZ_eval_bdry = scalar_view("valZ_eval_bdry", numCellsSide, numCubPerSide);
              (fe_bdry_[i][j])->evaluateValue(valU_eval_bdry, u_coeff_bdry);
              (fe_bdry_[i][j])->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
              // Compute Robin residual
              scalar_view robinVal_N("robinVal_N", numCellsSide, f, numCubPerSide);
              scalar_view robinJac("robinJac", numCellsSide, f, f);
              scalar_view robinVal("robinVal", numCellsSide, numCubPerSide);
              computeRobin(robinVal, valU_eval_bdry, valZ_eval_bdry, i, j, 1, 2);
              fst::scalarMultiplyDataField(robinVal_N, robinVal, (fe_bdry_[i][j])->N());
              fst::integrate(robinJac, robinVal_N, (fe_bdry_[i][j])->NdetJ(), false);
              // Add Stefan-Boltzmann residual to volume residual
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < f; ++l) {
                  for (int m = 0; m < f; ++m) {
                    jac(cidx, l, m) += robinJac(k, l, m);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  void Hessian_11(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN().extent_int(0);
    const int f = fe_vol_->gradN().extent_int(1);
    const int p = fe_vol_->gradN().extent_int(2);
    // INITIALIZE HESSIAN
    hess = scalar_view("hess11", c, f, f);
    // COMPUTE PDE COEFFICIENTS
    scalar_view lambda2("lambda2", c, p), delta2("delta2", c, p), scale("scale", c, p);
    computeCoefficients(lambda2,delta2,scale);
    // APPLY DIRICHLET CONDITIONS
    scalar_view l_dbc("l_dbc", c, p);
    Kokkos::deep_copy(l_dbc, l_coeff);
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      if ( i==1 || i==2 || i==3 ) {
        int numLocalSideIds = bdryCellLocIds_[i].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          int numCellsSide = bdryCellLocIds_[i][j].size();
          int numBdryDofs = fidx_[j].size();
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
            for (int l = 0; l < numBdryDofs; ++l) {
              l_dbc(cidx, fidx_[j][l]) = static_cast<Real>(0);
            }
          }
        }
      }
    }
    // COMPUTE NONLINEAR TERM
    scalar_view valU_eval("valU_eval", c, p);
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    scalar_view valL_eval("valL_eval", c, p);
    fe_vol_->evaluateValue(valL_eval, l_dbc);
    scalar_view d2phi_valU_eval("d2phi_valU_eval", c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        d2phi_valU_eval(i,j) = valL_eval(i,j)*delta2(i,j)
                              *(uScale_*std::exp(valU_eval(i,j))-vScale_*std::exp(-valU_eval(i,j)));
      }
    }
    scalar_view NLexpU("NLexpU", c, f, p);
    fst::scalarMultiplyDataField(NLexpU, d2phi_valU_eval, fe_vol_->N());
    fst::integrate(hess, NLexpU, fe_vol_->NdetJ(), false);
  }

  void Hessian_12(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Poisson_Boltzmann_ex02:Hessian_12: Hessian is zero.");
  }

  void Hessian_21(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Poisson_Boltzmann_ex02:Hessian_21: Hessian is zero.");
  }

  void Hessian_22(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Poisson_Boltzmann_ex02:Hessian_22: Hessian is zero.");
  }

  void RieszMap_1(scalar_view & riesz) override {
    // GET DIMENSIONS
    const int c = fe_vol_->N().extent_int(0);
    const int f = fe_vol_->N().extent_int(1);
    // INITIALIZE RIESZ
    riesz = scalar_view("riesz1", c, f, f);
    Kokkos::deep_copy(riesz, fe_vol_->stiffMat());
    rst::add(riesz, fe_vol_->massMat());
  }

  void RieszMap_2(scalar_view & riesz) override {
    // GET DIMENSIONS
    const int c = fe_vol_->N().extent_int(0);
    const int f = fe_vol_->N().extent_int(1);
    // INITIALIZE RIESZ
    riesz = scalar_view("riesz2", c, f, f);
    Kokkos::deep_copy(riesz, fe_vol_->stiffMat());
    rst::add(riesz, fe_vol_->massMat());
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
    fe_vol_ = ROL::makePtr<fe_type>(volCellNodes_, basisPtr_, cellCub_);
    // Set local boundary DOFs.
    fidx_ = fe_vol_->getBoundaryDofs();
    // Construct boundary FEs
    const int numSidesets = bdryCellNodes.size();
    fe_bdry_.resize(numSidesets);
    for(int i = 0; i < numSidesets; ++i) {
      int numLocSides = bdryCellNodes[i].size();
      fe_bdry_[i].resize(numLocSides);
      for (int j = 0; j < numLocSides; ++j) {
        if (bdryCellNodes_[i][j] != scalar_view()) {
          fe_bdry_[i][j] = ROL::makePtr<fe_type>(bdryCellNodes_[i][j], basisPtr_, bdryCub_, j);
        }
      }
    }
  }

  const ROL::Ptr<fe_type> getFE(void) const {
    return fe_vol_;
  }

  const std::vector<std::vector<ROL::Ptr<fe_type>>> getBdryFE(void) const {
    return fe_bdry_;
  }

  const scalar_view getCellNodes(void) const {
    return volCellNodes_;
  }

  const std::vector<std::vector<std::vector<int>>> getBdryCellLocIds(void) const {
    return bdryCellLocIds_;
  }

private:

  Real evaluateLambda2(const std::vector<Real> &x) const {
    const std::vector<Real> param = PDE<Real, DeviceType>::getParameter();
    Real p0 = (param.size()>0 ? param[0] : static_cast<Real>(-1.5));
    return static_cast<Real>(2.5)*std::pow(static_cast<Real>(10), p0);
  }

  Real evaluateDelta2(const std::vector<Real> &x) const {
    const std::vector<Real> param = PDE<Real, DeviceType>::getParameter();
    Real p1 = (param.size()>1 ? param[1] : static_cast<Real>(-0.5));
    return static_cast<Real>(1.45)*std::pow(static_cast<Real>(10), p1);
  }

  Real evaluateScale(const std::vector<Real> &x) const {
    return 1;
//    const std::vector<Real> param = PDE<Real>::getParameter();
//    const Real Y1 = static_cast<Real>(1) + static_cast<Real>(5)*(wy_-x[1])*param[2];
//    const Real Y2 = static_cast<Real>(1) + static_cast<Real>(5)*(wy_-x[1])*param[3];
//    const Real X1 = phi(x[0]);
//    const Real X2 = static_cast<Real>(1)-X1;
//    return Y1*X1 + Y2*X2;
  }

  Real phi(const Real x) const {
    const Real zero(0), one(1), two(2);
    const Real eps = std::sqrt(ROL::ROL_EPSILON<Real>());
    const Real m = one/(wx_-two*dx_);
    const Real v = (x < dx_+eps ? zero
                   : (x > wx_-dx_-eps ? one
                     : m*(x-dx_)));
    return v;
  }

  void computeCoefficients(scalar_view & lambda2,
                           scalar_view & delta2,
                           scalar_view & scale) const {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN().extent_int(0);
    const int p = fe_vol_->gradN().extent_int(2);
    const int d = fe_vol_->gradN().extent_int(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (fe_vol_->cubPts())(i,j,k);
        }
        // Compute scaled Debye length
        lambda2(i,j) = evaluateLambda2(pt);
        delta2(i,j)  = evaluateDelta2(pt);
        scale(i,j)   = -evaluateScale(pt);
      }
    }
  }

  Real evaluateRobin(const Real u, const Real z, const std::vector<Real> &x,
                     const int sideset, const int locSideId,
                     const int deriv = 0, const int component = 1) const {
    const Real four(4), two(2), one(1);
    Real C(0);
    if (sideset==1 || sideset==3) {
      C = static_cast<Real>(1);
    }
    if (sideset==2) {
      C = static_cast<Real>(0.3);
    }
    Real f0 = std::log((C+std::sqrt(C*C+four))/two);
    Real f1 = one/std::sqrt(C*C+four);
    if ( deriv == 1 ) {
      return (component==1) ? robinCoeff_ : -robinCoeff_ * f1;
    }
    if ( deriv > 1 ) {
      return static_cast<Real>(0);
    }
    return robinCoeff_ * (u - (f0 + f1*(z-C)));
  }

  void computeRobin(scalar_view & robin,
                    const scalar_view u,
                    const scalar_view z,
                    const int sideset,
                    const int locSideId,
                    const int deriv = 0,
                    const int component = 1) const {
    const int c = u.extent_int(0);
    const int p = u.extent_int(1);
    const int d = fe_vol_->gradN().extent_int(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (fe_bdry_[sideset][locSideId]->cubPts())(i, j, k);
        }
        robin(i,j) = evaluateRobin(u(i,j), z(i,j), pt, sideset, locSideId, deriv, component);
      }
    }
  }

  Real evaluateDirichlet(const std::vector<Real> & coords, int sideset, int locSideId) const {
    Real val(0);
    if (sideset==1) {
      val = static_cast<Real>(0.436);
    }
    if (sideset==2) {
      val = static_cast<Real>(0.407);
    }
    if (sideset==3) {
      val = static_cast<Real>(0.436);
    }
    return val;
  }

  void computeDirichlet(std::vector<std::vector<scalar_view>> & bdryCellDofValues) const {
    // Compute Dirichlet values at DOFs.
    int d = basisPtr_->getBaseCellTopology().getDimension();
    int numSidesets = bdryCellLocIds_.size();
    bdryCellDofValues.resize(numSidesets);
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellDofValues[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
        int f = basisPtr_->getCardinality();
        bdryCellDofValues[i][j] = scalar_view("bdryCellDofValues", c, f);
        scalar_view coords("coords", c, f, d);
        if (c > 0) {
          fe_vol_->computeDofCoords(coords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<f; ++l) {
            std::vector<Real> dofpoint(d);
            for (int m=0; m<d; ++m) {
              dofpoint[m] = coords(k, l, m);
            }
            (bdryCellDofValues[i][j])(k, l) = evaluateDirichlet(dofpoint, i, j);
          }
        }
      }
    }
  }

  scalar_view getBoundaryCoeff(
      const scalar_view cell_coeff,
      int sideSet, int cell) const {
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

}; // PDE_Poisson_Boltzmann

template <class Real, class DeviceType>
class PDE_Doping : public PDE<Real, DeviceType> {
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
  ROL::Ptr<fe_type> fe_vol_;
  std::vector<std::vector<ROL::Ptr<fe_type>>> fe_bdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fidx_;
  Real a_, b_;

public:
  PDE_Doping(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("Problem").get("Order of FE discretization",1);
    if (basisOrder == 1) {
      basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType, Real, Real>>();
    }
    else if (basisOrder == 2) {
      basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<DeviceType, Real, Real>>();
    }
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();                  // get the cell type from any basis
    Intrepid2::DefaultCubatureFactory cubFactory;                                      // create cubature factory
    int cubDegree = parlist.sublist("PDE Poisson Boltzmann").get("Cubature Degree",2); // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create<DeviceType, Real, Real>(cellType, cubDegree);         // create default cubature

    int d = cellType.getDimension();
    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2); // set cubature degree, e.g., 2
    bdryCub_ = cubFactory.create<DeviceType, Real, Real>(bdryCellType, bdryCubDegree);
    a_  = parlist.sublist("Problem").get("Desired Lower Doping Value", 0.3);
    b_  = parlist.sublist("Problem").get("Desired Upper Doping Value", 1.0);
  }

  void residual(scalar_view & res,
                const scalar_view u_coeff,
                const scalar_view z_coeff = scalar_view(),
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    const std::vector<Real> param = PDE<Real, DeviceType>::getParameter();
    Real cond = (param.size()>2 ? param[2] : static_cast<Real>(2e-5));
    // GET DIMENSIONS
    const int c = fe_vol_->gradN().extent_int(0);
    const int f = fe_vol_->gradN().extent_int(1);
    const int p = fe_vol_->gradN().extent_int(2);
    const int d = fe_vol_->gradN().extent_int(3);
    // INITIALIZE RESIDUAL
    res = scalar_view("res", c, f);
    // EVALUATE STATE AND CONTROL OF FEM BASIS
    scalar_view valU_eval, valZ_eval, gradU_eval;
    valU_eval  = scalar_view("valU_eval", c, p);
    valZ_eval  = scalar_view("valZ_eval", c, p);
    gradU_eval = scalar_view("gradU_eval", c, p, d);
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    // ADD STIFFNESS TERM TO RESIDUAL
    scalar_view lambda2_gradU_eval("lambda2_gradU_eval", c, p, d);
    rst::scale(lambda2_gradU_eval, gradU_eval, cond);
    fst::integrate(res, lambda2_gradU_eval, fe_vol_->gradNdetJ(), false);
    // ADD REACTION TERM TO RESIDUAL
    fst::integrate(res, valU_eval, fe_vol_->NdetJ(), true);
    // ADD CONTROL TERM TO RESIDUAL
    scalar_view valZ_scal("valZ_scal", c, p);
    rst::scale(valZ_scal, valZ_eval, static_cast<Real>(-1));
    fst::integrate(res, valZ_scal, fe_vol_->NdetJ(), true);
    // APPLY DIRICHLET CONDITIONS
    std::vector<std::vector<scalar_view>> bdryCellDofValues;
    computeDirichlet(bdryCellDofValues);
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      if ( i==1 || i==2 || i==3 ) {
        int numLocalSideIds = bdryCellLocIds_[i].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          int numCellsSide = bdryCellLocIds_[i][j].size();
          int numBdryDofs = fidx_[j].size();
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
            for (int l = 0; l < numBdryDofs; ++l) {
              res(cidx, fidx_[j][l]) = u_coeff(cidx, fidx_[j][l]) - (bdryCellDofValues[i][j])(k,l);
            }
          }
        }
      }
    }
  }

  void Jacobian_1(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    const std::vector<Real> param = PDE<Real, DeviceType>::getParameter();
    Real cond = (param.size()>2 ? param[2] : static_cast<Real>(2e-5));
    // GET DIMENSIONS
    const int c = fe_vol_->gradN().extent_int(0);
    const int f = fe_vol_->gradN().extent_int(1);
    const int p = fe_vol_->gradN().extent_int(2);
    const int d = fe_vol_->gradN().extent_int(3);
    // INITIALIZE JACOBIAN
    jac = scalar_view("jac1", c, f, f);
    // ADD STIFFNESS TERM TO JACOBIAN
    scalar_view lambda2_gradN_eval("lambda2_gradN_eval", c, f, p, d);
    rst::scale(lambda2_gradN_eval, fe_vol_->gradN(), cond);
    fst::integrate(jac, lambda2_gradN_eval, fe_vol_->gradNdetJ(), false);
    // ADD REACTION TERM TO JACOBIAN
    fst::integrate(jac, fe_vol_->N(), fe_vol_->NdetJ(), true);
    // APPLY DIRICHLET CONDITIONS
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      if ( i==1 || i==2 || i==3 ) {
        int numLocalSideIds = bdryCellLocIds_[i].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          int numCellsSide = bdryCellLocIds_[i][j].size();
          int numBdryDofs = fidx_[j].size();
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
            for (int l = 0; l < numBdryDofs; ++l) {
              //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
              for (int m = 0; m < f; ++m) {
                jac(cidx, fidx_[j][l], m) = static_cast<Real>(0);
              }
              jac(cidx, fidx_[j][l], fidx_[j][l]) = static_cast<Real>(1);
            }
          }
        }
      }
    }
  }

  void Jacobian_2(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // GET DIMENSIONS
    const int c = fe_vol_->gradN().extent_int(0);
    const int f = fe_vol_->gradN().extent_int(1);
    const int p = fe_vol_->gradN().extent_int(2);
    // INITIALIZE JACOBIAN
    jac = scalar_view("jac2", c, f, f);
    // ADD CONTROL TERM
    scalar_view valN_scal("valN_scal", c, f, p);
    rst::scale(valN_scal, fe_vol_->N(), static_cast<Real>(-1));
    fst::integrate(jac, valN_scal, fe_vol_->NdetJ(), false);
    // APPLY DIRICHLET CONDITIONS
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      if ( i==1 || i==2 || i==3 ) {
        int numLocalSideIds = bdryCellLocIds_[i].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          int numCellsSide = bdryCellLocIds_[i][j].size();
          int numBdryDofs = fidx_[j].size();
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
            for (int l = 0; l < numBdryDofs; ++l) {
              //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
              for (int m = 0; m < f; ++m) {
                jac(cidx, fidx_[j][l], m) = static_cast<Real>(0);
              }
            }
          }
        }
      }
    }
  }

  void Hessian_11(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Poisson_Boltzmann_ex02:Hessian_11: Hessian is zero.");
  }

  void Hessian_12(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Poisson_Boltzmann_ex02:Hessian_12: Hessian is zero.");
  }

  void Hessian_21(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Poisson_Boltzmann_ex02:Hessian_21: Hessian is zero.");
  }

  void Hessian_22(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Poisson_Boltzmann_ex02:Hessian_22: Hessian is zero.");
  }

  void RieszMap_1(scalar_view & riesz) override {
    // GET DIMENSIONS
    const int c = fe_vol_->N().extent_int(0);
    const int f = fe_vol_->N().extent_int(1);
    // INITIALIZE RIESZ
    riesz = scalar_view("riesz1", c, f, f);
    Kokkos::deep_copy(riesz, fe_vol_->stiffMat());
    rst::add(riesz, fe_vol_->massMat());
  }

  void RieszMap_2(scalar_view & riesz) override {
    // GET DIMENSIONS
    const int c = fe_vol_->N().extent_int(0);
    const int f = fe_vol_->N().extent_int(1);
    // INITIALIZE RIESZ
    riesz = scalar_view("riesz2", c, f, f);
    Kokkos::deep_copy(riesz, fe_vol_->stiffMat());
    rst::add(riesz, fe_vol_->massMat());
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
    fe_vol_ = ROL::makePtr<fe_type>(volCellNodes_, basisPtr_, cellCub_);
    // Set local boundary DOFs.
    fidx_ = fe_vol_->getBoundaryDofs();
    // Construct boundary FEs
    const int numSidesets = bdryCellNodes.size();
    fe_bdry_.resize(numSidesets);
    for(int i = 0; i < numSidesets; ++i) {
      int numLocSides = bdryCellNodes[i].size();
      fe_bdry_[i].resize(numLocSides);
      for (int j = 0; j < numLocSides; ++j) {
        if (bdryCellNodes_[i][j] != scalar_view()) {
          fe_bdry_[i][j] = ROL::makePtr<fe_type>(bdryCellNodes_[i][j], basisPtr_, bdryCub_, j);
        }
      }
    }
  }

  const ROL::Ptr<fe_type> getFE(void) const {
    return fe_vol_;
  }

  const std::vector<std::vector<fe_type>> getBdryFE(void) const {
    return fe_bdry_;
  }

  const scalar_view getCellNodes(void) const {
    return volCellNodes_;
  }

  const std::vector<std::vector<std::vector<int>>> getBdryCellLocIds(void) const {
    return bdryCellLocIds_;
  }

private:

  Real evaluateDirichlet(const std::vector<Real> & coords, int sideset, int locSideId) const {
    Real val(0);
    if (sideset==1) {
      val = a_;
    }
    if (sideset==2) {
      val = b_;
    }
    if (sideset==3) {
      val = a_;
    }
    return val;
  }

  void computeDirichlet(std::vector<std::vector<scalar_view>> & bdryCellDofValues) const {
    // Compute Dirichlet values at DOFs.
    int d = basisPtr_->getBaseCellTopology().getDimension();
    int numSidesets = bdryCellLocIds_.size();
    bdryCellDofValues.resize(numSidesets);
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellDofValues[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
        int f = basisPtr_->getCardinality();
        bdryCellDofValues[i][j] = scalar_view("bdryCellDofValues", c, f);
        scalar_view coords("coords", c, f, d);
        if (c > 0) {
          fe_vol_->computeDofCoords(coords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<f; ++l) {
            std::vector<Real> dofpoint(d);
            for (int m=0; m<d; ++m) {
              dofpoint[m] = coords(k, l, m);
            }
            (bdryCellDofValues[i][j])(k, l) = evaluateDirichlet(dofpoint, i, j);
          }
        }
      }
    }
  }

}; // PDE_Poisson_Boltzmann

#endif
