// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde_poisson_topOpt.hpp
    \brief Implements the local PDE interface for the Poisson
           topology optimization problem.
*/

#ifndef PDE_POISSON_TOPOPTK_HPP
#define PDE_POISSON_TOPOPTK_HPP

#include "../../TOOLS/pdeK.hpp"
#include "../../TOOLS/feK.hpp"

#include "Intrepid2_HVOL_C0_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_RealSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real, class DeviceType>
class PDE_Poisson_TopOpt : public PDE<Real,DeviceType> {
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
  std::vector<std::vector<scalar_view>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<fe_type> fe_vol_, fe_dens_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<scalar_view>> bdryCellDofValues_;
  // Force function evaluated at cubature points
  scalar_view force_eval_, nforce_eval_;
  // Inputs
  Real minConductivity_, SIMPpower_;
  bool getFields2called_;

  Real ForceFunc(const std::vector<Real> &x) const {
    return static_cast<Real>(0.01);
  }

  Real dirichletFunc(const std::vector<Real> & coords, const int sideset, const int locSideId) const {
    return static_cast<Real>(0);
  }

public:
  PDE_Poisson_TopOpt(Teuchos::ParameterList &parlist) : getFields2called_(false) {
    // Finite element fields.
    int basisOrder     = parlist.sublist("Problem").get("Order of FE discretization",1);
    int basisOrderDens = parlist.sublist("Problem").get("Density Basis Order",0);
    TEUCHOS_TEST_FOR_EXCEPTION(basisOrder > 2 || basisOrder < 1, std::invalid_argument,
      ">>> PDE-OPT/topopt/poisson/pde_poisson_topOpt.hpp: Basis order is not 1 or 2!");
    TEUCHOS_TEST_FOR_EXCEPTION(basisOrderDens > 1 || basisOrderDens < 0, std::invalid_argument,
      ">>> PDE-OPT/topopt/poisson/pde_poisson_topOpt.hpp: Basis order is not 0 or 1!");
    if (basisOrder == 1)
      basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType,Real,Real>>();
    else if (basisOrder == 2)
      basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<DeviceType,Real,Real>>();
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();        // get the cell type from any basis
    if (basisOrderDens == 1)
      basisPtrDens_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType,Real,Real>>();
    else
      basisPtrDens_ = ROL::makePtr<Intrepid2::Basis_HVOL_C0_FEM<DeviceType,Real,Real>>(cellType);
    basisPtrs_.clear();     basisPtrs_.push_back(basisPtr_);
    basisPtrsDens_.clear(); basisPtrsDens_.push_back(basisPtrDens_); // Density component
    // Quadrature rules.
    Intrepid2::DefaultCubatureFactory cubFactory;                            // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree",2);     // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create<DeviceType,Real,Real>(cellType, cubDegree); // create default cubature

    minConductivity_ = parlist.sublist("Problem").get("Minimum Conductivity",1.e-3);
    SIMPpower_       = parlist.sublist("Problem").get("SIMP Power",3.0);
  }

  void residual(scalar_view & res,
                const scalar_view u_coeff,
                const scalar_view z_coeff = scalar_view(),
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    const Real one(1);
    // Get dimensions
    const int c = fe_vol_->gradN().extent_int(0);
    const int f = fe_vol_->gradN().extent_int(1);
    const int p = fe_vol_->gradN().extent_int(2);
    const int d = fe_vol_->gradN().extent_int(3);
    // Initialize residual storage
    res = scalar_view("res",c,f);
    // Temporary storage
    scalar_view valZ_eval("valZ_eval",c,p);
    scalar_view gradU_eval("gradU_eval",c,p,d);
    scalar_view KgradU("KgradU",c,p,d);
    // Build SIMP density at cubature points
    fe_dens_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        valZ_eval(i,j) = minConductivity_
          + (one - minConductivity_) * std::pow(valZ_eval(i,j),SIMPpower_);
      }
    }
    // Build flux function
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    fst::scalarMultiplyDataData(KgradU,valZ_eval,gradU_eval);
    // Integrate stiffness term
    fst::integrate(res,KgradU,fe_vol_->gradNdetJ(),false);
    // Add force term
    fst::integrate(res,nforce_eval_,fe_vol_->NdetJ(),true);
    // Apply Dirichlet boundary conditions
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                res(cidx,fidx_[j][l]) = u_coeff(cidx,fidx_[j][l]) - (bdryCellDofValues_[i][j])(k,fidx_[j][l]);
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
    const Real one(1);
    // Get dimensions
    const int c = fe_vol_->gradN().extent_int(0);
    const int f = fe_vol_->gradN().extent_int(1);
    const int p = fe_vol_->gradN().extent_int(2);
    const int d = fe_vol_->gradN().extent_int(3);
    // Initialize Jacobian storage
    jac = scalar_view("jac", c, f, f);
    scalar_view valZ_eval("valZ_eval", c, p);
    scalar_view KgradN("KgradN", c, f, p, d);
    // Build density-dependent conductivity function
    fe_dens_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        valZ_eval(i,j) = minConductivity_
          + (one - minConductivity_) * std::pow(valZ_eval(i,j),SIMPpower_);
      }
    }
    // Build flux function
    fst::scalarMultiplyDataField(KgradN,valZ_eval,fe_vol_->gradN());
    // Integrate stiffness term
    fst::integrate(jac,KgradN,fe_vol_->gradNdetJ(),false);

    // Apply Dirichlet boundary conditions
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m = 0; m < f; ++m) {
                  jac(cidx,fidx_[j][l],m) = static_cast<Real>(0);
                }
                jac(cidx,fidx_[j][l],fidx_[j][l]) = static_cast<Real>(1);
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
    const Real one(1);
    // Get dimensions
    const int c  = fe_vol_->gradN().extent_int(0);
    const int f  = fe_vol_->gradN().extent_int(1);
    const int fd = fe_dens_->gradN().extent_int(1);
    const int p  = fe_vol_->gradN().extent_int(2);
    const int d  = fe_vol_->gradN().extent_int(3);
    // Initialize Jacobian storage
    jac = scalar_view("jac", c, f, fd);
    // Temporary storage
    scalar_view valZ_eval("valZ_eval", c, p);
    scalar_view gradU_eval("gradU_eval", c, p, d);
    scalar_view dKN("dKN", c, fd, p);
    scalar_view gradUgradN("gradUgradN", c, f, p);
    // Build density-dependent conductivity function
    fe_dens_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        valZ_eval(i,j) = SIMPpower_*(one - minConductivity_)
          * std::pow(valZ_eval(i,j),SIMPpower_-one);
      }
    }
    // Build derivative of conductivity function
    fst::scalarMultiplyDataField(dKN,valZ_eval,fe_dens_->N());
    // Integrate stiffness term
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    fst::dotMultiplyDataField(gradUgradN,gradU_eval,fe_vol_->gradNdetJ());
    fst::integrate(jac,gradUgradN,dKN,false);

    // Apply Dirichlet boundary conditions
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < f; ++m) {
                  jac(cidx,fidx_[j][l],m) = static_cast<Real>(0);
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
    throw Exception::Zero(">>> (PDE_Poisson_TopOpt::Hessian_11): Zero Hessian.");
  }

  void Hessian_12(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    const Real one(1);
    // Get dimensions
    const int c  = fe_vol_->gradN().extent_int(0);
    const int f  = fe_vol_->gradN().extent_int(1);
    const int fd = fe_dens_->gradN().extent_int(1);
    const int p  = fe_vol_->gradN().extent_int(2);
    const int d  = fe_vol_->gradN().extent_int(3);
    // Initialize Hessian storage
    hess = scalar_view("hess", c, fd, f);
    // Temporary storage
    scalar_view valZ_eval("valZ_eval", c, p);
    scalar_view l0_coeff("l0_coeff", c, f);
    scalar_view dKN("dKN", c, fd, p);
    scalar_view gradL_eval("gradL_eval", c, p, d);
    scalar_view gradLgradN("gradLgradN", c, f, p);
    // Build density-dependent conductivity function
    fe_dens_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        valZ_eval(i,j) = SIMPpower_*(one - minConductivity_)
          * std::pow(valZ_eval(i,j),SIMPpower_-one);
      }
    }
    // Apply Dirichlet conditions to the multipliers.
    Kokkos::deep_copy(l0_coeff,l_coeff);
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                l0_coeff(cidx,fidx_[j][l]) = static_cast<Real>(0);
              }
            }
          }
        }
      }
    }
    // Build derivative of conductivity function
    fst::scalarMultiplyDataField(dKN,valZ_eval,fe_dens_->N());
    // Integrate stiffness term
    fe_vol_->evaluateGradient(gradL_eval, l0_coeff);
    fst::dotMultiplyDataField(gradLgradN,gradL_eval,fe_vol_->gradNdetJ());
    fst::integrate(hess,dKN,gradLgradN,false);
  }

  void Hessian_21(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    const Real one(1);
    // Get dimensions
    const int c  = fe_vol_->gradN().extent_int(0);
    const int f  = fe_vol_->gradN().extent_int(1);
    const int fd = fe_dens_->gradN().extent_int(1);
    const int p  = fe_vol_->gradN().extent_int(2);
    const int d  = fe_vol_->gradN().extent_int(3);
    // Initialize Hessian storage
    hess = scalar_view("hess", c, f, fd);
    // Temporary storage
    scalar_view valZ_eval("valZ_eval", c, p);
    scalar_view l0_coeff("l0_coeff", c, f);
    scalar_view dKN("dKN", c, fd, p);
    scalar_view gradL_eval("gradL_eval", c, p, d);
    scalar_view gradLgradN("gradLgradN", c, f, p);
    // Build density-dependent conductivity function
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        valZ_eval(i,j) = SIMPpower_*(one - minConductivity_)
          * std::pow(valZ_eval(i,j),SIMPpower_-one);
      }
    }
    // Apply Dirichlet conditions to the multipliers.
    Kokkos::deep_copy(l0_coeff,l_coeff);
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                l0_coeff(cidx,fidx_[j][l]) = static_cast<Real>(0);
              }
            }
          }
        }
      }
    }
    // Build derivative of conductivity function
    fst::scalarMultiplyDataField(dKN,valZ_eval,fe_dens_->N());
    // Integrate stiffness term
    fe_vol_->evaluateGradient(gradL_eval, l0_coeff);
    fst::dotMultiplyDataField(gradLgradN,gradL_eval,fe_vol_->gradNdetJ());
    fst::integrate(hess,gradLgradN,dKN,false);
  }

  void Hessian_22(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    const Real one(1), two(2);
    // Get dimensions
    const int c  = fe_vol_->gradN().extent_int(0);
    const int f  = fe_vol_->gradN().extent_int(1);
    const int fd = fe_dens_->gradN().extent_int(1);
    const int p  = fe_vol_->gradN().extent_int(2);
    const int d  = fe_vol_->gradN().extent_int(3);
    // Initialize Hessian storage
    hess = scalar_view("hess", c, fd, fd);
    // Temporary storage
    scalar_view valZ_eval("valZ_eval", c, p);
    scalar_view l0_coeff("l0_coeff", c, f);
    scalar_view dKN("dKN", c, fd, p);
    scalar_view gradU_eval("gradU_eval", c, p, d);
    scalar_view gradL_eval("gradL_eval", c, p, d);
    scalar_view gradUgradL("gradUgradL", c, p);
    scalar_view NgradUgradL("NgradUgradL", c, fd, p);
    // Build density-dependent conductivity function
    fe_dens_->evaluateValue(valZ_eval, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        valZ_eval(i,j) = SIMPpower_ * (SIMPpower_ - one)*(one - minConductivity_)
          * std::pow(valZ_eval(i,j),SIMPpower_-two);
      }
    }
    // Apply Dirichlet conditions to the multipliers.
    Kokkos::deep_copy(l0_coeff,l_coeff);
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                l0_coeff(cidx,fidx_[j][l]) = static_cast<Real>(0);
              }
            }
          }
        }
      }
    }
    // Build derivative of conductivity function
    fst::scalarMultiplyDataField(dKN,valZ_eval,fe_dens_->N());
    // Integrate stiffness term
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    fe_vol_->evaluateGradient(gradL_eval, l0_coeff);
    fst::dotMultiplyDataData(gradUgradL,gradU_eval,gradL_eval);
    fst::scalarMultiplyDataField(NgradUgradL,gradUgradL,fe_dens_->NdetJ());
    fst::integrate(hess,dKN,NgradUgradL,false);
  }

  void RieszMap_1(scalar_view & riesz) override {
    // GET DIMENSIONS
    int c = fe_vol_->N().extent_int(0);
    int f = fe_vol_->N().extent_int(1);
    // INITIALIZE RIESZ
    riesz = scalar_view("riesz", c, f, f);
    Kokkos::deep_copy(riesz,fe_vol_->stiffMat());
    rst::add(riesz,fe_vol_->massMat());
  }

  void RieszMap_2(scalar_view & riesz) override {
    // GET DIMENSIONS
    int c = fe_dens_->N().extent_int(0);
    int f = fe_dens_->N().extent_int(1);
    // INITIALIZE RIESZ
    riesz = scalar_view("riesz", c, f, f);
    Kokkos::deep_copy(riesz,fe_dens_->massMat());
  }

  // This must be called before getFields2
  void setDensityFields(const std::vector<basis_ptr> &basisPtrs) {
    TEUCHOS_TEST_FOR_EXCEPTION(getFields2called_, std::invalid_argument,
      ">>> PDE-OPT/topo-opt/elasticity/src/pde_elasticity.hpp: Must call before getFields2!");

    basisPtrDens_  = basisPtrs[0];
    basisPtrsDens_ = basisPtrs;
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
    fe_vol_  = ROL::makePtr<fe_type>(volCellNodes_,basisPtr_,cellCub_);
    fe_dens_ = ROL::makePtr<fe_type>(volCellNodes_,basisPtrDens_,cellCub_,false);
    fidx_ = fe_vol_->getBoundaryDofs();
    computeForce();
    computeDirichlet();
  }

  void computeForce(void) {
    int c = fe_vol_->cubPts().extent_int(0);
    int p = fe_vol_->cubPts().extent_int(1);
    int d = fe_vol_->cubPts().extent_int(2);
    std::vector<Real> pt(d,0);
    force_eval_ = scalar_view("force_eval_",c,p);
    nforce_eval_ = scalar_view("nforce_eval_",c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (fe_vol_->cubPts())(i,j,k);
        }
        force_eval_(i,j)  = ForceFunc(pt);
        nforce_eval_(i,j) = -force_eval_(i,j);
      }
    }
  }

  void computeDirichlet(void) {
    // Compute Dirichlet values at DOFs.
    int d = basisPtr_->getBaseCellTopology().getDimension();
    int numSidesets = bdryCellLocIds_.size();
    bdryCellDofValues_.resize(numSidesets);
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellDofValues_[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
        int f = basisPtr_->getCardinality();
	std::stringstream name;
	name << "bdryCellDofValues" << i << j;
        bdryCellDofValues_[i][j] = scalar_view(name.str(), c, f);
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
            (bdryCellDofValues_[i][j])(k, l) = dirichletFunc(dofpoint, i, j);
          }
        }
      }
    }
  }

  const ROL::Ptr<fe_type> getFE(void) const {
    return fe_vol_;
  }

  const ROL::Ptr<fe_type> getDensityFE(void) const {
    return fe_dens_;
  }

  const scalar_view getForce(void) const {
    return force_eval_;
  }

}; // PDE_Poisson

#endif
