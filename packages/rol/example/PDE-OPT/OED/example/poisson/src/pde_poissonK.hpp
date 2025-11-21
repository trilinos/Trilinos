// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde_poisson.hpp
    \brief Implements the local PDE interface for the optimal control of
           Poisson.
*/

#ifndef PDE_OED_POISSONK_HPP
#define PDE_OED_POISSONK_HPP

#include "../../../../TOOLS/pdeK.hpp"
#include "../../../../TOOLS/feK.hpp"
#include "../../../../TOOLS/fieldhelperK.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"

#include "ROL_Ptr.hpp"


template <class Real, class DeviceType>
class PDE_OED_Poisson : public PDE<Real,DeviceType> {
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
  std::vector<std::vector<ROL::Ptr<fe_type>>> feBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<scalar_view>> bdryCellDofValues_;

  Real diffusivity(const std::vector<Real> &x) const {
    const Real half(0.5), rad(0.2), k0(10), k1(1e-1);
    Real norm(0);
    for (const auto xi : x) norm += (xi-half)*(xi-half);
    norm = std::sqrt(norm);
    return (norm <= rad ? k1 : k0);
  }

  void computeDiffusivity(scalar_view &kappa) const {
    const int c = fe_->gradN().extent_int(0);
    const int p = fe_->gradN().extent_int(2);
    const int d = fe_->gradN().extent_int(3);
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) x[k] = (fe_->cubPts())(i,j,k);
        kappa(i,j) = diffusivity(x);
      }
    }
  }

  void computeNeumannControl(scalar_view &Bz,
                             const ROL::Ptr<const std::vector<Real>> &zp,
                             const int sideset,
                             const int deriv = 0) const {
    const int c = Bz.extent_int(0);
    const int p = Bz.extent_int(1);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        Bz(i,j) = (deriv==0 ? (*zp)[sideset]
                : (deriv==1 ? static_cast<Real>(1)
                : static_cast<Real>(0)));
      }
    }
  }

  scalar_view getBoundaryCoeff(const scalar_view cell_coeff, int sideSet, int cell) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][cell];
    const int numCellsSide = bdryCellLocId.size();
    const int f = basisPtr_->getCardinality();
    
    scalar_view bdry_coeff("bdry_coeff", numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j)
        bdry_coeff(i, j) = cell_coeff(bdryCellLocId[i], j);
    }
    return bdry_coeff;
  }

public:
  PDE_OED_Poisson(ROL::ParameterList &parlist) {
    // Finite element fields.
    std::string elemtype = parlist.sublist("Geometry").get("Element Shape","quad");
    if (elemtype == "quad")
      basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType,Real,Real>>();
    else if (elemtype == "tri")
      basisPtr_ = ROL::makePtr<Intrepid2::Basis_HGRAD_TRI_C1_FEM<DeviceType,Real,Real>>();
    else
      throw ROL::Exception::NotImplemented(">>> Element type not implemented.");
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);  // acoustic pressure
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();        // get the cell type from the basis
    Intrepid2::DefaultCubatureFactory cubFactory;                            // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 2);    // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create<DeviceType,Real,Real>(cellType, cubDegree); // create default cubature
    int d = cellType.getDimension();

    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2); // set cubature degree, e.g., 2
    bdryCub_ = cubFactory.create<DeviceType,Real,Real>(bdryCellType, bdryCubDegree);
  }

  void residual(scalar_view & res,
                const scalar_view u_coeff,
                const scalar_view z_coeff = scalar_view(),
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);
    const int p = fe_->gradN().extent_int(2);
    const int d = fe_->gradN().extent_int(3);
    // Initialize residuals.
    res = scalar_view("res",c,f);
    // Initialize storage
    scalar_view valU_eval("valU_eval", c, p);
    scalar_view gradU_eval("gradU_eval", c, p, d);
    scalar_view kappa("kappa", c, p);
    scalar_view kappaU("kappaU", c, p, d);
    // Evaluate/interpolate finite element fields on cells.
    fe_->evaluateValue(valU_eval, u_coeff);
    fe_->evaluateGradient(gradU_eval, u_coeff);
    // Build wave number
    computeDiffusivity(kappa);
    fst::scalarMultiplyDataData(kappaU,kappa,gradU_eval);

    /*******************************************************************/
    /*** Evaluate weak form of the residual.****************************/
    /*******************************************************************/
    fst::integrate(res,kappaU,fe_->gradNdetJ(),false);

    // APPLY BOUNDARY CONDITIONS
    const int numSidesets = bdryCellLocIds_.size();
    for (int s = 0; s < numSidesets; ++s) {
      const int numLocalSideIds = bdryCellLocIds_[s].size();
      const int numCubPerSide = bdryCub_->getNumPoints();
      for (int j = 0; j < numLocalSideIds; ++j) {
        const int numCellsSide = bdryCellLocIds_[s][j].size();
        if (numCellsSide) {
          scalar_view BCres("BCres", numCellsSide, f);
          scalar_view BCcomp("BCcomp", numCellsSide, numCubPerSide);
          // Compute control operator
          computeNeumannControl(BCcomp,z_param,s,0);
          // Integrate residual
          fst::integrate(BCres,BCcomp,feBdry_[s][j]->NdetJ(),false);
          // Add Robin and Neumann control residual to volume residual
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[s][j][k];
            for (int l = 0; l < f; ++l)
              res(cidx,l) -= BCres(k,l);
          }
        }
      }
    }
  }

  void Jacobian_1(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);
    // Initialize Jacobians.
    jac = scalar_view("jac", c, f, f);
    // Initialize storage.
    scalar_view kappa("kappa", c, p);
    scalar_view kappaN("kappaN", c, f, p, d);
    // Build wave number
    computeDiffusivity(kappa);
    fst::scalarMultiplyDataField(kappaN,kappa,fe_->gradN());

    /*** Evaluate weak form of the Jacobian. ***/
    fst::integrate(jac,kappaN,fe_->gradNdetJ(),false);
  }


  void Jacobian_2(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Helmholtz_OCT::Jacobian_2): Jacobian is zero.");
  }

  void Jacobian_3(std::vector<scalar_view> & jac,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // GET DIMENSIONS
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);
    // ADD CONTROL TERM TO RESIDUAL
    const int numSidesets = bdryCellLocIds_.size();
    for (int s = 0; s < numSidesets; ++s) {
      jac[s] = scalar_view("jac", c, f);
      const int numLocalSideIds = bdryCellLocIds_[s].size();
      const int numCubPerSide = bdryCub_->getNumPoints();
      for (int j = 0; j < numLocalSideIds; ++j) {
        const int numCellsSide = bdryCellLocIds_[s][j].size();
        if (numCellsSide) {
          // Compute control operator
          scalar_view Bz("Bz", numCellsSide, numCubPerSide);
          computeNeumannControl(Bz,z_param,s,1);
          // Compute Neumann residual
          scalar_view neumJac("neumJac", numCellsSide, f);
          fst::integrate(neumJac,Bz,feBdry_[s][j]->NdetJ(),false);
          // Add Robin and Neumann control residual to volume residual
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[s][j][k];
            for (int l = 0; l < f; ++l) { 
              (jac[s])(cidx,l) -= neumJac(k,l);
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
    throw Exception::Zero(">>> (PDE_Helmholtz_OCT::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Helmholtz_OCT::Hessian_12): Hessian is zero.");
  }

  void Hessian_13(std::vector<scalar_view> & hess,
                          const scalar_view l_coeff,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Helmholtz_OCT::Hessian_13): Hessian is zero.");
  }

  void Hessian_21(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Helmholtz_OCT::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Helmholtz_OCT::Hessian_22): Hessian is zero.");
  }

  void Hessian_23(std::vector<scalar_view> & hess,
                          const scalar_view l_coeff,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Helmholtz_OCT::Hessian_23): Hessian is zero.");
  }

  void Hessian_31(std::vector<scalar_view> & hess,
                          const scalar_view l_coeff,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Helmholtz_OCT::Hessian_31): Hessian is zero.");
  }

  void Hessian_32(std::vector<scalar_view> & hess,
                          const scalar_view l_coeff,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Helmholtz_OCT::Hessian_32): Hessian is zero.");
  }

  void Hessian_33(std::vector<std::vector<scalar_view>> & hess,
                          const scalar_view l_coeff,
                          const scalar_view u_coeff,
                          const scalar_view z_coeff = scalar_view(),
                          const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Helmholtz_OCT::Hessian_33): Hessian is zero.");
  }

  void RieszMap_1(scalar_view & riesz) override {
    //throw Exception::NotImplemented(">>> (PDE_TopoOpt::RieszMap_1): Not implemented.");
    // Retrieve dimensions.
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    // Initialize Riesz Map.
    riesz = scalar_view("riesz1",c,f,f);
    Kokkos::deep_copy(riesz,fe_->stiffMat());
    rst::add(riesz,fe_->massMat());
  }

  void RieszMap_2(scalar_view & riesz) override {
    //throw Exception::NotImplemented(">>> (PDE_TopoOpt::RieszMap_2): Not implemented.");
    // Retrieve dimensions.
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    // Initialize Jacobians.
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
    fidx_ = fe_->getBoundaryDofs();
    // Construct boundary FE
    const int numSidesets = bdryCellNodes.size();
    feBdry_.resize(numSidesets);
    for (int s = 0; s < numSidesets; ++s) {
      const int numLocSides = bdryCellNodes[s].size();
      feBdry_[s].resize(numLocSides);
      for (int j = 0; j < numLocSides; ++j) {
        if (bdryCellNodes[s][j] != scalar_view()) {
          feBdry_[s][j] = ROL::makePtr<fe_type>(bdryCellNodes[s][j],basisPtr_,bdryCub_,j);
        }
      }
    }
  }

  const ROL::Ptr<fe_type> getFE(void) const {
    return fe_;
  }

  const std::vector<std::vector<ROL::Ptr<fe_type>>> getBdryFE(void) const {
    return feBdry_;
  }

  const std::vector<std::vector<int>> getBdryCellLocIds(const int sideset = 0) const {
    return bdryCellLocIds_[sideset];
  }

}; // PDE_Helmholtz_OCT


#endif
