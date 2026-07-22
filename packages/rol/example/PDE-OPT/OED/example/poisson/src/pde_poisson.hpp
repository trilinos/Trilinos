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

#ifndef PDE_OED_POISSON_HPP
#define PDE_OED_POISSON_HPP

#include "../../../../TOOLS/pde.hpp"
#include "../../../../TOOLS/fe.hpp"
#include "../../../../TOOLS/fieldhelper.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"


template <class Real>
class PDE_OED_Poisson : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>> basisPtr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real>> cellCub_;
  ROL::Ptr<Intrepid::Cubature<Real>> bdryCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real>> volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real>> fe_;
  std::vector<std::vector<ROL::Ptr<FE<Real>>>> feBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellDofValues_;

  Real diffusivity(const std::vector<Real> &x) const {
    const Real half(0.5), rad(0.2), k0(10), k1(1e-1);
    Real norm(0);
    for (const auto xi : x) norm += (xi-half)*(xi-half);
    norm = std::sqrt(norm);
    return (norm <= rad ? k1 : k0);
  }

  void computeDiffusivity(const ROL::Ptr<Intrepid::FieldContainer<Real>> &kappa) const {
    const int c = fe_->gradN()->dimension(0);
    const int p = fe_->gradN()->dimension(2);
    const int d = fe_->gradN()->dimension(3);
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) x[k] = (*fe_->cubPts())(i,j,k);
        (*kappa)(i,j) = diffusivity(x);
      }
    }
  }

  void computeNeumannControl(const ROL::Ptr<Intrepid::FieldContainer<Real>> &Bz,
                             const ROL::Ptr<const std::vector<Real>> &zp,
                             const int sideset,
                             const int deriv = 0) const {
    const int c = Bz->dimension(0);
    const int p = Bz->dimension(1);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*Bz)(i,j) = (deriv==0 ? (*zp)[sideset]
                   : (deriv==1 ? static_cast<Real>(1)
                   : static_cast<Real>(0)));
      }
    }
  }

  ROL::Ptr<Intrepid::FieldContainer<Real>> getBoundaryCoeff(
      const Intrepid::FieldContainer<Real> & cell_coeff,
      int sideSet, int cell) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][cell];
    const int numCellsSide = bdryCellLocId.size();
    const int f = basisPtr_->getCardinality();
    
    ROL::Ptr<Intrepid::FieldContainer<Real >> bdry_coeff = 
      ROL::makePtr<Intrepid::FieldContainer<Real >>(numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        (*bdry_coeff)(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

public:
  PDE_OED_Poisson(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    std::string elemtype = parlist.sublist("Geometry").get("Element Shape","quad");
    if (elemtype == "quad") {
      basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
    }
    else if (elemtype == "tri") {
      basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_TRI_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
    }
    else {
      throw ROL::Exception::NotImplemented(">>> Element type not implemented.");
    }
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);  // acoustic pressure
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();            // get the cell type from the basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                           // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 2);        // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                           // create default cubature
    int d = cellType.getDimension();

    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2); // set cubature degree, e.g., 2
    bdryCub_ = cubFactory.create(bdryCellType, bdryCubDegree);
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real>> & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    const int d = fe_->gradN()->dimension(3);
    // Initialize residuals.
    res = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f);
    // Initialize storage
    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval, gradU_eval, kappa, kappaU;
    valU_eval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    gradU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    kappa      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    kappaU     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    // Evaluate/interpolate finite element fields on cells.
    fe_->evaluateValue(valU_eval, u_coeff);
    fe_->evaluateGradient(gradU_eval, u_coeff);
    // Build wave number
    computeDiffusivity(kappa);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*kappaU,*kappa,*gradU_eval);

    /*******************************************************************/
    /*** Evaluate weak form of the residual.****************************/
    /*******************************************************************/
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *kappaU,           // kappa grad U
                                                  *fe_->gradNdetJ(), // grad N
                                                  Intrepid::COMP_CPP,
                                                  false);

    // APPLY BOUNDARY CONDITIONS
    const int numSidesets = bdryCellLocIds_.size();
    for (int s = 0; s < numSidesets; ++s) {
      const int numLocalSideIds = bdryCellLocIds_[s].size();
      const int numCubPerSide = bdryCub_->getNumPoints();
      for (int j = 0; j < numLocalSideIds; ++j) {
        const int numCellsSide = bdryCellLocIds_[s][j].size();
        if (numCellsSide) {
          ROL::Ptr<Intrepid::FieldContainer<Real>> BCres, BCcomp;
          BCres  = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f);
          BCcomp = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          // Compute control operator
          computeNeumannControl(BCcomp,z_param,s,0);
          // Integrate residual
          Intrepid::FunctionSpaceTools::integrate<Real>(*BCres,
                                                        *BCcomp,
                                                        *(feBdry_[s][j]->NdetJ()),
                                                        Intrepid::COMP_CPP, false);
          // Add Robin and Neumann control residual to volume residual
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[s][j][k];
            for (int l = 0; l < f; ++l) {
              (*res)(cidx,l) -= (*BCres)(k,l);
            }
          }
        }
      }
    }
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
    // Initialize Jacobians.
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // Initialize storage.
    ROL::Ptr<Intrepid::FieldContainer<Real>> kappa, kappaN;
    kappa  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    kappaN = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p, d);
    // Build wave number
    computeDiffusivity(kappa);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*kappaN,*kappa,*fe_->gradN());

    /*** Evaluate weak form of the Jacobian. ***/
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *kappaN,
                                                  *fe_->gradNdetJ(),
                                                  Intrepid::COMP_CPP,
                                                  false);
  }


  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Helmholtz_OCT::Jacobian_2): Jacobian is zero.");
  }

  void Jacobian_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & jac,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    // ADD CONTROL TERM TO RESIDUAL
    const int numSidesets = bdryCellLocIds_.size();
    for (int s = 0; s < numSidesets; ++s) {
      jac[s] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      const int numLocalSideIds = bdryCellLocIds_[s].size();
      const int numCubPerSide = bdryCub_->getNumPoints();
      for (int j = 0; j < numLocalSideIds; ++j) {
        const int numCellsSide = bdryCellLocIds_[s][j].size();
        if (numCellsSide) {
          ROL::Ptr<Intrepid::FieldContainer<Real>> neumJac, Bz;
          // Compute control operator
          Bz = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          computeNeumannControl(Bz,z_param,s,1);
          // Compute Neumann residual
          neumJac = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f);
          Intrepid::FunctionSpaceTools::integrate<Real>(*neumJac,
                                                        *Bz,
                                                        *(feBdry_[s][j]->NdetJ()),
                                                        Intrepid::COMP_CPP, false);
          // Add Robin and Neumann control residual to volume residual
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[s][j][k];
            for (int l = 0; l < f; ++l) { 
              (*jac[s])(cidx,l) -= (*neumJac)(k,l);
            }
          }
        }
      }
    }
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Helmholtz_OCT::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Helmholtz_OCT::Hessian_12): Hessian is zero.");
  }

  void Hessian_13(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Helmholtz_OCT::Hessian_13): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Helmholtz_OCT::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Helmholtz_OCT::Hessian_22): Hessian is zero.");
  }

  void Hessian_23(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Helmholtz_OCT::Hessian_23): Hessian is zero.");
  }

  void Hessian_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Helmholtz_OCT::Hessian_31): Hessian is zero.");
  }

  void Hessian_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Helmholtz_OCT::Hessian_32): Hessian is zero.");
  }

  void Hessian_33(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > & hess,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                          const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                          const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Helmholtz_OCT::Hessian_33): Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    //throw Exception::NotImplemented(">>> (PDE_TopoOpt::RieszMap_1): Not implemented.");
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    // Initialize Riesz Map.
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f);
    *riesz = *(fe_->stiffMat());
    Intrepid::RealSpaceTools<Real>::add(*riesz,*(fe_->massMat()));
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    //throw Exception::NotImplemented(">>> (PDE_TopoOpt::RieszMap_2): Not implemented.");
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    // Initialize Jacobians.
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f);
    *riesz = *(fe_->massMat());
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields() {
    return basisPtrs_;
  }

  void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real>> &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds) {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fe_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtr_,cellCub_);
    fidx_ = fe_->getBoundaryDofs();
    // Construct boundary FE
    const int numSidesets = bdryCellNodes.size();
    feBdry_.resize(numSidesets);
    for (int s = 0; s < numSidesets; ++s) {
      const int numLocSides = bdryCellNodes[s].size();
      feBdry_[s].resize(numLocSides);
      for (int j = 0; j < numLocSides; ++j) {
        if (bdryCellNodes[s][j] != ROL::nullPtr) {
          feBdry_[s][j] = ROL::makePtr<FE<Real>>(bdryCellNodes[s][j],basisPtr_,bdryCub_,j);
        }
      }
    }
  }

  const ROL::Ptr<FE<Real>> getFE(void) const {
    return fe_;
  }

  const std::vector<std::vector<ROL::Ptr<FE<Real>>>> getBdryFE(void) const {
    return feBdry_;
  }

  const std::vector<std::vector<int>> getBdryCellLocIds(const int sideset = 0) const {
    return bdryCellLocIds_[sideset];
  }

}; // PDE_Helmholtz_OCT


#endif
