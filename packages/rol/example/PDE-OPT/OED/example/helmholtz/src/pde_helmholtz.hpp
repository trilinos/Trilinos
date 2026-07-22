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

#ifndef PDE_HELMHOLTZ_OCT_REAL_HPP
#define PDE_HELMHOLTZ_OCT_REAL_HPP

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
class PDE_Helmholtz_OCT : public PDE<Real> {
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

  Real waveNumber_;
  Real RoiWaveNumber_;
  Real RoiIncWaveNumber_;
  Real RoiRadius_;
  Real impFactor_;
  Real freq_;
  Real density_;

  Real refractiveIndex(const std::vector<Real> &x) const {
    Real xnorm(0), xdnorm(0), val(0);
    const int d = x.size();
    const Real sqr2(std::sqrt(2.0)), half(0.5);
    const std::vector<Real> xd = {-half/sqr2*RoiRadius_, -half/sqr2*RoiRadius_};
    for (int i = 0; i < d; ++i) {
      xnorm += x[i]*x[i];
      xdnorm += (x[i]-xd[i])*(x[i]-xd[i]);
    }
    xnorm  = std::sqrt(xnorm);
    xdnorm = std::sqrt(xdnorm);
    val   = ((xnorm <= RoiRadius_)
          ? ((xdnorm <= half*RoiRadius_)
          ? RoiIncWaveNumber_ : RoiWaveNumber_) : waveNumber_);
    return val;
  }

  void computeRefractiveIndex(const ROL::Ptr<Intrepid::FieldContainer<Real>> &kappa) const {
    const int c = fe_->gradN()->dimension(0);
    const int p = fe_->gradN()->dimension(2);
    const int d = fe_->gradN()->dimension(3);
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (*fe_->cubPts())(i,j,k);
        }
        (*kappa)(i,j) = std::pow(refractiveIndex(x), 2);
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
  PDE_Helmholtz_OCT(Teuchos::ParameterList &parlist) {
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

    const Real pi2    = static_cast<Real>(2.0*M_PI);
    // Speed of sound of air
    Real speed        = parlist.sublist("Problem").get("Speed of Sound",               343.0);
    // Speed of sound of aluminum
    Real RoiSpeed     = parlist.sublist("Problem").get("ROI Speed of Sound",           5100.0);
    // Speed of sound of rubber
    Real RoiIncSpeed  = parlist.sublist("Problem").get("ROI Inclusion Speed of Sound", 1600.0);
    // Between middle C and the middle C on the treble clef
    freq_             = parlist.sublist("Problem").get("Frequency",                    400.0);
    // Density of air at 20C and standard atmospheric pressure
    density_          = parlist.sublist("Problem").get("Density",                      1.204);
    // Impedance for air v.s. wood (420.175 PA s/m / 0.5e6 PA s/m)
    impFactor_        = parlist.sublist("Problem").get("Impedance Factor",             8.4e-4);
    RoiRadius_        = parlist.sublist("Problem").get("ROI Radius",                   0.5);
    waveNumber_       = pi2*freq_/speed;
    RoiWaveNumber_    = pi2*freq_/RoiSpeed;
    RoiIncWaveNumber_ = pi2*freq_/RoiIncSpeed;
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
    kappaU     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    // Evaluate/interpolate finite element fields on cells.
    fe_->evaluateValue(valU_eval, u_coeff);
    fe_->evaluateGradient(gradU_eval, u_coeff);
    // Build wave number
    computeRefractiveIndex(kappa);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*kappaU,*kappa,*valU_eval);
    Intrepid::RealSpaceTools<Real>::scale(*kappaU,static_cast<Real>(-1));

    /*******************************************************************/
    /*** Evaluate weak form of the residual.****************************/
    /*******************************************************************/
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *gradU_eval,       // grad U
                                                  *fe_->gradNdetJ(), // grad N
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *kappaU,           // -kappa2 U
                                                  *fe_->NdetJ(),     // N
                                                  Intrepid::COMP_CPP,
                                                  true);

    // APPLY BOUNDARY CONDITIONS
    const int numSidesets = bdryCellLocIds_.size();
    for (int s = 0; s < numSidesets; ++s) {
      if (s != 8) { // s = 8 -> Speaker cabinents (0 impedance)
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
                (*res)(cidx,l) -= density_ * freq_ * (*BCres)(k,l);
              }
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
    // Initialize Jacobians.
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // Initialize storage.
    ROL::Ptr<Intrepid::FieldContainer<Real>> kappa, kappaN;
    kappa  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    kappaN = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    // Build wave number
    computeRefractiveIndex(kappa);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*kappaN,*kappa,*fe_->N());
    Intrepid::RealSpaceTools<Real>::scale(*kappaN,static_cast<Real>(-1));

    /*** Evaluate weak form of the Jacobian. ***/
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *kappaN,
                                                  *fe_->NdetJ(),
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::RealSpaceTools<Real>::add(*jac,*fe_->stiffMat());
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
      if (s != 8) { // s = 8 -> Speaker cabinent (zero impedance)
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
                (*jac[s])(cidx,l) -= density_*freq_*(*neumJac)(k,l);
              }
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
