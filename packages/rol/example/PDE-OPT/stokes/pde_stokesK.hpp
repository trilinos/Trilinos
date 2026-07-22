// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde_stokes.hpp
    \brief Implements the local PDE interface for the Stokes control problem.
*/

#ifndef PDE_STOKESK_HPP
#define PDE_STOKESK_HPP

#include "../TOOLS/pdeK.hpp"
#include "../TOOLS/feK.hpp"
#include "../TOOLS/fieldhelperK.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real, class DeviceType>
class PDE_Stokes : public PDE<Real, DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;
private:
  // Finite element basis information
  basis_ptr basisPtrVel_, basisPtrPrs_;
  std::vector<basis_ptr> basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid2::Cubature<DeviceType,Real,Real>> cellCub_, bdryCub_;
  // Cell node information
  scalar_view volCellNodes_;
  std::vector<std::vector<scalar_view>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<fe_type> feVel_, fePrs_;
  std::vector<std::vector<ROL::Ptr<fe_type>>> feVelBdry_, fePrsBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fvidx_, fpidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<scalar_view>> bdryCellVDofValues_;
  // Field pattern, offsets, etc.
  std::vector<std::vector<int>> fieldPattern_; // local Field/DOF pattern; set from DOF manager 
  int numFields_;                              // number of fields (equations in the PDE)
  int numDofs_;                                // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                    // for each field, a counting offset
  std::vector<int> numFieldDofs_;              // for each field, number of degrees of freedom
  
  // Problem parameters.
  Real Re_;
  bool pinPressure_, dirType_;

  ROL::Ptr<FieldHelper<Real,DeviceType>> fieldHelper_;

  Real velocityDirichletFunc(const std::vector<Real> & coords, int sideset, int locSideId, int dir) const {
    Real val(0);
    if (dirType_ == 0) {
      if ((sideset==3) && (dir==0)) {
        val = static_cast<Real>(1);
      }
    }
    else if (dirType_ == 1) {
      Real one(1);
      if ((sideset==1) && (dir==0)) {
        val = coords[1]*(one-coords[1]);
      }
      else if ((sideset==2) && (dir==0)) {
        val = coords[1]*(one-coords[1]);
      }
    }
    return val;
  }

  void computeDirichlet(void) {
    // Compute Dirichlet values at DOFs.
    int d  = basisPtrVel_->getBaseCellTopology().getDimension();
    int fv = basisPtrVel_->getCardinality();
    int numSidesets = bdryCellLocIds_.size();
    bdryCellVDofValues_.resize(numSidesets);
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellVDofValues_[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
	std::stringstream name;
	name << "bdryCellDofValues" << i << j;
        bdryCellVDofValues_[i][j] = scalar_view(name.str(), c, fv, d);
        scalar_view Vcoords;
        Vcoords = scalar_view("Vcoords", c, fv, d);
        if (c > 0) {
          feVel_->computeDofCoords(Vcoords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<fv; ++l) {
            std::vector<Real> dofpoint(d);
            //std::cout << "Sideset " << i << " LocalSide " << j << "  Cell " << k << "  Field " << l << "  Coord ";
            for (int m=0; m<d; ++m) {
              dofpoint[m] = Vcoords(k, l, m);
              //std::cout << dofpoint[m] << "  ";
            }

            for (int m=0; m<d; ++m) {
              (bdryCellVDofValues_[i][j])(k, l, m) = velocityDirichletFunc(dofpoint, i, j, m);
              //std::cout << "  " << m << "-Value " << DirichletFunc(dofpoint, i, j, m);
            }
            //std::cout << std::endl;
          }
        }
      }
    }
  }

public:
  PDE_Stokes(ROL::ParameterList &parlist) {
    // Finite element fields -- NOT DIMENSION INDEPENDENT!
    basisPtrVel_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<DeviceType,Real,Real>>();
    basisPtrPrs_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType,Real,Real>>();
    // Volume quadrature rules.
    shards::CellTopology cellType = basisPtrVel_->getBaseCellTopology();     // get the cell type from any basis
    Intrepid2::DefaultCubatureFactory cubFactory;                            // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 4);    // set cubature degree, e.g., 4
    cellCub_ = cubFactory.create<DeviceType,Real,Real>(cellType, cubDegree); // create default cubature
    // Boundary quadrature rules.
    int d = cellType.getDimension();
    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",4); // set cubature degree, e.g., 4
    bdryCub_ = cubFactory.create<DeviceType,Real,Real>(bdryCellType, bdryCubDegree);
    // Fill finite element basis container
    basisPtrs_.clear();
    basisPtrs_.resize(d,basisPtrVel_);  // Velocity
    basisPtrs_.push_back(basisPtrPrs_); // Pressure

    // Reynold's Number
    Re_ = parlist.sublist("Problem").get("Reynolds Number",100.0);

    // Pin pressure
    pinPressure_ = parlist.sublist("Problem").get("Pin Pressure",true);

    // Dirichlet Type
    dirType_ = parlist.sublist("Problem").get("Dirichlet Type",0);

    numDofs_ = 0;
    numFields_ = basisPtrs_.size();
    offset_.resize(numFields_);
    numFieldDofs_.resize(numFields_);
    for (int i=0; i<numFields_; ++i) {
      if (i==0)
        offset_[i] = 0;
      else
        offset_[i] = offset_[i-1] + basisPtrs_[i-1]->getCardinality();
      numFieldDofs_[i] = basisPtrs_[i]->getCardinality();
      numDofs_ += numFieldDofs_[i];
    }
  }

  void residual(scalar_view & res,
                const scalar_view u_coeff,
                const scalar_view z_coeff = scalar_view(),
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    int c  = u_coeff.extent_int(0);
    int p  = cellCub_->getNumPoints();
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize residuals.
    std::vector<scalar_view> R(d+1);
    for (int i = 0; i < d; ++i) {
      std::stringstream name;
      name << "res_vel" << i;
      R[i] = scalar_view(name.str(), c, fv);
    }
    R[d] = scalar_view("res_pres", c, fp);

    // Split u_coeff into components.
    std::vector<scalar_view> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);

    // Evaluate/interpolate finite element fields on cells.
    scalar_view valPres_eval("valPres_eval", c, p);
    fePrs_->evaluateValue(valPres_eval, U[d]);
    // Evaluate/interpolate gradient of finite element fields on cells.
    std::vector<scalar_view> gradVel_vec(d);
    for (int i = 0; i < d; ++i) {
      std::stringstream name;
      name << "gradVel_vec" << i;
      gradVel_vec[i] = scalar_view(name.str(), c, p, d);
      feVel_->evaluateGradient(gradVel_vec[i], U[i]);
    }

    // Assemble the velocity vector and its divergence.
    scalar_view divVel_eval("divVel_eval", c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        divVel_eval(i,j) = static_cast<Real>(0);
        for (int k = 0; k < d; ++k) {
          divVel_eval(i,j) += (gradVel_vec[k])(i,j,k);
        }
      }
    }

    // Scale pressure and velocity gradient
    const Real zero(0), one(1);
    for (int i = 0; i < d; ++i) {
      // Multiply velocity gradients with viscosity.
      rst::scale(gradVel_vec[i], one/Re_);
    }
    // Negative pressure
    rst::scale(valPres_eval, -one);

    /**************************************************************************/
    /*** EVALUATE WEAK FORM OF RESIDUAL ***************************************/
    /**************************************************************************/
    // Velocity equation.
    for (int i = 0; i < d; ++i) {
      fst::integrate(R[i], gradVel_vec[i], feVel_->gradNdetJ(), false);
      fst::integrate(R[i], valPres_eval, feVel_->DNDdetJ(i), true);
    }
    // Pressure equation.
    fst::integrate(R[d], divVel_eval, fePrs_->NdetJ(), false);
    rst::scale(R[d], -one);

    /**************************************************************************/
    /*** APPLY BOUNDARY CONDITIONS ********************************************/
    /**************************************************************************/
    // --> No slip boundaries: i=0,1,3
    // -->     Lid boundaries: i=2
    // -->       Pressure pin: i=4
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        // Velocity Boundary Conditions
        if (i!=4) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numVBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numVBdryDofs; ++l) {
                for (int m = 0; m < d; ++m) {
                  R[m](cidx,fvidx_[j][l]) = U[m](cidx,fvidx_[j][l]) - (bdryCellVDofValues_[i][j])(k,fvidx_[j][l],m);
                }
              }
            }
          }
        }
        // Pressure pinning
        if (i==4 && pinPressure_) {
          Real val = (dirType_ == 1) ? one/Re_ : zero;
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int l = 0, cidx = bdryCellLocIds_[i][j][k];
              R[d](cidx,fpidx_[j][l]) = U[d](cidx,fpidx_[j][l]) - val;
            }
          }
        }
      }
    }

    // Combine the residuals.
    fieldHelper_->combineFieldCoeff(res, R);
  }

  void Jacobian_1(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    int c  = u_coeff.extent_int(0);
    int p  = cellCub_->getNumPoints();
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<scalar_view>> J(d+1);
    for (int i = 0; i < d+1; ++i) J[i].resize(d+1);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        std::stringstream name;
        name << "jac" << i << j;
        J[i][j] = scalar_view(name.str(), c, fv, fv);
      }
      std::stringstream name1, name2;
      name1 << "jac" << d << i;
      name2 << "jac" << i << d;
      J[d][i] = scalar_view(name1.str(), c, fp, fv);
      J[i][d] = scalar_view(name2.str(), c, fv, fp);
    }
    std::stringstream name;
    name << "jac" << d << d;
    J[d][d] = scalar_view(name.str(), c, fp, fp);

    // Multiply velocity gradients with viscosity.
    const Real zero(0), one(1);
    scalar_view nuGradPhi_eval("nuGradPhi_eval", c, fv, p, d);
    Kokkos::deep_copy(nuGradPhi_eval,feVel_->gradN());
    rst::scale(nuGradPhi_eval, one/Re_);
    // Negative pressure basis.
    scalar_view negPrsPhi("negPrsPhi", c, fp, p);
    Kokkos::deep_copy(negPrsPhi,fePrs_->N());
    rst::scale(negPrsPhi, -one);

    /*** Evaluate weak form of the Jacobian. ***/
    for (int i = 0; i < d; ++i) {
      // Velocity components
      fst::integrate(J[i][i], nuGradPhi_eval, feVel_->gradNdetJ(), true);
      // Pressure components
      fst::integrate(J[i][d], feVel_->DNDdetJ(i), negPrsPhi, false);
      fst::integrate(J[d][i], fePrs_->NdetJ(), feVel_->DND(i), false);
      rst::scale(J[d][i], -one);
    }

    // APPLY BOUNDARY CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      // DIRICHLET CONDITIONS
      for (int i = 0; i < numSideSets; ++i) {
        // Velocity Boundary Conditions
        if (i!=4) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numVBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numVBdryDofs; ++l) {
                for (int m = 0; m < fv; ++m) {
                  for (int n = 0; n < d; ++n) {
                    for (int q = 0; q < d; ++q) {
                      (J[n][q])(cidx,fvidx_[j][l],m) = zero;
                    }
                    (J[n][n])(cidx,fvidx_[j][l],fvidx_[j][l]) = one;
                  }
                }
                for (int m = 0; m < d; ++m) {
                  for (int n = 0; n < fp; ++n) {
                    (J[m][d])(cidx,fvidx_[j][l],n) = zero;
                  }
                }
              }
            }
          }
        }
        // Pressure pinning
        if (i==4 && pinPressure_) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int l = 0, cidx = bdryCellLocIds_[i][j][k];
              for (int m = 0; m < fv; ++m) {
                for (int n = 0; n < d; ++n) {
                  (J[d][n])(cidx,fpidx_[j][l],m) = zero;
                }
              }
              for (int m = 0; m < fp; ++m) {
                (J[d][d])(cidx,fpidx_[j][l],m) = zero;
              }
              (J[d][d])(cidx,fpidx_[j][l],fpidx_[j][l]) = one;
            }
          }
        }
      }
    }
    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }


  void Jacobian_2(scalar_view & jac,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Retrieve dimensions.
    int c  = u_coeff.extent_int(0);
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<scalar_view>> J(d+1);
    for (int i = 0; i < d+1; ++i) J[i].resize(d+1);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        std::stringstream name;
	name << "jac" << i << j;
        J[i][j] = scalar_view(name.str(), c, fv, fv);
      }
      std::stringstream name1, name2;
      name1 << "jac" << d << i;
      name2 << "jac" << i << d;
      J[d][i] = scalar_view(name1.str(), c, fp, fv);
      J[i][d] = scalar_view(name2.str(), c, fv, fp);
    }
    std::stringstream name;
    name << "jac" << d << d;
    J[d][d] = scalar_view(name.str(), c, fp, fp);

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }

  void Hessian_11(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Stokes::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Stokes::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Stokes::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_Stokes::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(scalar_view & riesz) override {
    // Retrieve dimensions.
    int c  = feVel_->N().extent_int(0);
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<scalar_view>> J(d+1);
    for (int i = 0; i < d+1; ++i) J[i].resize(d+1);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        std::stringstream name;
	name << "jac" << i << j;
        J[i][j] = scalar_view(name.str(), c, fv, fv);
      }
      std::stringstream name1, name2;
      name1 << "jac" << d << i;
      name2 << "jac" << i << d;
      J[d][i] = scalar_view(name1.str(), c, fp, fv);
      J[i][d] = scalar_view(name2.str(), c, fv, fp);
    }
    std::stringstream name;
    name << "jac" << d << d;
    J[d][d] = scalar_view(name.str(), c, fp, fp);

    for (int i = 0; i < d; ++i) {
      Kokkos::deep_copy(J[i][i],feVel_->stiffMat());
      rst::add(J[i][i],feVel_->massMat());
    }
    Kokkos::deep_copy(J[d][d],fePrs_->massMat());

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  void RieszMap_2(scalar_view & riesz) override {
    // Retrieve dimensions.
    int c  = feVel_->N().extent_int(0);
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<scalar_view>> J(d+1);
    for (int i = 0; i < d+1; ++i) J[i].resize(d+1);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        std::stringstream name;
	name << "jac" << i << j;
        J[i][j] = scalar_view(name.str(), c, fv, fv);
      }
      std::stringstream name1, name2;
      name1 << "jac" << d << i;
      name2 << "jac" << i << d;
      J[d][i] = scalar_view(name1.str(), c, fp, fv);
      J[i][d] = scalar_view(name2.str(), c, fv, fp);
    }
    std::stringstream name;
    name << "jac" << d << d;
    J[d][d] = scalar_view(name.str(), c, fp, fp);

    for (int i = 0; i < d; ++i)
      Kokkos::deep_copy(J[i][i],feVel_->massMat());
    Kokkos::deep_copy(J[d][d],fePrs_->massMat());

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  std::vector<basis_ptr> getFields() override {
    return basisPtrs_;
  }

  void setCellNodes(const scalar_view &volCellNodes,
                    const std::vector<std::vector<scalar_view >> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds) override {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    feVel_ = ROL::makePtr<fe_type>(volCellNodes_,basisPtrVel_,cellCub_);
    fePrs_ = ROL::makePtr<fe_type>(volCellNodes_,basisPtrPrs_,cellCub_);
    // Get boundary degrees of freedom.
    fvidx_ = feVel_->getBoundaryDofs();
    fpidx_ = fePrs_->getBoundaryDofs();
    // Construct boundary FEs
    const int numSideSets = bdryCellNodes.size();
    feVelBdry_.resize(numSideSets);
    fePrsBdry_.resize(numSideSets);
    for (int i = 0; i < numSideSets; ++i) {
      int numLocSides = bdryCellNodes[i].size();
      feVelBdry_[i].resize(numLocSides);
      fePrsBdry_[i].resize(numLocSides);
      for (int j = 0; j < numLocSides; ++j) {
        if (bdryCellNodes[i][j] != scalar_view()) {
          feVelBdry_[i][j] = ROL::makePtr<fe_type>(bdryCellNodes[i][j],basisPtrVel_,bdryCub_,j);
          fePrsBdry_[i][j] = ROL::makePtr<fe_type>(bdryCellNodes[i][j],basisPtrPrs_,bdryCub_,j);
        }
      }
    }
    computeDirichlet();
  }

  void setFieldPattern(const std::vector<std::vector<int>> & fieldPattern) override {
    fieldPattern_ = fieldPattern;
    fieldHelper_ = ROL::makePtr<FieldHelper<Real,DeviceType>>(numFields_, numDofs_, numFieldDofs_, fieldPattern_);
  }

  const ROL::Ptr<fe_type> getVelocityFE(void) const {
    return feVel_;
  }

  const ROL::Ptr<fe_type> getPressureFE(void) const {
    return fePrs_;
  }

  const std::vector<std::vector<ROL::Ptr<fe_type>>> getVelocityBdryFE(void) const {
    return feVelBdry_;
  }

  const std::vector<std::vector<fe_type>> getPressureBdryFE(void) const {
    return fePrsBdry_;
  }

  const std::vector<std::vector<std::vector<int>>> getBdryCellLocIds(void) const {
    return bdryCellLocIds_;
  }

  const ROL::Ptr<FieldHelper<Real,DeviceType>> getFieldHelper(void) const {
    return fieldHelper_;
  }

}; // PDE_Stokes

#endif
