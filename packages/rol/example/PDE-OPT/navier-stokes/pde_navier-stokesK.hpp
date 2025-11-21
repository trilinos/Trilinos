// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde_navier-stokes.hpp
    \brief Implements the local PDE interface for the Navier-Stokes control problem.
*/

#ifndef PDE_NAVIERSTOKESK_HPP
#define PDE_NAVIERSTOKESK_HPP

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
class PDE_NavierStokes : public PDE<Real,DeviceType> {
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
  std::vector<ROL::Ptr<fe_type>> feVelBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fvidx_, fpidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<scalar_view>> bdryCellDofValues_;
  // Field pattern, offsets, etc.
  std::vector<std::vector<int>> fieldPattern_; // local Field/DOF pattern; set from DOF manager 
  int numFields_;                              // number of fields (equations in the PDE)
  int numDofs_;                                // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                    // for each field, a counting offset
  std::vector<int> numFieldDofs_;              // for each field, number of degrees of freedom
  
  // Problem parameters.
  Real viscosity_;
  bool horizontalControl_;
  Real DirichletControlPenalty_;
  bool useDirichletControl_;

  ROL::Ptr<FieldHelper<Real,DeviceType>> fieldHelper_;

  // Extract velocity coefficients on boundary.
  scalar_view getBoundaryCoeff(const scalar_view cell_coeff, int sideSet, int cell) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][cell];
    const int numCellsSide = bdryCellLocId.size();
    const int f = basisPtrVel_->getCardinality();
    
    scalar_view bdry_coeff("bdry_coeff", numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j)
        bdry_coeff(i, j) = cell_coeff(bdryCellLocId[i], j);
    }
    return bdry_coeff;
  }

  Real DirichletFunc(const std::vector<Real> &coords, int sideset, int locSideId, int dir) const {
    const std::vector<Real> param = PDE<Real,DeviceType>::getParameter();
    Real val(0);
    if ((sideset==4) && (dir==0)) {
      if ( param.size() ) {
        Real zero(0), half(0.5), one(1), two(2), pi(M_PI), four(4), c1(1), c2(-0.5);
        if ( param[0] == zero )
          val = half*std::sin(two*pi * (c1*coords[1] + c2));
        else {
          Real num = (std::exp(four*param[0]*(c1*coords[1] + c2)) - one);
          Real den = (std::exp(two*param[0])-one);
          val = half*std::sin(pi * num/den);
        }
      }
      else {
        val = static_cast<Real>(8)
             *(coords[1]-static_cast<Real>(0.5))
             *(static_cast<Real>(1)-coords[1]);
      }
    }
    return val;
  }

  void computeDirichlet(void) {
    // Compute Dirichlet values at DOFs.
    int d = basisPtrVel_->getBaseCellTopology().getDimension();
    int numSidesets = bdryCellLocIds_.size();
    bdryCellDofValues_.resize(numSidesets);
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellDofValues_[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
        int f = basisPtrVel_->getCardinality();
        bdryCellDofValues_[i][j] = scalar_view("bdryCellDofValues", c, f, d);
        scalar_view coords("coords", c, f, d);
        if (c > 0) {
          feVel_->computeDofCoords(coords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<f; ++l) {
            std::vector<Real> dofpoint(d);
            //std::cout << "Sideset " << i << " LocalSide " << j << "  Cell " << k << "  Field " << l << "  Coord ";
            for (int m=0; m<d; ++m) {
              dofpoint[m] = coords(k, l, m);
              //std::cout << dofpoint[m] << "  ";
            }

            for (int m=0; m<d; ++m) {
              (bdryCellDofValues_[i][j])(k, l, m) = DirichletFunc(dofpoint, i, j, m);
              //std::cout << "  " << m << "-Value " << DirichletFunc(dofpoint, i, j, m);
            }
            //std::cout << std::endl;
          }
        }
      }
    }
  }

  Real viscosityFunc(const std::vector<Real> & coords) const {
    const std::vector<Real> param = PDE<Real,DeviceType>::getParameter();
    if ( param.size() ) {
      Real five(5), three(3);
      return (five + three*param[1])*static_cast<Real>(1e-3);
    }
    return viscosity_;
  }

  void computeViscosity(scalar_view &nu) const {
    int c = feVel_->gradN().extent_int(0);
    int p = feVel_->gradN().extent_int(2);
    int d = feVel_->gradN().extent_int(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k)
          pt[k] = (feVel_->cubPts())(i,j,k);
        // Compute spatially dependent viscosity
        nu(i,j) = viscosityFunc(pt);
      }
    }
  }

public:
  PDE_NavierStokes(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    basisPtrVel_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<DeviceType,Real,Real>>();
    basisPtrPrs_ = ROL::makePtr<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<DeviceType,Real,Real>>();
    basisPtrs_.clear();
    basisPtrs_.push_back(basisPtrVel_);  // Velocity X
    basisPtrs_.push_back(basisPtrVel_);  // Velocity Y
    basisPtrs_.push_back(basisPtrPrs_);  // Pressure
    // Quadrature rules.
    shards::CellTopology cellType = basisPtrs_[0]->getBaseCellTopology();    // get the cell type from any basis
    Intrepid2::DefaultCubatureFactory cubFactory;                             // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 4);    // set cubature degree, e.g., 4
    cellCub_ = cubFactory.create<DeviceType,Real,Real>(cellType, cubDegree); // create default cubature

    int d = cellType.getDimension();
    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",4); // set cubature degree, e.g., 4
    bdryCub_ = cubFactory.create<DeviceType,Real,Real>(bdryCellType, bdryCubDegree);

    // Other problem parameters.
    viscosity_ = parlist.sublist("Problem").get("Viscosity", 5e-3);
    horizontalControl_ = parlist.sublist("Problem").get("Horizontal control",false);
    DirichletControlPenalty_ = parlist.sublist("Problem").get("Dirichlet control penalty",1e-5);
    useDirichletControl_ = parlist.sublist("Problem").get("Dirichlet control",false);

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
    scalar_view velX_res("velX_res", c, fv);
    scalar_view velY_res("velY_res", c, fv);
    scalar_view pres_res("pres_res", c, fp);
    std::vector<scalar_view> R;
    R.resize(numFields_);
    R[0] = velX_res;
    R[1] = velY_res;
    R[2] = pres_res;

    // Split u_coeff into components.
    std::vector<scalar_view> U;
    std::vector<scalar_view> Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // Evaluate/interpolate finite element fields on cells.
    scalar_view nu("nu", c, p);
    scalar_view nuGradVelX_eval("nuGradVelX_eval", c, p, d);
    scalar_view nuGradVelY_eval("nuGradVelY_eval", c, p, d);
    scalar_view valVel_eval("valVel_eval", c, p, d);
    scalar_view valVelX_eval("valVelX_eval", c, p);
    scalar_view valVelY_eval("valVelY_eval", c, p);
    scalar_view valPres_eval("valPres_eval", c, p);
    scalar_view gradVelX_eval("gradVelX_eval", c, p, d);
    scalar_view gradVelY_eval("gradVelY_eval", c, p, d);
    scalar_view divVel_eval("divVel_eval", c, p);
    scalar_view valVelDotgradVelX_eval("valVelDotgradVelX_eval", c, p);
    scalar_view valVelDotgradVelY_eval("valVelDotgradVelY_eval", c, p);
    feVel_->evaluateValue(valVelX_eval, U[0]);
    feVel_->evaluateValue(valVelY_eval, U[1]);
    fePrs_->evaluateValue(valPres_eval, U[2]);
    feVel_->evaluateGradient(gradVelX_eval, U[0]);
    feVel_->evaluateGradient(gradVelY_eval, U[1]);
    computeViscosity(nu);

    // Assemble the velocity vector and its divergence.
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        valVel_eval(i,j,0) = valVelX_eval(i,j);
        valVel_eval(i,j,1) = valVelY_eval(i,j);
        divVel_eval(i,j)   = gradVelX_eval(i,j,0) + gradVelY_eval(i,j,1);
      }
    }
    // Negative pressure
    rst::scale(valPres_eval,static_cast<Real>(-1));

    // Multiply velocity gradients with viscosity.
    fst::tensorMultiplyDataData(nuGradVelX_eval,nu,gradVelX_eval);
    fst::tensorMultiplyDataData(nuGradVelY_eval,nu,gradVelY_eval);

    // Compute nonlinear terms in the Navier-Stokes equations.
    fst::dotMultiplyDataData(valVelDotgradVelX_eval, valVel_eval, gradVelX_eval);
    fst::dotMultiplyDataData(valVelDotgradVelY_eval, valVel_eval, gradVelY_eval);

    /*** Evaluate weak form of the residual. ***/
    // X component of velocity equation.
    fst::integrate(velX_res,nuGradVelX_eval,feVel_->gradNdetJ(),false);
    fst::integrate(velX_res,valVelDotgradVelX_eval,feVel_->NdetJ(),true);
    fst::integrate(velX_res,valPres_eval,feVel_->DNDdetJ(0),true);
    // Y component of velocity equation.
    fst::integrate(velY_res,nuGradVelY_eval,feVel_->gradNdetJ(),false);
    fst::integrate(velY_res,valVelDotgradVelY_eval,feVel_->NdetJ(),true);
    fst::integrate(velY_res,valPres_eval,feVel_->DNDdetJ(1),true);
    // Pressure equation.
    fst::integrate(pres_res,divVel_eval,fePrs_->NdetJ(),false);
    rst::scale(pres_res,static_cast<Real>(-1));

    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      if ( !useDirichletControl_ ) {
        const int numCubPerSide = bdryCub_->getNumPoints();
        for (int i = 0; i < numSideSets; ++i) {
          if (i==10) {
            int numLocalSideIds = bdryCellLocIds_[i].size();
            for (int j = 0; j < numLocalSideIds; ++j) {
              int numCellsSide = bdryCellLocIds_[i][j].size();
              if (numCellsSide) {
                for (int k = 0; k < d; ++k) {
                  // Get U coefficients on Robin boundary
                  scalar_view u_coeff_bdry = getBoundaryCoeff(U[k], i, j);
                  // Evaluate U on FE basis
                  scalar_view valU_eval_bdry("valU_eval_bdry", numCellsSide, numCubPerSide);
                  feVelBdry_[j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
                  // Get Z coefficients on Robin boundary
                  scalar_view z_coeff_bdry = getBoundaryCoeff(Z[k], i, j);
                  // Evaluate Z on FE basis
                  scalar_view valZ_eval_bdry("valZ_eval_bdry", numCellsSide, numCubPerSide);
                  if (horizontalControl_) {
                    if (k==0) // control in horizontal direction only
                      feVelBdry_[j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
                    // else valZ_eval_bdry=0
                  }
                  else // control in both directions
                    feVelBdry_[j]->evaluateValue(valZ_eval_bdry, z_coeff_bdry);
                  // Compute Robin residual
                  scalar_view robinVal("robinVal", numCellsSide, numCubPerSide);
                  rst::subtract(robinVal,valU_eval_bdry,valZ_eval_bdry);
                  rst::scale(robinVal,static_cast<Real>(1)/DirichletControlPenalty_);
                  scalar_view robinRes("robinRes", numCellsSide, fv);
                  fst::integrate(robinRes,robinVal,feVelBdry_[j]->NdetJ(),false);
                  // Add Robin residual to volume residual
                  for (int l = 0; l < numCellsSide; ++l) {
                    int cidx = bdryCellLocIds_[i][j][l];
                    for (int m = 0; m < fv; ++m)
                      (R[k])(cidx,m) += robinRes(l,m);
                  }
                }
              }
            }
          }
        }
      }
      // APPLY DIRICHLET CONDITIONS
      computeDirichlet();
      for (int i = 0; i < numSideSets; ++i) {
        if ( useDirichletControl_ ) {
          // Apply Dirichlet control
          if (i==6) {
            int numLocalSideIds = bdryCellLocIds_[i].size();
            for (int j = 0; j < numLocalSideIds; ++j) {
              int numCellsSide = bdryCellLocIds_[i][j].size();
              int numBdryDofs = fvidx_[j].size();
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < numBdryDofs; ++l) {
                  //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
                  if ( !horizontalControl_ ) {
                    for (int m=0; m < d; ++m) {
                      (R[m])(cidx,fvidx_[j][l]) = (U[m])(cidx,fvidx_[j][l]) - (Z[m])(cidx,fvidx_[j][l]);
                    }
                  }
                  else {
                    (R[0])(cidx,fvidx_[j][l]) = (U[0])(cidx,fvidx_[j][l]) - (Z[0])(cidx,fvidx_[j][l]);
                    (R[1])(cidx,fvidx_[j][l]) = (U[1])(cidx,fvidx_[j][l]);
                  }
                }
              }
            }
          }
          // Apply Dirichlet conditions
          if (i==8) {
            int numLocalSideIds = bdryCellLocIds_[i].size();
            for (int j = 0; j < numLocalSideIds; ++j) {
              int numCellsSide = bdryCellLocIds_[i][j].size();
              int numBdryDofs = fvidx_[j].size();
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < numBdryDofs; ++l) {
                  //std::cout << "\n  i=" << i << "   cidx=" << cidx << "   j=" << j << "  l=" << l << "  " << fvidx_[j][l] << " " << (*bdryCellDofValues_[i][j])(k,fvidx_[j][l],0);
                  //std::cout << "\n  i=" << i << "   cidx=" << cidx << "   j=" << j << "  l=" << l << "  " << fvidx_[j][l] << " " << (*bdryCellDofValues_[i][j])(k,fvidx_[j][l],1);
                  for (int m=0; m < d; ++m) {
                    (R[m])(cidx,fvidx_[j][l]) = (U[m])(cidx,fvidx_[j][l]) - (bdryCellDofValues_[i][j])(k,fvidx_[j][l],m);
                  }
                }
              }
            }
          }
        }
        // Apply Dirichlet conditions
        if ((i==0) || (i==3) || (i==4) || (i==5)) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                //std::cout << "\n  i=" << i << "   cidx=" << cidx << "   j=" << j << "  l=" << l << "  " << fvidx_[j][l] << " " << (*bdryCellDofValues_[i][j])(k,fvidx_[j][l],0);
                //std::cout << "\n  i=" << i << "   cidx=" << cidx << "   j=" << j << "  l=" << l << "  " << fvidx_[j][l] << " " << (*bdryCellDofValues_[i][j])(k,fvidx_[j][l],1);
                for (int m=0; m < d; ++m) {
                  (R[m])(cidx,fvidx_[j][l]) = (U[m])(cidx,fvidx_[j][l]) - (bdryCellDofValues_[i][j])(k,fvidx_[j][l],m);
                }
              }
            }
          }
        }
        // Step corner
        if (i==7) {
          int j = 0, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) { 
            int cidx = bdryCellLocIds_[i][0][0];
            for (int m=0; m < d; ++m) {
              (R[m])(cidx,fvidx_[j][l]) = (U[m])(cidx,fvidx_[j][l]) - (bdryCellDofValues_[i][j])(0,fvidx_[j][l],m);
            }
          }
        }
        // Pressure pinning
        if (i==9) {
          int j = 2, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) { 
            int cidx = bdryCellLocIds_[i][0][0];
            (R[d])(cidx,fpidx_[j][l]) = (U[d])(cidx,fpidx_[j][l]);
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
 
    // Initialize residuals.
    scalar_view velXvelX_jac("velXvelX_jac", c, fv, fv);
    scalar_view velXvelY_jac("velXvelY_jac", c, fv, fv);
    scalar_view velXpres_jac("velXpres_jac", c, fv, fp);
    scalar_view velYvelX_jac("velYvelX_jac", c, fv, fv);
    scalar_view velYvelY_jac("velYvelY_jac", c, fv, fv);
    scalar_view velYpres_jac("velYpres_jac", c, fv, fp);
    scalar_view presvelX_jac("presvelX_jac", c, fp, fv);
    scalar_view presvelY_jac("presvelY_jac", c, fp, fv);
    scalar_view prespres_jac("prespres_jac", c, fp, fp);
    std::vector<std::vector<scalar_view>> J;
    J.resize(numFields_);
    J[0].resize(numFields_);
    J[1].resize(numFields_);
    J[2].resize(numFields_);
    J[0][0] = velXvelX_jac; J[0][1] = velXvelY_jac; J[0][2] = velXpres_jac;  
    J[1][0] = velYvelX_jac; J[1][1] = velYvelY_jac; J[1][2] = velYpres_jac;  
    J[2][0] = presvelX_jac; J[2][1] = presvelY_jac; J[2][2] = prespres_jac;  

    // Split u_coeff into components.
    std::vector<scalar_view> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);

    // Evaluate/interpolate finite element fields on cells.
    scalar_view nu("nu", c, p);
    scalar_view nuGradPhiX_eval("nuGradPhiX", c, fv, p, d);
    scalar_view nuGradPhiY_eval("nuGradPhiY", c, fv, p, d);
    scalar_view valVel_eval("valVel_eval", c, p, d);
    scalar_view valVelX_eval("valVelX_eval", c, p);
    scalar_view valVelY_eval("valVelY_eval", c, p);
    scalar_view gradVelX_eval("gradVelX_eval", c, p, d);
    scalar_view ddxVelX_eval("ddxVelX_eval", c, p);
    scalar_view ddyVelX_eval("ddyVelX_eval", c, p);
    scalar_view ddxVelXPhiX_eval("ddxVelXPhiX_eval", c, fv, p);
    scalar_view ddyVelXPhiY_eval("ddyVelXPhiY_eval", c, fv, p);
    scalar_view gradVelY_eval("gradVelY_eval", c, p, d);
    scalar_view ddxVelY_eval("ddxVelY_eval", c, p);
    scalar_view ddyVelY_eval("ddyVelY_eval", c, p);
    scalar_view ddyVelYPhiY_eval("ddyVelYPhiY_eval", c, fv, p);
    scalar_view ddxVelYPhiX_eval("ddxVelYPhiX_eval", c, fv, p);
    scalar_view valVelDotgradPhiX_eval("valVelDotgradPhiX_eval", c, fv, p);
    scalar_view valVelDotgradPhiY_eval("valVelDotgradPhiY_eval", c, fv, p);
    feVel_->evaluateValue(valVelX_eval, U[0]);
    feVel_->evaluateValue(valVelY_eval, U[1]);
    feVel_->evaluateGradient(gradVelX_eval, U[0]);
    feVel_->evaluateGradient(gradVelY_eval, U[1]);
    computeViscosity(nu);

    // Assemble the velocity vector and its divergence.
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        valVel_eval(i,j,0) = valVelX_eval(i,j);
        valVel_eval(i,j,1) = valVelY_eval(i,j);
        ddxVelX_eval(i,j)  = gradVelX_eval(i,j,0);
        ddyVelX_eval(i,j)  = gradVelX_eval(i,j,1);
        ddxVelY_eval(i,j)  = gradVelY_eval(i,j,0);
        ddyVelY_eval(i,j)  = gradVelY_eval(i,j,1);
      }
    }

    // Multiply velocity gradients with viscosity.
    fst::tensorMultiplyDataField(nuGradPhiX_eval,nu,feVel_->gradN());
    fst::tensorMultiplyDataField(nuGradPhiY_eval,nu,feVel_->gradN());

    // Compute nonlinear terms in the Navier-Stokes equations.
    fst::dotMultiplyDataField(valVelDotgradPhiX_eval, valVel_eval, feVel_->gradN());
    fst::dotMultiplyDataField(valVelDotgradPhiY_eval, valVel_eval, feVel_->gradN());
    fst::scalarMultiplyDataField(ddxVelXPhiX_eval, ddxVelX_eval, feVel_->N());
    fst::scalarMultiplyDataField(ddyVelXPhiY_eval, ddyVelX_eval, feVel_->N());
    fst::scalarMultiplyDataField(ddxVelYPhiX_eval, ddxVelY_eval, feVel_->N());
    fst::scalarMultiplyDataField(ddyVelYPhiY_eval, ddyVelY_eval, feVel_->N());

    // Negative pressure basis.
    scalar_view negPrsN("negPrsN", c, fp, p);
    Kokkos::deep_copy(negPrsN,fePrs_->N());
    rst::scale(negPrsN,static_cast<Real>(-1));

    /*** Evaluate weak form of the Jacobian. ***/
    // X component of velocity equation.
    fst::integrate(velXvelX_jac,nuGradPhiX_eval,feVel_->gradNdetJ(),false);
    fst::integrate(velXvelX_jac,feVel_->NdetJ(),valVelDotgradPhiX_eval,true);
    fst::integrate(velXvelX_jac,feVel_->NdetJ(),ddxVelXPhiX_eval,true);
    fst::integrate(velXvelY_jac,feVel_->NdetJ(),ddyVelXPhiY_eval,false);
    fst::integrate(velXpres_jac,feVel_->DNDdetJ(0),negPrsN,false);
    // Y component of velocity equation.
    fst::integrate(velYvelY_jac,nuGradPhiY_eval,feVel_->gradNdetJ(),false);
    fst::integrate(velYvelY_jac,feVel_->NdetJ(),valVelDotgradPhiY_eval,true);
    fst::integrate(velYvelY_jac,feVel_->NdetJ(),ddyVelYPhiY_eval,true);
    fst::integrate(velYvelX_jac,feVel_->NdetJ(),ddxVelYPhiX_eval,false);
    fst::integrate(velYpres_jac,feVel_->DNDdetJ(1),negPrsN,false);
    // Pressure equation.
    fst::integrate(presvelX_jac,fePrs_->NdetJ(),feVel_->DND(0),false);
    fst::integrate(presvelY_jac,fePrs_->NdetJ(),feVel_->DND(1),false);
    rst::scale(presvelX_jac,static_cast<Real>(-1));
    rst::scale(presvelY_jac,static_cast<Real>(-1));

    // APPLY DIRICHLET CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      if (!useDirichletControl_) {
        for (int i = 0; i < numSideSets; ++i) {
          if (i==10) {
            int numLocalSideIds = bdryCellLocIds_[i].size();
            const int numCubPerSide = bdryCub_->getNumPoints();
            for (int j = 0; j < numLocalSideIds; ++j) {
              int numCellsSide = bdryCellLocIds_[i][j].size();
              if (numCellsSide) {
                for (int k = 0; k < d; ++k) {
                  // Compute Robin residual
                  scalar_view RobinDerivUN("RobinDerivUN", numCellsSide, fv, numCubPerSide);
                  rst::scale(RobinDerivUN,feVelBdry_[j]->N(),static_cast<Real>(1)/DirichletControlPenalty_);
                  scalar_view RobinJac("RobinJac", numCellsSide, fv, fv);
                  fst::integrate(RobinJac,RobinDerivUN,feVelBdry_[j]->NdetJ(),false);
                  // Add Stefan-Boltzmann residual to volume residual
                  for (int l = 0; l < numCellsSide; ++l) {
                    int cidx = bdryCellLocIds_[i][j][l];
                    for (int m = 0; m < fv; ++m) { 
                      for (int n = 0; n < fv; ++n) { 
                        (J[k][k])(cidx,m,n) += RobinJac(l,m,n);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      for (int i = 0; i < numSideSets; ++i) {
        if ( useDirichletControl_ ) {
          if ((i==6) || (i==8)) {
            int numLocalSideIds = bdryCellLocIds_[i].size();
            for (int j = 0; j < numLocalSideIds; ++j) {
              int numCellsSide = bdryCellLocIds_[i][j].size();
              int numBdryDofs = fvidx_[j].size();
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < numBdryDofs; ++l) {
                  for (int m=0; m < fv; ++m) {
                    for (int n=0; n < d; ++n) {
                      for (int q=0; q < d; ++q) {
                        (J[n][q])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                      }
                      (J[n][n])(cidx,fvidx_[j][l],fvidx_[j][l]) = static_cast<Real>(1);
                    }
                  }
                  for (int m=0; m < fp; ++m) {
                    for (int n=0; n < d; ++n) {
                      (J[n][2])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                    }
                  }
                }
              }
            }
          }
        }
        if ((i==0) || (i==3) || (i==4) || (i==5)) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < fv; ++m) {
                  for (int n=0; n < d; ++n) {
                    for (int q=0; q < d; ++q) {
                      (J[n][q])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                    }
                    (J[n][n])(cidx,fvidx_[j][l],fvidx_[j][l]) = static_cast<Real>(1);
                  }
                }
                for (int m=0; m < fp; ++m) {
                  for (int n=0; n < d; ++n) {
                    (J[n][2])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
              }
            }
          }
        }
        // Step corner
        if (i==7) {
          int j = 0, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) { 
            int cidx = bdryCellLocIds_[i][0][0];
            for (int m=0; m < fv; ++m) {
              for (int n=0; n < d; ++n) {
                for (int q=0; q < d; ++q) {
                  (J[n][q])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                }
                (J[n][n])(cidx,fvidx_[j][l],fvidx_[j][l]) = static_cast<Real>(1);
              }
            }
            for (int m=0; m < fp; ++m) {
              for (int n=0; n < d; ++n) {
                (J[n][d])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
              }
            }
          }
        }
        // Pressure pinning
        if (i==9) {
          int j = 2, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) { 
            int cidx = bdryCellLocIds_[i][0][0];
            for (int m = 0; m < fv; ++m) {
              for (int n = 0; n < d; ++n) {
                (J[d][n])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
              }
            }
            for (int m = 0; m < fp; ++m) {
              (J[d][d])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
            }
            (J[d][d])(cidx,fpidx_[j][l],fpidx_[j][l]) = static_cast<Real>(1);
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

    // Initialize residuals.
    scalar_view velXvelX_jac("velXvelX_jac", c, fv, fv);
    scalar_view velXvelY_jac("velXvelY_jac", c, fv, fv);
    scalar_view velXpres_jac("velXpres_jac", c, fv, fp);
    scalar_view velYvelX_jac("velYvelX_jac", c, fv, fv);
    scalar_view velYvelY_jac("velYvelY_jac", c, fv, fv);
    scalar_view velYpres_jac("velYpres_jac", c, fv, fp);
    scalar_view presvelX_jac("presvelX_jac", c, fp, fv);
    scalar_view presvelY_jac("presvelY_jac", c, fp, fv);
    scalar_view prespres_jac("prespres_jac", c, fp, fp);
    std::vector<std::vector<scalar_view>> J;
    J.resize(numFields_);
    J[0].resize(numFields_);
    J[1].resize(numFields_);
    J[2].resize(numFields_);
    J[0][0] = velXvelX_jac; J[0][1] = velXvelY_jac; J[0][2] = velXpres_jac;
    J[1][0] = velYvelX_jac; J[1][1] = velYvelY_jac; J[1][2] = velYpres_jac;
    J[2][0] = presvelX_jac; J[2][1] = presvelY_jac; J[2][2] = prespres_jac;

    // APPLY DIRICHLET CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      if (!useDirichletControl_) {
        for (int i = 0; i < numSideSets; ++i) {
          if (i==10) {
            int numLocalSideIds = bdryCellLocIds_[i].size();
            const int numCubPerSide = bdryCub_->getNumPoints();
            for (int j = 0; j < numLocalSideIds; ++j) {
              int numCellsSide = bdryCellLocIds_[i][j].size();
              if (numCellsSide) {
                for (int k = 0; k < d; ++k) {
                  // Compute Robin Jacobian
                  scalar_view RobinDerivZN("RobinDerivZN", numCellsSide, fv, numCubPerSide);
                  if (horizontalControl_) {
                    if (k==0) { // control in horizontal direction only
                      rst::scale(RobinDerivZN,feVelBdry_[j]->N(),static_cast<Real>(-1)/DirichletControlPenalty_);
                    }
                    // else RobinDerivZN=0
                  }
                  else { // control in both directions
                    rst::scale(RobinDerivZN,feVelBdry_[j]->N(),static_cast<Real>(-1)/DirichletControlPenalty_);
                  }
                  scalar_view RobinJac("RobinJac", numCellsSide, fv, fv);
                  fst::integrate(RobinJac,RobinDerivZN,feVelBdry_[j]->NdetJ(),false);
                  // Add Robin Jacobian to volume Jacobian
                  for (int l = 0; l < numCellsSide; ++l) {
                    int cidx = bdryCellLocIds_[i][j][l];
                    for (int m = 0; m < fv; ++m) { 
                      for (int n = 0; n < fv; ++n) { 
                        (J[k][k])(cidx,m,n) += RobinJac(l,m,n);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      for (int i = 0; i < numSideSets; ++i) {
        if ( useDirichletControl_ ) {
          // Apply Dirichlet controls
          if (i==6) {
            int numLocalSideIds = bdryCellLocIds_[i].size();
            for (int j = 0; j < numLocalSideIds; ++j) {
              int numCellsSide = bdryCellLocIds_[i][j].size();
              int numBdryDofs = fvidx_[j].size();
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < numBdryDofs; ++l) {
                  //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
                  for (int m=0; m < fv; ++m) {
                    if ( !horizontalControl_ ) {
                      for (int n=0; n < d; ++n) {
                        (J[n][n])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                      }
                    }
                    else {
                      (J[0][0])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                    }
                  }
                  if ( !horizontalControl_ ) {
                    for (int n=0; n < d; ++n) {
                      (J[n][n])(cidx,fvidx_[j][l],fvidx_[j][l]) = static_cast<Real>(-1);
                    }
                  }
                  else {
                    (J[0][0])(cidx,fvidx_[j][l],fvidx_[j][l]) = static_cast<Real>(-1);
                  }
                }
              }
            }
          }
          // Apply Dirichlet conditions
          if (i==8) {
            int numLocalSideIds = bdryCellLocIds_[i].size();
            for (int j = 0; j < numLocalSideIds; ++j) {
              int numCellsSide = bdryCellLocIds_[i][j].size();
              int numBdryDofs = fvidx_[j].size();
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < numBdryDofs; ++l) {
                  //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
                  for (int m=0; m < fv; ++m) {
                    for (int n=0; n < d; ++n) {
                      (J[n][n])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                    }
                  }
                }
              }
            }
          }
        }
        // Apply Dirichlet conditions
        if ((i==0) || (i==3) || (i==4) || (i==5)) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
                for (int m=0; m < fv; ++m) {
                  for (int n=0; n < d; ++n) {
                    (J[n][n])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
              }
            }
          }
        }
        // Step corner
        if (i==7) {
          int j = 0, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) { 
            int cidx = bdryCellLocIds_[i][0][0];
            for (int m=0; m < fv; ++m) {
              for (int n=0; n < d; ++n) {
                for (int p=0; p < d; ++p) {
                  (J[n][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                }
              }
            }
            for (int m=0; m < fp; ++m) {
              for (int n=0; n < d; ++n) {
                (J[n][d])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
              }
            }
          }
        }
        // Pressure pinning
        if (i==9) {
          int j = 2, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) { 
            int cidx = bdryCellLocIds_[i][0][0];
            for (int m = 0; m < fv; ++m) {
              for (int n = 0; n < d; ++n) {
                (J[d][n])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
              }
            }
            for (int m = 0; m < fp; ++m) {
              (J[d][d])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
            }
          }
        }
      }
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);

  }

  void Hessian_11(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
//    throw Exception::NotImplemented("");
    // Retrieve dimensions.
    int c  = u_coeff.extent_int(0);
    int p  = cellCub_->getNumPoints();
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize residuals.
    scalar_view velXvelX_jac("velXvelX_jac", c, fv, fv);
    scalar_view velXvelY_jac("velXvelY_jac", c, fv, fv);
    scalar_view velXpres_jac("velXpres_jac", c, fv, fp);
    scalar_view velYvelX_jac("velYvelX_jac", c, fv, fv);
    scalar_view velYvelY_jac("velYvelY_jac", c, fv, fv);
    scalar_view velYpres_jac("velYpres_jac", c, fv, fp);
    scalar_view presvelX_jac("presvelX_jac", c, fp, fv);
    scalar_view presvelY_jac("presvelY_jac", c, fp, fv);
    scalar_view prespres_jac("prespres_jac", c, fp, fp);
    std::vector<std::vector<scalar_view>> J;
    J.resize(numFields_);
    J[0].resize(numFields_);
    J[1].resize(numFields_);
    J[2].resize(numFields_);
    J[0][0] = velXvelX_jac; J[0][1] = velXvelY_jac; J[0][2] = velXpres_jac;  
    J[1][0] = velYvelX_jac; J[1][1] = velYvelY_jac; J[1][2] = velYpres_jac;  
    J[2][0] = presvelX_jac; J[2][1] = presvelY_jac; J[2][2] = prespres_jac;  

    // Split u_coeff into components.
    std::vector<scalar_view> L;
    fieldHelper_->splitFieldCoeff(L, l_coeff);

    // Apply Dirichlet conditions
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if ( useDirichletControl_ ) {
          if ((i==6) || (i==8)) {
            int numLocalSideIds = bdryCellLocIds_[i].size();
            for (int j = 0; j < numLocalSideIds; ++j) {
              int numCellsSide = bdryCellLocIds_[i][j].size();
              int numBdryDofs = fvidx_[j].size();
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < numBdryDofs; ++l) {
                  //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
                  for (int m=0; m < d; ++m) {
                    (L[m])(cidx,fvidx_[j][l]) = static_cast<Real>(0);
                  }
                }
              }
            }
          }
        }
        if ((i==0) || (i==3) || (i==4) || (i==5)) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
                for (int m=0; m < d; ++m) {
                  (L[m])(cidx,fvidx_[j][l]) = static_cast<Real>(0);
                }
              }
            }
          }
        }
        // Step corner
        if (i==7) {
          int j = 0, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) { 
            int cidx = bdryCellLocIds_[i][0][0];
            for (int m=0; m < d; ++m) {
              (L[m])(cidx,fvidx_[j][l]) = static_cast<Real>(0);
            }
          }
        }
        // Pressure pinning
        if (i==9) {
          int j = 2, l = 0;
          if (bdryCellLocIds_[i][0].size() > 0) { 
            int cidx = bdryCellLocIds_[i][0][0];
            (L[d])(cidx,fpidx_[j][l]) = static_cast<Real>(0);
          }
        }
      }
    }

    // Evaluate/interpolate finite element fields on cells.
    scalar_view valVelX_eval("valVelX_eval", c, p);
    scalar_view valVelY_eval("valVelY_eval", c, p);
    scalar_view valVelXPhi_eval("valVelXPhi_eval", c, fv, p);
    scalar_view valVelYPhi_eval("valVelYPhi_eval", c, fv, p);
    feVel_->evaluateValue(valVelX_eval, L[0]);
    feVel_->evaluateValue(valVelY_eval, L[1]);

    // Compute nonlinear terms in the Navier-Stokes equations.
    fst::scalarMultiplyDataField(valVelXPhi_eval, valVelX_eval, feVel_->N());
    fst::scalarMultiplyDataField(valVelYPhi_eval, valVelY_eval, feVel_->N());

    /*** Evaluate weak form of the Hessian. ***/
    // X component of velocity equation.
    fst::integrate(velXvelX_jac,valVelXPhi_eval,feVel_->DNDdetJ(0),false);
    fst::integrate(velXvelX_jac,feVel_->DNDdetJ(0),valVelXPhi_eval,true);
    fst::integrate(velXvelY_jac,feVel_->DNDdetJ(1),valVelXPhi_eval,false);
    fst::integrate(velXvelY_jac,valVelYPhi_eval,feVel_->DNDdetJ(0),true);
    // Y component of velocity equation.
    fst::integrate(velYvelY_jac,valVelYPhi_eval,feVel_->DNDdetJ(1),false);
    fst::integrate(velYvelY_jac,feVel_->DNDdetJ(1),valVelYPhi_eval,true);
    fst::integrate(velYvelX_jac,feVel_->DNDdetJ(0),valVelYPhi_eval,false);
    fst::integrate(velYvelX_jac,valVelXPhi_eval,feVel_->DNDdetJ(1),true);

    // Combine the Hessians.
    fieldHelper_->combineFieldCoeff(hess, J);
  }

  void Hessian_12(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_NavierStokes::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_NavierStokes::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(scalar_view & hess,
                  const scalar_view l_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> (PDE_NavierStokes::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(scalar_view & riesz) override {
    // Retrieve dimensions.
    int c  = feVel_->N().extent_int(0);
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
 
    // Initialize residuals.
    scalar_view velXvelX_jac("velXvelX_jac", c, fv, fv);
    scalar_view velXvelY_jac("velXvelY_jac", c, fv, fv);
    scalar_view velXpres_jac("velXpres_jac", c, fv, fp);
    scalar_view velYvelX_jac("velYvelX_jac", c, fv, fv);
    scalar_view velYvelY_jac("velYvelY_jac", c, fv, fv);
    scalar_view velYpres_jac("velYpres_jac", c, fv, fp);
    scalar_view presvelX_jac("presvelX_jac", c, fp, fv);
    scalar_view presvelY_jac("presvelY_jac", c, fp, fv);
    scalar_view prespres_jac("prespres_jac", c, fp, fp);
    std::vector<std::vector<scalar_view>> J;
    J.resize(numFields_);
    J[0].resize(numFields_);
    J[1].resize(numFields_);
    J[2].resize(numFields_);
    J[0][0] = velXvelX_jac; J[0][1] = velXvelY_jac; J[0][2] = velXpres_jac;  
    J[1][0] = velYvelX_jac; J[1][1] = velYvelY_jac; J[1][2] = velYpres_jac;  
    J[2][0] = presvelX_jac; J[2][1] = presvelY_jac; J[2][2] = prespres_jac;  

    Kokkos::deep_copy(J[0][0],feVel_->stiffMat());
    rst::add(J[0][0],feVel_->massMat());
    Kokkos::deep_copy(J[1][1],feVel_->stiffMat());
    rst::add(J[1][1],feVel_->massMat());
    Kokkos::deep_copy(J[2][2],fePrs_->massMat());

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  void RieszMap_2(scalar_view & riesz) override {
    // Retrieve dimensions.
    int c  = feVel_->N().extent_int(0);
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
 
    // Initialize residuals.
    scalar_view velXvelX_jac("velXvelX_jac", c, fv, fv);
    scalar_view velXvelY_jac("velXvelY_jac", c, fv, fv);
    scalar_view velXpres_jac("velXpres_jac", c, fv, fp);
    scalar_view velYvelX_jac("velYvelX_jac", c, fv, fv);
    scalar_view velYvelY_jac("velYvelY_jac", c, fv, fv);
    scalar_view velYpres_jac("velYpres_jac", c, fv, fp);
    scalar_view presvelX_jac("presvelX_jac", c, fp, fv);
    scalar_view presvelY_jac("presvelY_jac", c, fp, fv);
    scalar_view prespres_jac("prespres_jac", c, fp, fp);
    std::vector<std::vector<scalar_view>> J;
    J.resize(numFields_);
    J[0].resize(numFields_);
    J[1].resize(numFields_);
    J[2].resize(numFields_);
    J[0][0] = velXvelX_jac; J[0][1] = velXvelY_jac; J[0][2] = velXpres_jac;  
    J[1][0] = velYvelX_jac; J[1][1] = velYvelY_jac; J[1][2] = velYpres_jac;  
    J[2][0] = presvelX_jac; J[2][1] = presvelY_jac; J[2][2] = prespres_jac;  

    Kokkos::deep_copy(J[0][0],feVel_->stiffMat());
    rst::add(J[0][0],feVel_->massMat());
    Kokkos::deep_copy(J[1][1],feVel_->stiffMat());
    rst::add(J[1][1],feVel_->massMat());
    Kokkos::deep_copy(J[2][2],fePrs_->massMat());

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
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
    feVel_ = ROL::makePtr<fe_type>(volCellNodes_,basisPtrVel_,cellCub_);
    fePrs_ = ROL::makePtr<fe_type>(volCellNodes_,basisPtrPrs_,cellCub_);
    fvidx_ = feVel_->getBoundaryDofs();
    fpidx_ = fePrs_->getBoundaryDofs();
    // Construct control boundary FE
    int sideset = (useDirichletControl_) ? 6 : 10;
    int numLocSides = bdryCellNodes[sideset].size();
    feVelBdry_.resize(numLocSides);
    for (int j = 0; j < numLocSides; ++j) {
      if (bdryCellNodes[sideset][j] != scalar_view()) {
        feVelBdry_[j] = ROL::makePtr<fe_type>(bdryCellNodes[sideset][j],basisPtrVel_,bdryCub_,j);
      }
    }
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

  const std::vector<ROL::Ptr<fe_type>> getVelocityBdryFE(void) const {
    return feVelBdry_;
  }

  const std::vector<std::vector<int>> getBdryCellLocIds(const int sideset = -1) const {
    int side = sideset;
    if ( sideset < 0 ) {
      side = (useDirichletControl_) ? 6 : 10;
    }
    return bdryCellLocIds_[side];
  }

  const ROL::Ptr<FieldHelper<Real,DeviceType>> getFieldHelper(void) const {
    return fieldHelper_;
  }

}; // PDE_NavierStokes

#endif
