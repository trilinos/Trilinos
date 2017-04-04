// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  pde_stokes.hpp
    \brief Implements the local PDE interface for the Stokes control problem.
*/

#ifndef PDE_STOKES_HPP
#define PDE_STOKES_HPP

#include "../TOOLS/pde.hpp"
#include "../TOOLS/fe.hpp"
#include "../TOOLS/fieldhelper.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "Teuchos_RCP.hpp"

template <class Real>
class PDE_Stokes : public PDE<Real> {
private:
  // Finite element basis information
  Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtrVel_;
  Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtrPrs_;
  std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;
  // Cell cubature information
  Teuchos::RCP<Intrepid::Cubature<Real> > cellCub_;
  Teuchos::RCP<Intrepid::Cubature<Real> > bdryCub_;
  // Cell node information
  Teuchos::RCP<Intrepid::FieldContainer<Real> > volCellNodes_;
  std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;
  std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
  // Finite element definition
  Teuchos::RCP<FE<Real> > feVel_;
  Teuchos::RCP<FE<Real> > fePrs_;
  std::vector<std::vector<Teuchos::RCP<FE<Real> > > > feVelBdry_;
  std::vector<std::vector<Teuchos::RCP<FE<Real> > > > fePrsBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int> > fvidx_;
  std::vector<std::vector<int> > fpidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > bdryCellVDofValues_;
  // Field pattern, offsets, etc.
  std::vector<std::vector<int> > fieldPattern_;  // local Field/DOF pattern; set from DOF manager 
  int numFields_;                                // number of fields (equations in the PDE)
  int numDofs_;                                  // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                      // for each field, a counting offset
  std::vector<int> numFieldDofs_;                // for each field, number of degrees of freedom
  
  // Problem parameters.
  Real Re_;
  bool pinPressure_, dirType_;

  Teuchos::RCP<FieldHelper<Real> > fieldHelper_;

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
        bdryCellVDofValues_[i][j] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, d));
        Teuchos::RCP<Intrepid::FieldContainer<Real> > Vcoords, Tcoords;
        Vcoords = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, d));
        if (c > 0) {
          feVel_->computeDofCoords(Vcoords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<fv; ++l) {
            std::vector<Real> dofpoint(d);
            //std::cout << "Sideset " << i << " LocalSide " << j << "  Cell " << k << "  Field " << l << "  Coord ";
            for (int m=0; m<d; ++m) {
              dofpoint[m] = (*Vcoords)(k, l, m);
              //std::cout << dofpoint[m] << "  ";
            }

            for (int m=0; m<d; ++m) {
              (*bdryCellVDofValues_[i][j])(k, l, m) = velocityDirichletFunc(dofpoint, i, j, m);
              //std::cout << "  " << m << "-Value " << DirichletFunc(dofpoint, i, j, m);
            }
            //std::cout << std::endl;
          }
        }
      }
    }
  }

public:
  PDE_Stokes(Teuchos::ParameterList &parlist) {
    // Finite element fields -- NOT DIMENSION INDEPENDENT!
    basisPtrVel_ = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >);
    basisPtrPrs_ = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >);
    // Volume quadrature rules.
    shards::CellTopology cellType = basisPtrVel_->getBaseCellTopology();         // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                           // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 4);        // set cubature degree, e.g., 4
    cellCub_ = cubFactory.create(cellType, cubDegree);                           // create default cubature
    // Boundary quadrature rules.
    int d = cellType.getDimension();
    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",4); // set cubature degree, e.g., 4
    bdryCub_ = cubFactory.create(bdryCellType, bdryCubDegree);
    // Fill finite element basis container
    basisPtrs_.clear();
    for (int i = 0; i < d; ++i) {
      basisPtrs_.push_back(basisPtrVel_); // Velocity
    }
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
      if (i==0) {
        offset_[i] = 0;
      }
      else {
        offset_[i] = offset_[i-1] + basisPtrs_[i-1]->getCardinality();
      }
      numFieldDofs_[i] = basisPtrs_[i]->getCardinality();
      numDofs_ += numFieldDofs_[i];
    }
  }

  void residual(Teuchos::RCP<Intrepid::FieldContainer<Real> > & res,
                const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Retrieve dimensions.
    int c  = u_coeff->dimension(0);
    int p  = cellCub_->getNumPoints();
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize residuals.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > R(d+1);
    for (int i = 0; i < d; ++i) {
      R[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv));
    }
    R[d]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp));

    // Split u_coeff into components.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);

    // Evaluate/interpolate finite element fields on cells.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valPres_eval;
    valPres_eval = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fePrs_->evaluateValue(valPres_eval, U[d]);
    // Evaluate/interpolate gradient of finite element fields on cells.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > gradVel_vec(d);
    for (int i = 0; i < d; ++i) {
      gradVel_vec[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
      feVel_->evaluateGradient(gradVel_vec[i], U[i]);
    }

    // Assemble the velocity vector and its divergence.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > divVel_eval;
    divVel_eval = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*divVel_eval)(i,j) = static_cast<Real>(0);
        for (int k = 0; k < d; ++k) {
          (*divVel_eval)(i,j)  += (*gradVel_vec[k])(i,j,k);
        }
      }
    }

    // Scale pressure and velocity gradient
    const Real zero(0), one(1);
    for (int i = 0; i < d; ++i) {
      // Multiply velocity gradients with viscosity.
      Intrepid::RealSpaceTools<Real>::scale(*gradVel_vec[i], one/Re_);
    }
    // Negative pressure
    Intrepid::RealSpaceTools<Real>::scale(*valPres_eval, -one);

    /**************************************************************************/
    /*** EVALUATE WEAK FORM OF RESIDUAL ***************************************/
    /**************************************************************************/
    // Velocity equation.
    for (int i = 0; i < d; ++i) {
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *gradVel_vec[i],        // nu gradUX
                                                    *(feVel_->gradNdetJ()), // gradPhi
                                                    Intrepid::COMP_CPP,
                                                    false);
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *valPres_eval,          // p
                                                    *(feVel_->DNDdetJ(i)),  // dPhi/dx
                                                    Intrepid::COMP_CPP,
                                                    true);
    }
    // Pressure equation.
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[d],
                                                  *divVel_eval,             // divU
                                                  *(fePrs_->NdetJ()),       // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::RealSpaceTools<Real>::scale(*R[d], -one);

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
                  (*R[m])(cidx,fvidx_[j][l]) = (*U[m])(cidx,fvidx_[j][l]) - (*bdryCellVDofValues_[i][j])(k,fvidx_[j][l],m);
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
              (*R[d])(cidx,fpidx_[j][l]) = (*U[d])(cidx,fpidx_[j][l]) - val;
            }
          }
        }
      }
    }

    // Combine the residuals.
    fieldHelper_->combineFieldCoeff(res, R);
  }

  void Jacobian_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & jac,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Retrieve dimensions.
    int c  = u_coeff->dimension(0);
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(d+1);
    for (int i = 0; i < d+1; ++i) {
      J[i].resize(d+1);
    }
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        J[i][j] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fv));
      }
      J[d][i]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fv));
      J[i][d]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fp));
    }
    J[d][d]     = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fp));

    // Multiply velocity gradients with viscosity.
    const Real zero(0), one(1);
    Intrepid::FieldContainer<Real> nuGradPhi_eval(*feVel_->gradN());
    Intrepid::RealSpaceTools<Real>::scale(nuGradPhi_eval, one/Re_);
    // Negative pressure basis.
    Intrepid::FieldContainer<Real> negPrsPhi(*fePrs_->N());
    Intrepid::RealSpaceTools<Real>::scale(negPrsPhi, -one);

    /*** Evaluate weak form of the Jacobian. ***/
    for (int i = 0; i < d; ++i) {
      // Velocity components
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][i],
                                                    nuGradPhi_eval,       // nu gradPhi
                                                    *feVel_->gradNdetJ(), // gradPhi
                                                    Intrepid::COMP_CPP,
                                                    true);
      // Pressure components
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][d],
                                                    *feVel_->DNDdetJ(i),  // dPhi/dx
                                                    negPrsPhi,            // -Phi
                                                    Intrepid::COMP_CPP,
                                                    false);
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[d][i],
                                                    *fePrs_->NdetJ(),     // Phi
                                                    *feVel_->DND(i),      // dPhi/dx
                                                    Intrepid::COMP_CPP,
                                                    false);
      Intrepid::RealSpaceTools<Real>::scale(*J[d][i], -one);
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
                    for (int p = 0; p < d; ++p) {
                      (*J[n][p])(cidx,fvidx_[j][l],m) = zero;
                    }
                    (*J[n][n])(cidx,fvidx_[j][l],fvidx_[j][l]) = one;
                  }
                }
                for (int m = 0; m < d; ++m) {
                  for (int n = 0; n < fp; ++n) {
                    (*J[m][d])(cidx,fvidx_[j][l],n) = zero;
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
                  (*J[d][n])(cidx,fpidx_[j][l],m) = zero;
                }
              }
              for (int m = 0; m < fp; ++m) {
                (*J[d][d])(cidx,fpidx_[j][l],m) = zero;
              }
              (*J[d][d])(cidx,fpidx_[j][l],fpidx_[j][l]) = one;
            }
          }
        }
      }
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);

  }


  void Jacobian_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & jac,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Retrieve dimensions.
    int c  = u_coeff->dimension(0);
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(d+1);
    for (int i = 0; i < d+1; ++i) {
      J[i].resize(d+1);
    }
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        J[i][j] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fv));
      }
      J[d][i]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fv));
      J[i][d]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fp));
    }
    J[d][d]     = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fp));

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);

  }

  void Hessian_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Stokes::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Stokes::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Stokes::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Stokes::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & riesz) {
    // Retrieve dimensions.
    int c  = feVel_->N()->dimension(0);
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(d+1);
    for (int i = 0; i < d+1; ++i) {
      J[i].resize(d+1);
    }
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        J[i][j] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fv));
      }
      J[d][i]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fv));
      J[i][d]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fp));
    }
    J[d][d]     = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fp));

    for (int i = 0; i < d; ++i) {
      *(J[i][i]) = *(feVel_->stiffMat());
      Intrepid::RealSpaceTools<Real>::add(*(J[i][i]),*(feVel_->massMat()));
    }
    *(J[d][d]) = *(fePrs_->massMat());

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  void RieszMap_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & riesz) {
    // Retrieve dimensions.
    int c  = feVel_->N()->dimension(0);
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(d+1);
    for (int i = 0; i < d+1; ++i) {
      J[i].resize(d+1);
    }
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        J[i][j] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fv));
      }
      J[d][i]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fv));
      J[i][d]   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fv, fp));
    }
    J[d][d]     = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, fp, fp));

    for (int i = 0; i < d; ++i) {
      *(J[i][i]) = *(feVel_->massMat());
    }
    *(J[d][d]) = *(fePrs_->massMat());

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > getFields() {
    return basisPtrs_;
  }

  void setCellNodes(const Teuchos::RCP<Intrepid::FieldContainer<Real> > &volCellNodes,
                    const std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds) {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    feVel_ = Teuchos::rcp(new FE<Real>(volCellNodes_,basisPtrVel_,cellCub_));
    fePrs_ = Teuchos::rcp(new FE<Real>(volCellNodes_,basisPtrPrs_,cellCub_));
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
        if (bdryCellNodes[i][j] != Teuchos::null) {
          feVelBdry_[i][j] = Teuchos::rcp(new FE<Real>(bdryCellNodes[i][j],basisPtrVel_,bdryCub_,j));
          fePrsBdry_[i][j] = Teuchos::rcp(new FE<Real>(bdryCellNodes[i][j],basisPtrPrs_,bdryCub_,j));
        }
      }
    }
    computeDirichlet();
  }

  void setFieldPattern(const std::vector<std::vector<int> > & fieldPattern) {
    fieldPattern_ = fieldPattern;
    fieldHelper_ = Teuchos::rcp(new FieldHelper<Real>(numFields_, numDofs_, numFieldDofs_, fieldPattern_));
  }

  const Teuchos::RCP<FE<Real> > getVelocityFE(void) const {
    return feVel_;
  }

  const Teuchos::RCP<FE<Real> > getPressureFE(void) const {
    return fePrs_;
  }

  const std::vector<std::vector<Teuchos::RCP<FE<Real> > > > getVelocityBdryFE(void) const {
    return feVelBdry_;
  }

  const std::vector<std::vector<Teuchos::RCP<FE<Real> > > > getPressureBdryFE(void) const {
    return fePrsBdry_;
  }

  const std::vector<std::vector<std::vector<int> > > getBdryCellLocIds(void) const {
    return bdryCellLocIds_;
  }

  const Teuchos::RCP<FieldHelper<Real> > getFieldHelper(void) const {
    return fieldHelper_;
  }

}; // PDE_Stokes

#endif
