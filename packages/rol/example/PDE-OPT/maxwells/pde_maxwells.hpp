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

/*! \file  pde_maxwells.hpp
    \brief Implements the local PDE interface for the optimal control of
           Maxwells.
*/

#ifndef PDE_MAXWELLS_HPP
#define PDE_MAXWELLS_HPP

#include "../TOOLS/pde.hpp"
#include "../TOOLS/fe_curl.hpp"
#include "../TOOLS/fieldhelper.hpp"

#include "Intrepid_HCURL_HEX_I1_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"


template <class Real>
class PDE_Maxwells : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real> > cellCub_;
  ROL::Ptr<Intrepid::Cubature<Real> > bdryCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real> > volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;
  std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE_CURL<Real> > fe_;
  std::vector<ROL::Ptr<FE_CURL<Real> > > feBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int> > fidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellDofValues_;
  // Field pattern, offsets, etc.
  std::vector<std::vector<int> > fieldPattern_;  // local Field/DOF pattern; set from DOF manager 
  int numFields_;                                // number of fields (equations in the PDE)
  int numDofs_;                                  // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                      // for each field, a counting offset
  std::vector<int> numFieldDofs_;                // for each field, number of degrees of freedom
  ROL::Ptr<FieldHelper<Real> > fieldHelper_;

  void computeMuInv(const ROL::Ptr<Intrepid::FieldContainer<Real> > &muInv) const {
    int c = fe_->curlN()->dimension(0);
    int p = fe_->curlN()->dimension(2);
    int d = fe_->curlN()->dimension(3);
   
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (*fe_->cubPts())(i,j,k);
        }
        (*muInv)(i,j) = static_cast<Real>(1)/evaluateMu(x);
      }
    } 
  }

  void computeKappa(const ROL::Ptr<Intrepid::FieldContainer<Real> > &kappa, const int component) const {
    int c = fe_->curlN()->dimension(0);
    int p = fe_->curlN()->dimension(2);
    int d = fe_->curlN()->dimension(3);
   
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (*fe_->cubPts())(i,j,k);
        }
        (*kappa)(i,j) = evaluateKappa(x,component);
      }
    } 
  }

  void computeRHS(const ROL::Ptr<Intrepid::FieldContainer<Real> > &F, const int component) const {
    int c = fe_->curlN()->dimension(0);
    int p = fe_->curlN()->dimension(2);
    int d = fe_->curlN()->dimension(3);
   
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (*fe_->cubPts())(i,j,k);
        }
        for (int k = 0; k < d; ++k) {
          (*F)(i,j,k) = -evaluateRHS(x,k,component);
        }
      }
    }
  }

  ROL::Ptr<Intrepid::FieldContainer<Real> > getBoundaryCoeff(
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

public:
  PDE_Maxwells(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    basisPtr_ = ROL::makePtr<Intrepid::Basis_HCURL_HEX_I1_FEM<Real, Intrepid::FieldContainer<Real> >>();
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();            // get the cell type from the basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                           // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 2);        // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                           // create default cubature
    int d = cellType.getDimension();

    basisPtrs_.clear();
    for (int i=0; i<2; ++i) {
      basisPtrs_.push_back(basisPtr_);  // Displacement component
    }

    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2); // set cubature degree, e.g., 2
    bdryCub_ = cubFactory.create(bdryCellType, bdryCubDegree);

    numDofs_ = 0;
    numFields_ = basisPtrs_.size();
    offset_.resize(numFields_);
    numFieldDofs_.resize(numFields_);
    for (int i=0; i<numFields_; ++i) {
      if (i==0) {
        offset_[i]  = 0;
      }
      else {
        offset_[i]  = offset_[i-1] + basisPtrs_[i-1]->getCardinality();
      }
      numFieldDofs_[i] = basisPtrs_[i]->getCardinality();
      numDofs_ += numFieldDofs_[i];
    }
  }

  virtual Real evaluateMu(const std::vector<Real> &x) const {
    return static_cast<Real>(1);
  }

  virtual Real evaluateKappa(const std::vector<Real> &x, const int component) const {
    Real val(0);
    if (component == 0) {
      val = static_cast<Real>(1);
    }
    return val;
  }

  virtual Real evaluateRHS(const std::vector<Real> &x, const int d, const int component) const {
    const Real two(2), four(4);
    const Real pi = M_PI, pi2 = pi*pi;
    const Real mu = evaluateMu(x);
    const Real kappa = evaluateKappa(x, component);
    const Real cx = std::cos(pi*x[0]), sx = std::sin(pi*x[0]); 
    const Real cy = std::cos(pi*x[1]), sy = std::sin(pi*x[1]); 
    const Real cz = std::cos(pi*x[2]), sz = std::sin(pi*x[2]);
    const Real f  = sx*sy*sz;
    Real val(0);
    if (d == 0) {
      val = -pi2/mu * (cx*cy*sz + two*cx*sy*cz + two*f) + kappa*f;
    }
    else if (d == 1) {
      val = -pi2/mu * (cx*cy*sz + two*sx*cy*cz - two*f) - kappa*f;
    }
    else {
      val = -pi2/mu * (cx*sy*cz - sx*cy*cz + four*f) + two*kappa*f;
    }
    return val;
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real> > & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c = fe_->curlN()->dimension(0);
    int f = fe_->curlN()->dimension(1);
    int p = fe_->curlN()->dimension(2);
    int d = fe_->curlN()->dimension(3);
 
    // Initialize residuals.
    int nr = 2;
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > R(nr);
    for (int i=0; i<nr; ++i) {
      R[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f);
    }

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > valZ_eval(nr);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > valU_eval(nr);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > curlU_eval(nr);
    for (int i=0; i<nr; ++i) {
      valZ_eval[i]  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      fe_->evaluateValue(valZ_eval[i], Z[i]);
      valU_eval[i]  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      fe_->evaluateValue(valU_eval[i], U[i]);
      curlU_eval[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      fe_->evaluateCurl(curlU_eval[i], U[i]);
    }

    // Build 1/mu and 1/mu * curl U
    ROL::Ptr<Intrepid::FieldContainer<Real> > muInv
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    computeMuInv(muInv);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > muInvCurlU(nr);
    for (int i=0; i<nr; ++i) {
      muInvCurlU[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*muInvCurlU[i],*muInv,*curlU_eval[i]);
    }

    // Build kappa and kappa * U
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > kappa(nr);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > kappaU(2*nr);
    for (int i=0; i<nr; ++i) {
      kappa[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      computeKappa(kappa[i],i);
      kappaU[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*kappaU[i],*kappa[i],*valU_eval[i]);
      if (i==1) {
        Intrepid::RealSpaceTools<Real>::scale(*kappaU[i],static_cast<Real>(-1));
      }
      kappaU[i+nr] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*kappaU[i+nr],*kappa[i],*valU_eval[(i+1)%2]);
    }

    // Build right hand side
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > F(nr);
    for (int i=0; i<nr; ++i) {
      F[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      computeRHS(F[i],i);
    }

    /*******************************************************************/
    /*** Evaluate weak form of the residual.****************************/
    /*** a(ur,N) + br(ur,N) - bi(ui,N) = zr(N) + fr(N) *****************/
    /*** a(ui,N) + br(ui,N) + bi(ur,N) = zi(N) + fi(N) *****************/
    /*******************************************************************/
    for (int i=0; i<nr; ++i) {
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *muInvCurlU[i],      // 1/mu * curl U
                                                    *fe_->curlNdetJ(),   // curl N
                                                    Intrepid::COMP_CPP,
                                                    false);
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *kappaU[nr*i],       // kappa U
                                                    *fe_->NdetJ(),       // N
                                                    Intrepid::COMP_CPP,
                                                    true);
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *kappaU[nr*i+1],     // kappa U
                                                    *fe_->NdetJ(),       // N
                                                    Intrepid::COMP_CPP,
                                                    true);
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *F[i],              // -F
                                                    *fe_->NdetJ(),      // N
                                                    Intrepid::COMP_CPP,
                                                    true);
      Intrepid::RealSpaceTools<Real>::scale(*valZ_eval[i],static_cast<Real>(-1));
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *valZ_eval[i],      // -Z
                                                    *fe_->NdetJ(),      // N
                                                    Intrepid::COMP_CPP,
                                                    true);
    }

//    // APPLY ROBIN CONTROLS: Sideset 0
//    int sideset = 0;
//    int numLocalSideIds = bdryCellLocIds_[sideset].size();
//    const int numCubPerSide = bdryCub_->getNumPoints();
//    for (int j = 0; j < numLocalSideIds; ++j) {
//      int numCellsSide = bdryCellLocIds_[sideset][j].size();
//      if (numCellsSide) {
//        std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > robinRes(2);
//        for (int i = 0; i < 2; ++i) {
//          robinRes[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f);
//          // Get U coefficients on Robin boundary
//          ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff_bdry
//            = getBoundaryCoeff(*U[i], sideset, j);
//          // Evaluate U on FE basis
//          ROL::Ptr<Intrepid::FieldContainer<Real> > valU_eval_bdry
//            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
//          feBdry_[j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
//          // Compute Neumann residual
//          Intrepid::FunctionSpaceTools::integrate<Real>(*robinRes[i],
//                                                        *valU_eval_bdry,
//                                                        *(feBdry_[j]->NdetJ()),
//                                                        Intrepid::COMP_CPP, false);
//        }
//        // Add Neumann control residual to volume residual
//        for (int k = 0; k < numCellsSide; ++k) {
//          int cidx = bdryCellLocIds_[sideset][j][k];
//          for (int l = 0; l < f; ++l) { 
//            (*R[0])(cidx,l) += waveNumber_ * ((*robinRes[0])(k,l) + (*robinRes[1])(k,l));
//            (*R[1])(cidx,l) += waveNumber_ * ((*robinRes[1])(k,l) - (*robinRes[0])(k,l));
//          }
//        }
//      }
//    }

    // Combine the residuals.
    fieldHelper_->combineFieldCoeff(res, R);
  }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c = fe_->curlN()->dimension(0);
    int f = fe_->curlN()->dimension(1);
    int p = fe_->curlN()->dimension(2);
    int d = fe_->curlN()->dimension(3);
 
    // Initialize Jacobians.
    const int nr = 2;
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(nr);
    for (int i=0; i<nr; ++i) {
      for (int j=0; j<nr; ++j) {
        J[i].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f));
      }
    }

    // Build 1/mu and 1/mu * curl N
    ROL::Ptr<Intrepid::FieldContainer<Real> > muInv
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    computeMuInv(muInv);
    ROL::Ptr<Intrepid::FieldContainer<Real> > muInvCurlN
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p, d);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*muInvCurlN,*muInv,*fe_->curlN());

    // Build kappa and kappa * N
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > kappa(nr);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > kappaN(nr);
    for (int i=0; i<nr; ++i) {
      kappa[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      computeKappa(kappa[i],i);
      kappaN[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p, d);
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*kappaN[i],*kappa[i],*fe_->N());
    }

    /*** Evaluate weak form of the Jacobian. ***/
    // REAL-REAL BLOCK
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0],
                                                  *muInvCurlN,
                                                  *fe_->curlNdetJ(),
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0],
                                                  *kappaN[0],
                                                  *fe_->NdetJ(),
                                                  Intrepid::COMP_CPP,
                                                  true);
    // IMAG-IMAG BLOCK
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][1],
                                                  *muInvCurlN,
                                                  *fe_->curlNdetJ(),
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][1],
                                                  *kappaN[0],
                                                  *fe_->NdetJ(),
                                                  Intrepid::COMP_CPP,
                                                  true);
    // IMAG-REAL BLOCK
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][0],
                                                  *kappaN[1],
                                                  *fe_->NdetJ(),
                                                  Intrepid::COMP_CPP,
                                                  false);
    // REAL-IMAG BLOCK
    Intrepid::RealSpaceTools<Real>::scale(*kappaN[1],static_cast<Real>(-1));
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][1],
                                                 *kappaN[1],
                                                 *fe_->NdetJ(),
                                                 Intrepid::COMP_CPP,
                                                 false);

//    // APPLY ROBIN CONTROL: Sideset 0
//    int sideset = 0;
//    int numLocalSideIds = bdryCellLocIds_[sideset].size();
//    for (int j = 0; j < numLocalSideIds; ++j) {
//      int numCellsSide = bdryCellLocIds_[sideset][j].size();
//      if (numCellsSide) {
//        // Add Neumann control Jacobian to volume residual
//        for (int k = 0; k < numCellsSide; ++k) {
//          int cidx = bdryCellLocIds_[sideset][j][k];
//          for (int l = 0; l < f; ++l) { 
//            for (int m = 0; m < f; ++m) { 
//              (*J[0][0])(cidx,l,m) += waveNumber_ * (*feBdry_[j]->massMat())(k,l,m);
//              (*J[0][1])(cidx,l,m) += waveNumber_ * (*feBdry_[j]->massMat())(k,l,m);
//              (*J[1][0])(cidx,l,m) -= waveNumber_ * (*feBdry_[j]->massMat())(k,l,m);
//              (*J[1][1])(cidx,l,m) += waveNumber_ * (*feBdry_[j]->massMat())(k,l,m);
//            }
//          }
//        }
//      }
//    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }


  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c = fe_->curlN()->dimension(0);
    int f = fe_->curlN()->dimension(1);

    // Initialize Jacobians.
    const int nr = 2;
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(nr);
    for (int i=0; i<nr; ++i) {
      for (int j=0; j<nr; ++j) {
        J[i].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f));
      }
    }

    // Compute integral of -N * N
    ROL::Ptr<Intrepid::FieldContainer<Real> > negMass
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*negMass,
                                                 *fe_->N(),
                                                 *fe_->NdetJ(),
                                                 Intrepid::COMP_CPP,
                                                 false);
    Intrepid::RealSpaceTools<Real>::scale(*negMass,static_cast<Real>(-1));

    // REAL-REAL BLOCK
    *J[0][0] = *negMass;
    // IMAG-IMAG BLOCK
    *J[1][1] = *negMass;

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Maxwells::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Maxwells::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Maxwells::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Maxwells::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    //throw Exception::NotImplemented(">>> (PDE_TopoOpt::RieszMap_1): Not implemented.");

    // Retrieve dimensions.
    int c = fe_->curlN()->dimension(0);
    int f = fe_->curlN()->dimension(1);
 
    // Initialize Jacobians.
    const int nr = 2;
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(nr);
    for (int i=0; i<nr; ++i) {
      for (int j=0; j<nr; ++j) {
        J[i].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f));
      }
    }

    for (int i=0; i<nr; ++i) {
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][i],
                                                    *fe_->curlN(),
                                                    *fe_->curlNdetJ(),
                                                    Intrepid::COMP_CPP,
                                                    false);
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][i],
                                                    *fe_->N(),
                                                    *fe_->NdetJ(),
                                                    Intrepid::COMP_CPP,
                                                    true);
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    //throw Exception::NotImplemented(">>> (PDE_TopoOpt::RieszMap_2): Not implemented.");

    // Retrieve dimensions.
    int c = fe_->curlN()->dimension(0);
    int f = fe_->curlN()->dimension(1);
 
    // Initialize Jacobians.
    const int nr = 2;
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(nr);
    for (int i=0; i<nr; ++i) {
      for (int j=0; j<nr; ++j) {
        J[i].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f));
      }
    }

    for (int i=0; i<nr; ++i) {
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][i],
                                                    *fe_->N(),
                                                    *fe_->NdetJ(),
                                                    Intrepid::COMP_CPP,
                                                    false);
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > getFields() {
    return basisPtrs_;
  }

  void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real> > &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds) {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fe_ = ROL::makePtr<FE_CURL<Real>>(volCellNodes_,basisPtr_,cellCub_);
    fidx_ = fe_->getBoundaryDofs();
//    // Construct boundary FE
//    int sideset = 0;
//    int numLocSides = bdryCellNodes[sideset].size();
//    feBdry_.resize(numLocSides);
//    for (int j = 0; j < numLocSides; ++j) {
//      if (bdryCellNodes[sideset][j] != ROL::nullPtr) {
//        feBdry_[j] = ROL::makePtr<FE_CURL<Real>>(bdryCellNodes[sideset][j],basisPtr_,bdryCub_,j);
//      }
//    }
    // Compute control weight
    //computeControlWeight();
    //buildControlJacobian();
  }

  void setFieldPattern(const std::vector<std::vector<int> > & fieldPattern) {
    fieldPattern_ = fieldPattern;
    fieldHelper_ = ROL::makePtr<FieldHelper<Real>>(numFields_, numDofs_, numFieldDofs_, fieldPattern_);
  }

  const ROL::Ptr<FE_CURL<Real> > getFE(void) const {
    return fe_;
  }

  const std::vector<ROL::Ptr<FE_CURL<Real> > > getBdryFE(void) const {
    return feBdry_;
  }

  const std::vector<std::vector<int> > getBdryCellLocIds(const int sideset = 0) const {
    return bdryCellLocIds_[sideset];
  }

  const ROL::Ptr<FieldHelper<Real> > getFieldHelper(void) const {
    return fieldHelper_;
  }

}; // PDE_Maxwells


#endif
