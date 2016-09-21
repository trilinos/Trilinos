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

/*! \file  obj.hpp
    \brief Provides the interface for local (cell-based) objective function computations.
*/

#ifndef PDEOPT_QOI_L2TRACKING_NAVIERSTOKES_HPP
#define PDEOPT_QOI_L2TRACKING_NAVIERSTOKES_HPP

#include "../TOOLS/qoi.hpp"
#include "pde_navier-stokes.hpp"

template <class Real>
class QoI_L2Tracking_NavierStokes : public QoI<Real> {
private:
  const Teuchos::RCP<FE<Real> > feVel_;
  const Teuchos::RCP<FE<Real> > fePrs_;
  const Teuchos::RCP<FieldHelper<Real> > fieldHelper_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > weight_;

  const Real eps_;

  Real weightFunc(const std::vector<Real> & x) const {
    return (((x[1] <= 0.5+eps_)&&(x[0] >= 1.0-eps_)&&(x[0] <= 4.0+eps_)) ?
                static_cast<Real>(1) : static_cast<Real>(0));
  }

public:
  QoI_L2Tracking_NavierStokes(const Teuchos::RCP<FE<Real> > &feVel,
                              const Teuchos::RCP<FE<Real> > &fePrs,
                              const Teuchos::RCP<FieldHelper<Real> > &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), fieldHelper_(fieldHelper), eps_(std::sqrt(ROL::ROL_EPSILON<Real>())) {
    int c = feVel_->cubPts()->dimension(0);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    std::vector<Real> pt(d);
    weight_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,p));
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (*feVel_->cubPts())(i,j,k);
        }
        (*weight_)(i,j) = weightFunc(pt);
      }
    }
  }

  Real value(Teuchos::RCP<Intrepid::FieldContainer<Real> > & val,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
             const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    int c = u_coeff->dimension(0);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    // Initialize output val
    val = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c));
    // Get components of the control
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradUX_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradUY_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    feVel_->evaluateGradient(gradUX_eval, U[0]);
    feVel_->evaluateGradient(gradUY_eval, U[1]);
    // Compute curl
    Teuchos::RCP<Intrepid::FieldContainer<Real> > curlU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*curlU_eval)(i,j)   = (*gradUY_eval)(i,j,0) - (*gradUX_eval)(i,j,1);
      }
    }
    // Multiply by weight
    Teuchos::RCP<Intrepid::FieldContainer<Real> > weighted_curlU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*weighted_curlU_eval,
                                                               *weight_,
                                                               *curlU_eval);
    // Compute L2 norm squared
    feVel_->computeIntegral(val,curlU_eval,weighted_curlU_eval,false);
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    int c = u_coeff->dimension(0);
    int fv = feVel_->N()->dimension(1);
    int fp = fePrs_->N()->dimension(1);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    // Initialize output grad
    Intrepid::FieldContainer<Real> velUX_grad(c, fv);
    Intrepid::FieldContainer<Real> velUY_grad(c, fv);
    Intrepid::FieldContainer<Real> presU_grad(c, fp);
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > G;
    G.resize(fieldHelper_->numFields());
    G[0] = Teuchos::rcpFromRef(velUX_grad);
    G[1] = Teuchos::rcpFromRef(velUY_grad);
    G[2] = Teuchos::rcpFromRef(presU_grad);
    // Get components of the control
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradUX_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradUY_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    feVel_->evaluateGradient(gradUX_eval, U[0]);
    feVel_->evaluateGradient(gradUY_eval, U[1]);
    // Compute curl
    Teuchos::RCP<Intrepid::FieldContainer<Real> > curlU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*curlU_eval)(i,j) = (*weight_)(i,j)*((*gradUY_eval)(i,j,0) - (*gradUX_eval)(i,j,1));
      }
    }
    // Build local gradient of state tracking term
    Intrepid::FunctionSpaceTools::integrate<Real>(velUX_grad,
                                                  *curlU_eval,
                                                  *(feVel_->DNDdetJ(1)),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(velUX_grad,static_cast<Real>(-1));
    Intrepid::FunctionSpaceTools::integrate<Real>(velUY_grad,
                                                  *curlU_eval,
                                                  *(feVel_->DNDdetJ(0)),
                                                  Intrepid::COMP_CPP, false);

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_L2Tracking_NavierStokes::gradient_2 is zero.");
  }

  void HessVec_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    int c  = z_coeff->dimension(0);
    int p  = feVel_->cubPts()->dimension(1);
    int fv = feVel_->N()->dimension(1);
    int fp = fePrs_->N()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    // Initialize output grad
    Intrepid::FieldContainer<Real> velVX_grad(c, fv);
    Intrepid::FieldContainer<Real> velVY_grad(c, fv);
    Intrepid::FieldContainer<Real> presV_grad(c, fp);
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > G;
    G.resize(fieldHelper_->numFields());
    G[0] = Teuchos::rcpFromRef(velVX_grad);
    G[1] = Teuchos::rcpFromRef(velVY_grad);
    G[2] = Teuchos::rcpFromRef(presV_grad);
    // Get components of the control
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Evaluate on FE basis
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradVX_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradVY_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    feVel_->evaluateGradient(gradVX_eval, V[0]);
    feVel_->evaluateGradient(gradVY_eval, V[1]);
    // Compute curl
    Teuchos::RCP<Intrepid::FieldContainer<Real> > curlV_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*curlV_eval)(i,j) = (*weight_)(i,j)*((*gradVY_eval)(i,j,0) - (*gradVX_eval)(i,j,1));
      }
    }
    // Build local gradient of state tracking term
    Intrepid::FunctionSpaceTools::integrate<Real>(velVX_grad,
                                                  *curlV_eval,
                                                  *(feVel_->DNDdetJ(1)),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(velVX_grad,static_cast<Real>(-1));
    Intrepid::FunctionSpaceTools::integrate<Real>(velVY_grad,
                                                  *curlV_eval,
                                                  *(feVel_->DNDdetJ(0)),
                                                  Intrepid::COMP_CPP, false);

    fieldHelper_->combineFieldCoeff(hess, G);
  }

  void HessVec_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_L2Tracking_NavierStokes::HessVec_12 is zero.");
  }

  void HessVec_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_L2Tracking_NavierStokes::HessVec_21 is zero.");
  }

  void HessVec_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_L2Tracking_NavierStokes::HessVec_22 is zero.");
  }

}; // QoI_L2Tracking_NavierStokes

template <class Real>
class QoI_L2Penalty_NavierStokes : public QoI<Real> {
private:
  const Teuchos::RCP<FE<Real> > feVel_;
  const Teuchos::RCP<FE<Real> > fePrs_;
  const std::vector<Teuchos::RCP<FE<Real> > > feVelBdry_;
  const std::vector<std::vector<int> > bdryCellLocIds_;
  const Teuchos::RCP<FieldHelper<Real> > fieldHelper_;

  Teuchos::RCP<Intrepid::FieldContainer<Real> > getBoundaryCoeff(
      const Intrepid::FieldContainer<Real> & cell_coeff,
      int locSideId) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[locSideId];
    const int numCellsSide = bdryCellLocId.size();
    const int f = feVel_->N()->dimension(1);
    
    Teuchos::RCP<Intrepid::FieldContainer<Real > > bdry_coeff = 
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, f));
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        (*bdry_coeff)(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

public:
  QoI_L2Penalty_NavierStokes(const Teuchos::RCP<FE<Real> > &feVel,
                             const Teuchos::RCP<FE<Real> > &fePrs,
                             const std::vector<Teuchos::RCP<FE<Real> > > & feVelBdry,
                             const std::vector<std::vector<int> > bdryCellLocIds,
                             const Teuchos::RCP<FieldHelper<Real> > &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), feVelBdry_(feVelBdry),
      bdryCellLocIds_(bdryCellLocIds), fieldHelper_(fieldHelper) {}

  Real value(Teuchos::RCP<Intrepid::FieldContainer<Real> > & val,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
             const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    const int c = feVel_->gradN()->dimension(0);
    // Get components of the control
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Initialize output val
    val = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c));
    // Compute cost integral
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = feVelBdry_[l]->cubPts()->dimension(1);
        // Evaluate x-component of control on FE basis
        Teuchos::RCP<Intrepid::FieldContainer<Real> > zx_coeff_bdry
          = getBoundaryCoeff(*Z[0], l);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > valZX_eval
          = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
        feVelBdry_[l]->evaluateValue(valZX_eval, zx_coeff_bdry);
        // Evaluate y-component of control on FE basis
        Teuchos::RCP<Intrepid::FieldContainer<Real> > zy_coeff_bdry
          = getBoundaryCoeff(*Z[1], l);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > valZY_eval
          = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
        feVelBdry_[l]->evaluateValue(valZY_eval, zy_coeff_bdry);
        // Integrate cell L2 norm squared
        Teuchos::RCP<Intrepid::FieldContainer<Real> > intVal
          = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide));
        feVelBdry_[l]->computeIntegral(intVal,valZX_eval,valZX_eval,false);
        feVelBdry_[l]->computeIntegral(intVal,valZY_eval,valZY_eval,true);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          (*val)(cidx) += static_cast<Real>(0.5)*(*intVal)(i);
        }
      }
    }
    return static_cast<Real>(0);
  }

  void gradient_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_L2Tracking_NavierStokes::gradient_1 is zero.");
  }

  void gradient_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    // Initialize output grad
    Intrepid::FieldContainer<Real> velZX_grad(c, fv);
    Intrepid::FieldContainer<Real> velZY_grad(c, fv);
    Intrepid::FieldContainer<Real> presZ_grad(c, fp);
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > G;
    G.resize(fieldHelper_->numFields());
    G[0] = Teuchos::rcpFromRef(velZX_grad);
    G[1] = Teuchos::rcpFromRef(velZY_grad);
    G[2] = Teuchos::rcpFromRef(presZ_grad);
    // Get components of the control
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Compute cost integral
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = feVelBdry_[l]->cubPts()->dimension(1);
        // Evaluate x-component of control on FE basis
        Teuchos::RCP<Intrepid::FieldContainer<Real> > zx_coeff_bdry
          = getBoundaryCoeff(*Z[0], l);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > valZX_eval
          = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
        feVelBdry_[l]->evaluateValue(valZX_eval, zx_coeff_bdry);
        // Evaluate y-component of control on FE basis
        Teuchos::RCP<Intrepid::FieldContainer<Real> > zy_coeff_bdry
          = getBoundaryCoeff(*Z[1], l);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > valZY_eval
          = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
        feVelBdry_[l]->evaluateValue(valZY_eval, zy_coeff_bdry);
        // Compute gradient of squared L2-norm of diff
        Teuchos::RCP<Intrepid::FieldContainer<Real> > intGradX
          = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, fv));
        Intrepid::FunctionSpaceTools::integrate<Real>(*intGradX,
                                                      *valZX_eval,
                                                      *(feVelBdry_[l]->NdetJ()),
                                                      Intrepid::COMP_CPP, false);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > intGradY
          = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, fv));
        Intrepid::FunctionSpaceTools::integrate<Real>(*intGradY,
                                                      *valZY_eval,
                                                      *(feVelBdry_[l]->NdetJ()),
                                                      Intrepid::COMP_CPP, false);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          for (int j = 0; j < fv; ++j) {
            (*G[0])(cidx,j) += (*intGradX)(i,j);
            (*G[1])(cidx,j) += (*intGradY)(i,j);
          }
        }
      }
    }

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void HessVec_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_L2Penalty_NavierStokes::HessVec_11 is zero.");
  }

  void HessVec_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_L2Penalty_NavierStokes::HessVec_12 is zero.");
  }

  void HessVec_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_L2Penalty_NavierStokes::HessVec_21 is zero.");
  }

  void HessVec_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    // Initialize output grad
    Intrepid::FieldContainer<Real> velVX_grad(c, fv);
    Intrepid::FieldContainer<Real> velVY_grad(c, fv);
    Intrepid::FieldContainer<Real> presV_grad(c, fp);
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > G;
    G.resize(fieldHelper_->numFields());
    G[0] = Teuchos::rcpFromRef(velVX_grad);
    G[1] = Teuchos::rcpFromRef(velVY_grad);
    G[2] = Teuchos::rcpFromRef(presV_grad);
    // Get components of the control
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Compute cost integral
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = feVelBdry_[l]->cubPts()->dimension(1);
        // Evaluate x-component of control on FE basis
        Teuchos::RCP<Intrepid::FieldContainer<Real> > vx_coeff_bdry
          = getBoundaryCoeff(*V[0], l);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > valVX_eval
          = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
        feVelBdry_[l]->evaluateValue(valVX_eval, vx_coeff_bdry);
        // Evaluate y-component of control on FE basis
        Teuchos::RCP<Intrepid::FieldContainer<Real> > vy_coeff_bdry
          = getBoundaryCoeff(*V[1], l);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > valVY_eval
          = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
        feVelBdry_[l]->evaluateValue(valVY_eval, vy_coeff_bdry);
        // Compute gradient of squared L2-norm of diff
        Teuchos::RCP<Intrepid::FieldContainer<Real> > intHessX
          = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, fv));
        Intrepid::FunctionSpaceTools::integrate<Real>(*intHessX,
                                                      *valVX_eval,
                                                      *(feVelBdry_[l]->NdetJ()),
                                                      Intrepid::COMP_CPP, false);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > intHessY
          = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, fv));
        Intrepid::FunctionSpaceTools::integrate<Real>(*intHessY,
                                                      *valVY_eval,
                                                      *(feVelBdry_[l]->NdetJ()),
                                                      Intrepid::COMP_CPP, false);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          for (int j = 0; j < fv; ++j) {
            (*G[0])(cidx,j) += (*intHessX)(i,j);
            (*G[1])(cidx,j) += (*intHessY)(i,j);
          }
        }
      }
    }

    fieldHelper_->combineFieldCoeff(hess, G);
  }

}; // QoI_L2Penalty_NavierStokes

template <class Real>
class StdObjective_NavierStokes : public ROL::ParametrizedStdObjective<Real> {
private:
  Real alpha_;

public:
  StdObjective_NavierStokes(Teuchos::ParameterList &parlist) {
    alpha_ = parlist.sublist("Problem").get("Control penalty parameter",1.e-4);
  }

  Real value(const std::vector<Real> &x, Real &tol) {
    return x[0] + alpha_*x[1];
  }

  void gradient(std::vector<Real> &g, const std::vector<Real> &x, Real &tol) {
    const Real one(1);
    g[0] = one;
    g[1] = alpha_;
  }

  void hessVec(std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    const Real zero(0);
    hv[0] = zero;
    hv[1] = zero;
  }

}; // OBJ_SCALAR

#endif
