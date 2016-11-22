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

#ifndef PDEOPT_QOI_GINZBURGLANDAU_HPP
#define PDEOPT_QOI_GINZBURGLANDAU_HPP

#include "../TOOLS/qoi.hpp"
#include "pde_ginzburg-landau.hpp"

template <class Real>
class QoI_GinzburgLandau_StateTracking : public QoI<Real> {
private:
  const Teuchos::RCP<FE<Real> > fe_;
  const Teuchos::RCP<FieldHelper<Real> > fieldHelper_;

  std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > target_;

  Real epsilon0_;
  Real lambda_;

protected:
  void computeTarget(void) {
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);

    target_.clear(); target_.resize(2);
    target_[0] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,p));
    target_[1] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,p));
   
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (*fe_->cubPts())(i,j,k);
        }
        (*target_[0])(i,j) = evaluateRealTarget(x);
        (*target_[1])(i,j) = evaluateImagTarget(x);
      }
    } 
  }
  
public:
  QoI_GinzburgLandau_StateTracking(const Teuchos::RCP<FE<Real> > &fe,
                                   const Teuchos::RCP<FieldHelper<Real> > &fieldHelper,
                                   Teuchos::ParameterList &parlist)
    : fe_(fe), fieldHelper_(fieldHelper) {
    lambda_   = parlist.sublist("Problem").get("Current Loading",1.0);
    epsilon0_ = parlist.sublist("Problem").get("State Scaling",1.0);
  }

  virtual Real evaluateRealTarget(const std::vector<Real> &x) const = 0;

  virtual Real evaluateImagTarget(const std::vector<Real> &x) const = 0;

  Real value(Teuchos::RCP<Intrepid::FieldContainer<Real> > & val,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
             const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    // Initialize output val
    val = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c));
    // Get components of the control
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    for (int i=0; i<2; ++i) {
      Teuchos::RCP<Intrepid::FieldContainer<Real> > valU_eval
        = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      fe_->evaluateValue(valU_eval, U[i]);
      Intrepid::RealSpaceTools<Real>::subtract(*valU_eval,*target_[i]);
      fe_->computeIntegral(val,valU_eval,valU_eval,true);
    }
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5)*lambda_/epsilon0_);
    return static_cast<Real>(0);
  }

  void gradient_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    // Initialize output grad
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > G(2);
    // Get components of the control
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    for (int i=0; i<2; ++i) {
      Teuchos::RCP<Intrepid::FieldContainer<Real> > valU_eval
        = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      fe_->evaluateValue(valU_eval, U[i]);
      Intrepid::RealSpaceTools<Real>::subtract(*valU_eval,*target_[i]);
      G[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
      Intrepid::FunctionSpaceTools::integrate<Real>(*G[i],
                                                    *valU_eval,
                                                    *(fe_->NdetJ()),
                                                    Intrepid::COMP_CPP,
                                                    false);
      Intrepid::RealSpaceTools<Real>::scale(*G[i],lambda_/epsilon0_);
    }

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_GinzburgLandau_StateTracking::gradient_2 is zero.");
  }

  void HessVec_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    // Initialize output hessvec
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > H(2);
    // Get components of the control
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Evaluate on FE basis
    for (int i=0; i<2; ++i) {
      Teuchos::RCP<Intrepid::FieldContainer<Real> > valV_eval
        = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      fe_->evaluateValue(valV_eval, V[i]);
      H[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
      Intrepid::FunctionSpaceTools::integrate<Real>(*H[i],
                                                    *valV_eval,
                                                    *(fe_->NdetJ()),
                                                    Intrepid::COMP_CPP,
                                                    false);
      Intrepid::RealSpaceTools<Real>::scale(*H[i],lambda_/epsilon0_);
    }

    fieldHelper_->combineFieldCoeff(hess, H);
  }

  void HessVec_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_TopoOpt::HessVec_12 is zero.");
  }

  void HessVec_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_TopoOpt::HessVec_21 is zero.");
  }

  void HessVec_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_TopoOpt::HessVec_22 is zero.");
  }

}; // QoI_TopoOpt


template <class Real>
class QoI_GinzburgLandau_ControlPenalty : public QoI<Real> {
private:
  const Teuchos::RCP<FE<Real> > fe_;
  const std::vector<Teuchos::RCP<FE<Real> > > feBdry_;
  const std::vector<std::vector<int> > bdryCellLocIds_;
  const Teuchos::RCP<FieldHelper<Real> > fieldHelper_;
  Real delta0_;
  Real lambda_;

  Teuchos::RCP<Intrepid::FieldContainer<Real> > getBoundaryCoeff(
      const Intrepid::FieldContainer<Real> & cell_coeff,
      int locSideId) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[locSideId];
    const int numCellsSide = bdryCellLocId.size();
    const int f = fe_->N()->dimension(1);
    
    Teuchos::RCP<Intrepid::FieldContainer<Real > > bdry_coeff = 
      Teuchos::rcp(new Intrepid::FieldContainer<Real > (numCellsSide, f));
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        (*bdry_coeff)(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

public:
  QoI_GinzburgLandau_ControlPenalty(const Teuchos::RCP<FE<Real> > &fe,
                                    const std::vector<Teuchos::RCP<FE<Real> > > &feBdry,
                                    const std::vector<std::vector<int> > &bdryCellLocIds,
                                    const Teuchos::RCP<FieldHelper<Real> > &fieldHelper,
                                    Teuchos::ParameterList &parlist)
  : fe_(fe), feBdry_(feBdry), bdryCellLocIds_(bdryCellLocIds), fieldHelper_(fieldHelper) {
    delta0_ = parlist.sublist("Problem").get("Control Scaling",1.0);
    lambda_ = parlist.sublist("Problem").get("Current Loading",1.0);
  }

  Real value(Teuchos::RCP<Intrepid::FieldContainer<Real> > & val,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
             const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    const int c = fe_->gradN()->dimension(0);
    // Initialize output val
    val = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c));
    // Get components of the control
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Compute cost integral
    const int numLocSides = bdryCellLocIds_.size();
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < numLocSides; ++j) {
        const int numCellsSide = bdryCellLocIds_[j].size();
        if ( numCellsSide ) {
          const int numCubPerSide = feBdry_[j]->cubPts()->dimension(1);
          // Evaluate control on FE basis
          Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff_bdry
            = getBoundaryCoeff(*Z[i], j);
          Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval
            = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
          feBdry_[j]->evaluateValue(valZ_eval, z_coeff_bdry);
          // Integrate cell cost
          Teuchos::RCP<Intrepid::FieldContainer<Real> > intVal
            = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide));
          feBdry_[j]->computeIntegral(intVal,valZ_eval,valZ_eval,false);
          // Add to integral value
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[j][k];
            (*val)(cidx) += static_cast<Real>(0.5)*delta0_*lambda_*(*intVal)(k);
          }
        }
      }
    }
    return static_cast<Real>(0);
  }

  void gradient_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_GinzburgLandau_ControlCost::gradient_1 is zero.");
  }

  void gradient_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    // Initialize output grad
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > G(2);
    for (int i=0; i<2; ++i) {
      G[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
    }
    // Get components of the control
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Compute cost integral
    const int numLocSides = bdryCellLocIds_.size();
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < numLocSides; ++j) {
        const int numCellsSide  = bdryCellLocIds_[j].size();
        if ( numCellsSide ) {
          const int numCubPerSide = feBdry_[j]->cubPts()->dimension(1);
          // Evaluate control on FE basis
          Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff_bdry
            = getBoundaryCoeff(*Z[i], j);
          Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval
            = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
          feBdry_[j]->evaluateValue(valZ_eval, z_coeff_bdry);
          // Compute gradient of squared L2-norm
          Intrepid::FieldContainer<Real> intGrad(numCellsSide, f);
          Intrepid::FunctionSpaceTools::integrate<Real>(intGrad,
                                                        *valZ_eval,
                                                        *(feBdry_[j]->NdetJ()),
                                                        Intrepid::COMP_CPP, false);
          // Add to integral value
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[j][k];
            for (int l = 0; l < f; ++l) {
              (*G[i])(cidx,l) += lambda_*delta0_*intGrad(k,l);
            }
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
    throw Exception::Zero(">>> QoI_GinzburgLandau_ControlCost::HessVec_11 is zero.");
  }

  void HessVec_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_GinzburgLandau_ControlCost::HessVec_12 is zero.");
  }

  void HessVec_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_GinzburgLandau_ControlCost::HessVec_21 is zero.");
  }

  void HessVec_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    // Initialize output grad
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > H(2);
    for (int i=0; i<2; ++i) {
      H[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
    }
    // Get components of the control
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Compute cost integral
    const int numLocSides = bdryCellLocIds_.size();
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < numLocSides; ++j) {
        const int numCellsSide  = bdryCellLocIds_[j].size();
        if ( numCellsSide ) {
          const int numCubPerSide = feBdry_[j]->cubPts()->dimension(1);
          // Evaluate control on FE basis
          Teuchos::RCP<Intrepid::FieldContainer<Real> > v_coeff_bdry
            = getBoundaryCoeff(*V[i], j);
          Teuchos::RCP<Intrepid::FieldContainer<Real> > valV_eval
            = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
          feBdry_[j]->evaluateValue(valV_eval, v_coeff_bdry);
          // Compute gradient of squared L2-norm of diff
          Intrepid::FieldContainer<Real> intHess(numCellsSide, f);
          Intrepid::FunctionSpaceTools::integrate<Real>(intHess,
                                                        *valV_eval,
                                                        *(feBdry_[j]->NdetJ()),
                                                        Intrepid::COMP_CPP, false);
          // Add to integral value
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[j][k];
            for (int l = 0; l < f; ++l) {
              (*H[i])(cidx,l) += lambda_*delta0_*intHess(k,l);
            }
          }
        }
      }
    }

    fieldHelper_->combineFieldCoeff(hess, H);
  }

}; // QoI_GinzburgLandau_ControlPenalty

#endif
