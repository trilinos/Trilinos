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

#ifndef PDEOPT_QOI_L2TRACKING_POISSON_BOLTZMANN_EX02_HPP
#define PDEOPT_QOI_L2TRACKING_POISSON_BOLTZMANN_EX02_HPP

#include "../TOOLS/qoi.hpp"
#include "doping.hpp"
#include "pde_poisson_boltzmann_ex02.hpp"

template <class Real>
class QoI_State_Cost_1_Poisson_Boltzmann : public QoI<Real> {
private:
  const Teuchos::RCP<FE<Real> > fe_vol_;
  const std::vector<std::vector<Teuchos::RCP<FE<Real> > > > fe_bdry_;
  const std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;

  Teuchos::RCP<Intrepid::FieldContainer<Real> > getBoundaryCoeff(
      const Intrepid::FieldContainer<Real> & cell_coeff,
      const int sideSet, const int locSideId) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][locSideId];
    const int numCellsSide = bdryCellLocId.size();
    const int f = fe_vol_->N()->dimension(1);
    
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
  QoI_State_Cost_1_Poisson_Boltzmann(const Teuchos::RCP<FE<Real> > &fe_vol,
                                     const std::vector<std::vector<Teuchos::RCP<FE<Real> > > > &fe_bdry,
                                     const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds,
                                     Teuchos::ParameterList &parlist)
    : fe_vol_(fe_vol), fe_bdry_(fe_bdry), bdryCellLocIds_(bdryCellLocIds) {}

  Real value(Teuchos::RCP<Intrepid::FieldContainer<Real> > & val,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
             const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    const int c = fe_vol_->gradN()->dimension(0);
    const int d = fe_vol_->gradN()->dimension(3);
    // Initialize output val
    val = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c));
    // Compute cost integral
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      if ( i == 3 ) {
        const int numLocSides = bdryCellLocIds_[i].size();
        for (int j = 0; j < numLocSides; ++j) {
          const int numCellsSide  = bdryCellLocIds_[i][j].size();
          if ( numCellsSide ) {
            const int numCubPerSide = fe_bdry_[i][j]->cubPts()->dimension(1);
            // Evaluate control on FE basis
            Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff_bdry
              = getBoundaryCoeff(*u_coeff, i, j);
            Teuchos::RCP<Intrepid::FieldContainer<Real> > gradU_eval
              = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide, d));
            fe_bdry_[i][j]->evaluateGradient(gradU_eval, u_coeff_bdry);
            // Compute cost
            Teuchos::RCP<Intrepid::FieldContainer<Real> > costU
              = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
            Teuchos::RCP<Intrepid::FieldContainer<Real> > weightU
              = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
            for (int k = 0; k < numCellsSide; ++k) {
              for (int l = 0; l < numCubPerSide; ++l) {
                (*costU)(k,l) = (*gradU_eval)(k,l,1);
                (*weightU)(k,l) = static_cast<Real>(1);
              }
            }
            // Integrate cell cost
            Teuchos::RCP<Intrepid::FieldContainer<Real> > intVal
              = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide));
            fe_bdry_[i][j]->computeIntegral(intVal,costU,weightU);
            // Add to integral value
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              (*val)(cidx) += (*intVal)(k);
            }
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
    const int c = fe_vol_->gradN()->dimension(0);
    const int f = fe_vol_->gradN()->dimension(1);
    // Initialize output val
    grad = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
    // Compute cost integral
    const int numSideSets = bdryCellLocIds_.size();
    for (int i = 0; i < numSideSets; ++i) {
      if ( i == 3 ) {
        const int numLocSides = bdryCellLocIds_[i].size();
        for (int j = 0; j < numLocSides; ++j) {
          const int numCellsSide  = bdryCellLocIds_[i][j].size();
          if ( numCellsSide ) {
            const int numCubPerSide = fe_bdry_[i][j]->cubPts()->dimension(1);
            // Compute cost
            Teuchos::RCP<Intrepid::FieldContainer<Real> > weightU
              = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, numCubPerSide));
            for (int k = 0; k < numCellsSide; ++k) {
              for (int l = 0; l < numCubPerSide; ++l) {
                (*weightU)(k,l) = static_cast<Real>(1);
              }
            }
            // Compute gradient of squared L2-norm of diff
            Teuchos::RCP<Intrepid::FieldContainer<Real> > intGrad
              = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCellsSide, f));
            Intrepid::FunctionSpaceTools::integrate<Real>(*intGrad,
                                                          *weightU,
                                                          *(fe_bdry_[i][j]->DNDdetJ(1)),
                                                          Intrepid::COMP_CPP, false);
            // Add to integral value
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < f; ++l) {
                (*grad)(cidx,l) += (*intGrad)(k,l);
              }
            }
          }
        }
      }
    }
  }

  void gradient_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_ControlCost::gradient_2 is zero.");
  }

  void HessVec_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_ControlCost::HessVec_11 is zero.");
  }

  void HessVec_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_ControlCost::HessVec_12 is zero.");
  }

  void HessVec_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_ControlCost::HessVec_21 is zero.");
  }

  void HessVec_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_ControlCost::HessVec_22 is zero.");
  }

}; // QoI_State_Cost_2_Poisson_Boltzmann

template <class Real>
class QoI_State_Cost_2_Poisson_Boltzmann : public QoI<Real> {
private:
  Teuchos::RCP<FE<Real> > fe_;

  Teuchos::RCP<Intrepid::FieldContainer<Real> > target_;

  Real targetFunc(const std::vector<Real> & x) const {
    return static_cast<Real>(1);
//    int size = x.size();
//    Real val(0);
//    for (int i = 0; i < size; ++i) {
//      val += x[i]*x[i];
//    }
//    return val;
  }

public:
  QoI_State_Cost_2_Poisson_Boltzmann(const Teuchos::RCP<FE<Real> > &fe) : fe_(fe) {
    int c = fe_->cubPts()->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int d = fe_->cubPts()->dimension(2);
    std::vector<Real> pt(d);
    target_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,p));
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (*fe_->cubPts())(i,j,k);
        }
        (*target_)(i,j) = targetFunc(pt);
      }
    }
  }

  Real value(Teuchos::RCP<Intrepid::FieldContainer<Real> > & val,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
             const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    int c = u_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    // Initialize output val
    val = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c));
    // Evaluate state on FE basis
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute difference between state and target
    Intrepid::RealSpaceTools<Real>::subtract(*valU_eval,*target_);
    // Compute squared L2-norm of diff
    fe_->computeIntegral(val,valU_eval,valU_eval);
    // Scale by one half
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    int c = u_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int f = fe_->N()->dimension(1);
    // Initialize output grad
    grad = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
    // Evaluate state on FE basis
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute difference between state and target
    Intrepid::RealSpaceTools<Real>::subtract(*valU_eval,*target_);
    // Compute gradient of squared L2-norm of diff
    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  *valU_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void gradient_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_State_Cost_Poisson_Boltzmann::gradient_2 is zero.");
  }

  void HessVec_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    int c = v_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int f = fe_->N()->dimension(1);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valV_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    hess = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
    fe_->evaluateValue(valV_eval, v_coeff);
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *valV_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void HessVec_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_State_Cost_Poisson_Boltzmann::HessVec_12 is zero.");
  }

  void HessVec_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_State_Cost_Poisson_Boltzmann::HessVec_21 is zero.");
  }

  void HessVec_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_State_Cost_Poisson_Boltzmann::HessVec_22 is zero.");
  }

}; // QoI_State_Cost

template <class Real>
class QoI_Control_Cost_1_Poisson_Boltzmann : public QoI<Real> {
private:
  const Teuchos::RCP<FE<Real> > fe_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > target_;

public:
  QoI_Control_Cost_1_Poisson_Boltzmann(const Teuchos::RCP<FE<Real> > &fe,
                                       const Teuchos::RCP<Doping<Real> > &dope)
    : fe_(fe) {
    // Compute reference doping profile
    const int c = fe_->gradN()->dimension(0);
    const int p = fe_->gradN()->dimension(2);
    const int d = fe_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    target_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,p));
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (*fe_->cubPts())(i,j,k);
        }
        (*target_)(i,j) = dope->evaluate(pt);
      }
    }
  }

  Real value(Teuchos::RCP<Intrepid::FieldContainer<Real> > & val,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
             const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int p = fe_->gradN()->dimension(2);
    // Initialize output val
    val = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c));
    // Evaluate state on FE basis
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_->evaluateValue(valZ_eval, z_coeff);
    // Compute difference between control and target
    Intrepid::RealSpaceTools<Real>::subtract(*valZ_eval,*target_);
    // Compute squared L2-norm of diff
    fe_->computeIntegral(val,valZ_eval,valZ_eval);
    // Scale by one half
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_Control_Cost_1_Poisson_Boltzmann::gradient_1 is zero.");
  }

  void gradient_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    // Initialize output grad
    grad = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
    // Evaluate state on FE basis
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_->evaluateValue(valZ_eval, z_coeff);
    // Compute difference between control and target
    Intrepid::RealSpaceTools<Real>::subtract(*valZ_eval,*target_);
    // Compute gradient of squared L2-norm of diff
    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  *valZ_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void HessVec_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_Control_Cost_1_Poisson_Boltzmann::HessVec_11 is zero.");
  }

  void HessVec_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_Control_Cost_1_Poisson_Boltzmann::HessVec_12 is zero.");
  }

  void HessVec_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_Control_Cost_1_Poisson_Boltzmann::HessVec_21 is zero.");
  }

  void HessVec_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valV_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    hess = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
    fe_->evaluateValue(valV_eval, v_coeff);
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *valV_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

}; // QoI_Control_Cost_1

template <class Real>
class QoI_Control_Cost_2_Poisson_Boltzmann : public QoI<Real> {
private:
  Teuchos::RCP<FE<Real> > fe_;

public:
  QoI_Control_Cost_2_Poisson_Boltzmann(const Teuchos::RCP<FE<Real> > &fe) : fe_(fe) {}

  Real value(Teuchos::RCP<Intrepid::FieldContainer<Real> > & val,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
             const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int p = fe_->gradN()->dimension(2);
    const int d = fe_->gradN()->dimension(3);
    // Initialize output val
    val = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c));
    // Evaluate state on FE basis
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradZ_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    fe_->evaluateGradient(gradZ_eval, z_coeff);
    // Compute squared L2-norm of grad
    fe_->computeIntegral(val,gradZ_eval,gradZ_eval);
    // Scale by one half
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_Control_Cost_2_Poisson_Boltzmann::gradient_1 is zero.");
  }

  void gradient_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    const int d = fe_->gradN()->dimension(3);
    // Initialize output grad
    grad = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
    // Evaluate state on FE basis
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradZ_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    fe_->evaluateGradient(gradZ_eval, z_coeff);
    // Compute gradient of squared L2-norm of gradient
    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  *gradZ_eval,
                                                  *(fe_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void HessVec_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_Control_Cost_2_Poisson_Boltzmann::HessVec_11 is zero.");
  }

  void HessVec_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_Control_Cost_2_Poisson_Boltzmann::HessVec_12 is zero.");
  }

  void HessVec_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_Control_Cost_2_Poisson_Boltzmann::HessVec_21 is zero.");
  }

  void HessVec_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    const int d = fe_->gradN()->dimension(3);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradV_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    hess = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
    fe_->evaluateGradient(gradV_eval, v_coeff);
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *gradV_eval,
                                                  *(fe_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

}; // QoI_Control_Cost_2

template <class Real>
class StdObjective_Poisson_Boltzmann : public ROL::StdObjective<Real> {
private:
  const Real J_, w1_, w2_, w3_;

public:
  StdObjective_Poisson_Boltzmann(const Real J, const Real w1, const Real w2, const Real w3)
    : J_(J), w1_(w1), w2_(w2), w3_(w3) {}

  Real value(const std::vector<Real> &x, Real &tol) {
    const Real half(0.5);
    return half*w1_*std::pow(x[0]-J_,2) + w2_*x[1] + w3_*x[2];
  }

  void gradient(std::vector<Real> &g, const std::vector<Real> &x, Real &tol) {
    g[0] = w1_*(x[0]-J_);
    g[1] = w2_;
    g[2] = w3_;
  }

  void hessVec(std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    const Real zero(0);
    hv[0] = w1_*v[0];
    hv[1] = zero;
    hv[2] = zero;
  }

}; // OBJ_SCALAR

#endif
