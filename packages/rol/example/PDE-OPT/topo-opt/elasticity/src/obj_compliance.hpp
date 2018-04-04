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

#ifndef PDEOPT_QOI_TOPOOPT_HPP
#define PDEOPT_QOI_TOPOOPT_HPP

#include "../../../TOOLS/qoi.hpp"
#include "pde_elasticity.hpp"

template <class Real>
class QoI_Compliance_TopoOpt : public QoI<Real> {
private:
  // Volumetric Information
  const ROL::Ptr<FE<Real> > fe_;
  const ROL::Ptr<Load<Real> > load_;
  // Boundary Information
  const std::vector<std::vector<ROL::Ptr<FE<Real> > > > feBdry_;
  const std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
  const ROL::Ptr<Traction<Real> > traction_;
  // Field Information
  const ROL::Ptr<FieldHelper<Real> > fieldHelper_;
  // QoI Scaling
  const Real scale_;

  ROL::Ptr<Intrepid::FieldContainer<Real> > getBoundaryCoeff(
      const Intrepid::FieldContainer<Real> & cell_coeff,
      const int sideset, const int locSideId) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideset][locSideId];
    const int numCellsSide = bdryCellLocId.size();
    const int f = fe_->N()->dimension(1);
    
    ROL::Ptr<Intrepid::FieldContainer<Real > > bdry_coeff = 
      ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        (*bdry_coeff)(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

public:
  QoI_Compliance_TopoOpt(const ROL::Ptr<FE<Real>>          &fe,
                         const ROL::Ptr<Load<Real>>        &load,
                         const ROL::Ptr<FieldHelper<Real>> &fieldHelper,
                         const Real scale = 1.0)
    : fe_(fe), load_(load), traction_(ROL::makePtr<Traction<Real>>()),
      fieldHelper_(fieldHelper), scale_(scale) {}

  QoI_Compliance_TopoOpt(const std::vector<std::vector<ROL::Ptr<FE<Real>>>> &feBdry,
                         const std::vector<std::vector<std::vector<int>>>   &bdryCellLocIds,
                         const ROL::Ptr<Traction<Real>>                     &traction,
                         const ROL::Ptr<FieldHelper<Real>>                  &fieldHelper,
                         const Real scale = 1.0)
    : feBdry_(feBdry), load_(ROL::makePtr<Load<Real>>()),
      bdryCellLocIds_(bdryCellLocIds), traction_(traction), 
      fieldHelper_(fieldHelper), scale_(scale) {}

  QoI_Compliance_TopoOpt(const ROL::Ptr<FE<Real>>                           &fe,
                         const ROL::Ptr<Load<Real>>                         &load,
                         const std::vector<std::vector<ROL::Ptr<FE<Real>>>> &feBdry,
                         const std::vector<std::vector<std::vector<int>>>   &bdryCellLocIds,
                         const ROL::Ptr<Traction<Real>>                     &traction,
                         const ROL::Ptr<FieldHelper<Real>>                  &fieldHelper,
                         const Real scale = 1.0)
    : fe_(fe), load_(load), feBdry_(feBdry), bdryCellLocIds_(bdryCellLocIds),
      traction_(traction), fieldHelper_(fieldHelper), scale_(scale) {}

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);

    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);

    // Get components of the displacement
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);

    // Evaluate on FE basis
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > valU_eval(d);
    for (int i=0; i<d; ++i) {
      valU_eval[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      fe_->evaluateValue(valU_eval[i], U[i]);
    }

    // Compute load integral
    if (!load_->isNull()) {
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > load_eval(d);
      for (int i=0; i<d; ++i) {
        load_eval[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      }
      load_->compute(load_eval, fe_, QoI<Real>::getParameter());
      for (int i=0; i<d; ++i) {
        fe_->computeIntegral(val,load_eval[i],valU_eval[i],true);
      }
    }

    // Compute traction integral
    if (!traction_->isNull()) {
      std::vector<std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > > traction;
      traction_->compute(traction, feBdry_, QoI<Real>::getParameter());
      const int numSideSets = bdryCellLocIds_.size();
      if (numSideSets > 0) {
        for (int i = 0; i < numSideSets; ++i) {
          if (traction[i].size() > 0) {
            const int numLocalSideIds = bdryCellLocIds_[i].size();
            for (int j = 0; j < numLocalSideIds; ++j) {
              const int numCellsSide  = bdryCellLocIds_[i][j].size();
              if (numCellsSide > 0 && traction[i][j].size() > 0) {
                const int numCubPerSide = feBdry_[i][j]->cubPts()->dimension(1);
                for (int k = 0; k < d; ++k) {
                  // Evaluate i-component of state on FE basis
                  ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff_bdry
                    = getBoundaryCoeff(*U[k], i, j);
                  ROL::Ptr<Intrepid::FieldContainer<Real> > valU_eval
                    = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
                  feBdry_[i][j]->evaluateValue(valU_eval, u_coeff_bdry);
                  // Integrate cell traction compliance
                  ROL::Ptr<Intrepid::FieldContainer<Real> > intVal
                    = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide);
                  feBdry_[i][j]->computeIntegral(intVal,valU_eval,traction[i][j][k],false);
                  // Add to integral value
                  for (int l = 0; l < numCellsSide; ++l) {
                    int cidx = bdryCellLocIds_[i][j][l];
                    (*val)(cidx) += (*intVal)(l);
                  }
                }
              }
            }
          }
        }
      }
    }

    Intrepid::RealSpaceTools<Real>::scale(*val, scale_);
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);

    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > G(d);
    for (int i=0; i<d; ++i) {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    }

    // Compute load integral
    if (!load_->isNull()) {
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > load_eval(d);
      for (int i=0; i<d; ++i) {
        load_eval[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      }
      load_->compute(load_eval, fe_, QoI<Real>::getParameter());
      for (int i=0; i<d; ++i) {
        Intrepid::FunctionSpaceTools::integrate<Real>(*G[i],
                                                      *load_eval[i],
                                                      *(fe_->NdetJ()),
                                                      Intrepid::COMP_CPP,
                                                      false);
      }
    }

    // Compute traction integral
    if (!traction_->isNull()) {
      traction_->apply(G, feBdry_, QoI<Real>::getParameter());
    }

    // Multiply gradient by compliance scale
    for (int i=0; i<d; ++i) {
      Intrepid::RealSpaceTools<Real>::scale(*G[i], scale_);
    }
    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Compliance_TopoOpt::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Compliance_TopoOpt::HessVec_12 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Compliance_TopoOpt::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Compliance_TopoOpt::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Compliance_TopoOpt::HessVec_22 is zero.");
  }

}; // QoI_Compliance_TopoOpt


template <class Real>
class QoI_Energy_TopoOpt : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real> > fe_;
  const ROL::Ptr<MaterialTensor<Real> > matTensor_;
  const ROL::Ptr<FieldHelper<Real> > fieldHelper_;
  const Real scale_;

public:
  QoI_Energy_TopoOpt(const ROL::Ptr<FE<Real>>             &fe,
                     const ROL::Ptr<MaterialTensor<Real>> &matTensor,
                     const ROL::Ptr<FieldHelper<Real>>    &fieldHelper,
                     const Real scale = 1.0)
    : fe_(fe), matTensor_(matTensor), fieldHelper_(fieldHelper), scale_(scale) {}

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    int matd = matTensor_->getMatrixDim();
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U, Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval, UMat, rho, rhoUMat, CUMat;
    valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    rho       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    rhoUMat   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
    CUMat     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
    fe_->evaluateValue(valZ_eval, Z[0]);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > gradU_eval(d);
    for (int i=0; i<d; ++i) {
      gradU_eval[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      fe_->evaluateGradient(gradU_eval[i], U[i]);
    }
    // Compute matrices at density rho(z)
    matTensor_->computeUmat(UMat, gradU_eval);
    matTensor_->applyTensor(CUMat, UMat);
    matTensor_->computeDensity(rho, valZ_eval);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*rhoUMat, *rho, *UMat);
    // Integrate
    fe_->computeIntegral(val, rhoUMat, CUMat, false);
    Intrepid::RealSpaceTools<Real>::scale(*val, static_cast<Real>(0.5)*scale_);
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    int matd = matTensor_->getMatrixDim();
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > G(d);
    for (int i=0; i<d; ++i) {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    }
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U, Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval, UMat, rho, rhoUMat;
    valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    rho       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    rhoUMat   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
    fe_->evaluateValue(valZ_eval, Z[0]);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > gradU_eval(d);
    for (int i=0; i<d; ++i) {
      gradU_eval[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      fe_->evaluateGradient(gradU_eval[i], U[i]);
    }
    // Compute matrices at density rho(z)
    matTensor_->computeUmat(UMat, gradU_eval);
    matTensor_->computeDensity(rho, valZ_eval);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*rhoUMat, *rho, *UMat);
    // Integrate
    for (int i=0; i<d; ++i) {
      Intrepid::FunctionSpaceTools::integrate<Real>(*G[i],
                                                    *rhoUMat,               // rho B U
                                                    *matTensor_->CBdetJ(i), // B' C
                                                    Intrepid::COMP_CPP,
                                                    false);
    }
    fieldHelper_->combineFieldCoeff(grad, G);
    Intrepid::RealSpaceTools<Real>::scale(*grad, scale_);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    int matd = matTensor_->getMatrixDim();
    // Initialize Gradients.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > G(d);
    for (int i=0; i<d; ++i) {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    }
    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U, Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval, rho, UMat, rhoUMat, CUMat, UUMat;
    valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    rho       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    rhoUMat   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
    CUMat     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
    UUMat     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valZ_eval, Z[0]);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > gradU_eval(d);
    for (int i=0; i<d; ++i) {
      gradU_eval[i] =  ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      fe_->evaluateGradient(gradU_eval[i], U[i]);
    }
    matTensor_->computeUmat(UMat, gradU_eval);
    matTensor_->applyTensor(CUMat, UMat);
    matTensor_->computeDensity(rho, valZ_eval, 1);  // first derivative
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*rhoUMat, *rho, *UMat);
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(*UUMat, *rhoUMat, *CUMat);

    /*** Evaluate weak form of the residual. ***/
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[0],
                                                  *UUMat,        // B' C drho B U
                                                  *fe_->NdetJ(), // N
                                                  Intrepid::COMP_CPP,
                                                  false);
    fieldHelper_->combineFieldCoeff(grad, G);
    Intrepid::RealSpaceTools<Real>::scale(*grad, static_cast<Real>(0.5)*scale_);
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    int matd = matTensor_->getMatrixDim();
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > H(d);
    for (int i=0; i<d; ++i) {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    }
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > V, Z;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval, VMat, rho, rhoVMat;
    valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    rho       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    rhoVMat   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
    fe_->evaluateValue(valZ_eval, Z[0]);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > gradV_eval(d);
    for (int i=0; i<d; ++i) {
      gradV_eval[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      fe_->evaluateGradient(gradV_eval[i], V[i]);
    }
    // Compute matrices at density rho(z)
    matTensor_->computeUmat(VMat, gradV_eval);
    matTensor_->computeDensity(rho, valZ_eval);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*rhoVMat, *rho, *VMat);
    // Integrate
    for (int i=0; i<d; ++i) {
      Intrepid::FunctionSpaceTools::integrate<Real>(*H[i],
                                                    *rhoVMat,               // rho B U
                                                    *matTensor_->CBdetJ(i), // B' C
                                                    Intrepid::COMP_CPP,
                                                    false);
    }

    fieldHelper_->combineFieldCoeff(hess, H);
    Intrepid::RealSpaceTools<Real>::scale(*hess, scale_);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    int matd = matTensor_->getMatrixDim();
    // Initialize Gradients.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > H(d);
    for (int i=0; i<d; ++i) {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    }
    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U, Z, V;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval, valV_eval,
      rho, rhoV, UMat, rhoUMat, CUMat, UUMat;
    valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    valV_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    rho       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    rhoV      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    rhoUMat   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
    CUMat     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
    UUMat     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valZ_eval, Z[0]);
    fe_->evaluateValue(valV_eval, V[0]);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > gradU_eval(d);
    for (int i=0; i<d; ++i) {
      gradU_eval[i] =  ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      fe_->evaluateGradient(gradU_eval[i], U[i]);
    }
    matTensor_->computeUmat(UMat, gradU_eval);
    matTensor_->applyTensor(CUMat, UMat);
    matTensor_->computeDensity(rho, valZ_eval, 1);  // first derivative
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*rhoV, *rho, *valV_eval);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*rhoUMat, *rhoV, *UMat);
    // Evaluate Hessian-times-a-vector.
    for (int i = 0; i < d; ++i) {
      Intrepid::FunctionSpaceTools::integrate<Real>(*H[i],
                                                    *rhoUMat,               // drho B U
                                                    *matTensor_->CBdetJ(i), // B' C
                                                    Intrepid::COMP_CPP,
                                                    false);
    }
    fieldHelper_->combineFieldCoeff(hess, H);
    Intrepid::RealSpaceTools<Real>::scale(*hess, scale_);
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    int matd = matTensor_->getMatrixDim();
    // Initialize Gradients.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > H(d);
    for (int i=0; i<d; ++i) {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    }
    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U, Z, V;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval,
      rho, UMat, VMat, rhoUMat, CVMat, UVMat;
    valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    rho       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    rhoUMat   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
    CVMat     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
    UVMat     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valZ_eval, Z[0]);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > gradU_eval(d);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > gradV_eval(d);
    for (int i=0; i<d; ++i) {
      gradU_eval[i] =  ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      gradV_eval[i] =  ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      fe_->evaluateGradient(gradU_eval[i], U[i]);
      fe_->evaluateGradient(gradV_eval[i], V[i]);
    }
    matTensor_->computeUmat(UMat, gradU_eval);
    matTensor_->computeUmat(VMat, gradV_eval);
    matTensor_->applyTensor(CVMat, VMat);
    matTensor_->computeDensity(rho, valZ_eval, 1);  // first derivative
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*rhoUMat, *rho, *UMat);
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(*UVMat, *rhoUMat, *CVMat);
    // Evaluate Hessian-times-a-vector.
    Intrepid::FunctionSpaceTools::integrate<Real>(*H[0],
                                                  *UVMat,        // B' C drho B U
                                                  *fe_->NdetJ(), // N
                                                  Intrepid::COMP_CPP,
                                                  false);
    fieldHelper_->combineFieldCoeff(hess, H);
    Intrepid::RealSpaceTools<Real>::scale(*hess, scale_);
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    int matd = matTensor_->getMatrixDim();
    // Initialize Hessians.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > H(d);
    for (int i=0; i<d; ++i) {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    }
    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U, Z, V;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval, valV_eval,
      rho, rhoV, UMat, rhoUMat, CUMat, UUMat;
    valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    valV_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    rho       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    rhoV      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    rhoUMat   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
    CUMat     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
    UUMat     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valZ_eval, Z[0]);
    fe_->evaluateValue(valV_eval, V[0]);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > gradU_eval(d);
    for (int i=0; i<d; ++i) {
      gradU_eval[i] =  ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      fe_->evaluateGradient(gradU_eval[i], U[i]);
    }
    matTensor_->computeUmat(UMat, gradU_eval);
    matTensor_->applyTensor(CUMat, UMat);
    matTensor_->computeDensity(rho, valZ_eval, 2);  // second derivative
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*rhoV, *rho, *valV_eval);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*rhoUMat, *rhoV, *UMat);
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(*UUMat, *rhoUMat, *CUMat);

    /*** Evaluate weak form of the residual. ***/
    Intrepid::FunctionSpaceTools::integrate<Real>(*H[0],
                                                  *UUMat,        // B' C drho B U
                                                  *fe_->NdetJ(), // N
                                                  Intrepid::COMP_CPP,
                                                  false);
    fieldHelper_->combineFieldCoeff(hess, H);
    Intrepid::RealSpaceTools<Real>::scale(*hess, static_cast<Real>(0.5)*scale_);
  }

}; // QoI_Energy_TopoOpt

#endif
