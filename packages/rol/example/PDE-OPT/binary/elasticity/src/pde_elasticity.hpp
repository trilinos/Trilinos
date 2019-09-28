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

/*! \file  pde_elasticity.hpp
    \brief Implements the local PDE interface for the structural topology
           optimization problem.
*/

#ifndef PDE_TOPO_OPT_ELASTICITY_HPP
#define PDE_TOPO_OPT_ELASTICITY_HPP

#include "../../../TOOLS/pde.hpp"
#include "../../../TOOLS/fe.hpp"
#include "../../../TOOLS/fieldhelper.hpp"

#include "dirichlet.hpp"
#include "traction.hpp"
#include "load.hpp"
#include "materialtensor.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real>
class PDE_Elasticity : public PDE<Real> {
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
  // Field pattern, offsets, etc.
  std::vector<std::vector<int>> fieldPattern_;  // local Field/DOF pattern; set from DOF manager 
  int numFields_;                                // number of fields (equations in the PDE)
  int numDofs_;                                  // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                      // for each field, a counting offset
  std::vector<int> numFieldDofs_;                // for each field, number of degrees of freedom

  ROL::Ptr<Load<Real>>                        load_; 
  std::vector<ROL::Ptr<MaterialTensor<Real>>> matTensor_;
  ROL::Ptr<Dirichlet<Real>>                   dirichlet_;
  ROL::Ptr<Traction<Real>>                    traction_;
  ROL::Ptr<FieldHelper<Real>>                 fieldHelper_;
  Real minDensity_, maxDensity_;

  int M_, N_;
  Real xlbnd_, xubnd_, ylbnd_, yubnd_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> weights_;

  void computeWeight(void) {
    const int c = fe_->gradN()->dimension(0);
    const int p = fe_->gradN()->dimension(2);
    const int d = fe_->gradN()->dimension(3);

    Real x0(0), x1(0), y0(0), y1(0);
    std::vector<Real> pt(d);
    weights_.clear(); weights_.resize(M_);
    for (int i = 0; i < M_; ++i) {
      x0 = (xubnd_-xlbnd_)*static_cast<Real>(i)/static_cast<Real>(M_)+xlbnd_;
      x1 = (xubnd_-xlbnd_)*static_cast<Real>(i+1)/static_cast<Real>(M_)+xlbnd_;
      weights_[i].resize(N_);
      for (int j = 0; j < N_; ++j) {
        y0 = (yubnd_-ylbnd_)*static_cast<Real>(j)/static_cast<Real>(N_)+ylbnd_;
        y1 = (yubnd_-ylbnd_)*static_cast<Real>(j+1)/static_cast<Real>(N_)+ylbnd_;
        weights_[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
        for (int k = 0; k < c; ++k) {
          for (int l = 0; l < p; ++l) {
            for (int m = 0; m < d; ++m) {
              pt[m] = (*fe_->cubPts())(k,l,m);
            }
            // is x in ij region?
            (*weights_[i][j])(k,l) = static_cast<Real>(0);
            if ((pt[0] >= x0 && pt[0] < x1) || (i == M_-1 && pt[0] == x1)) {
              if ((pt[1] >= y0 && pt[1] < y1) || (j == N_-1 && pt[1] == y1)) {
                (*weights_[i][j])(k,l) = static_cast<Real>(1);
              }
            }
          }
        }
      }
    }
  }

public:
  PDE_Elasticity(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder    = parlist.sublist("Problem").get("Basis Order",1);
    int cubDegree     = parlist.sublist("Problem").get("Cubature Degree",4);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2);
    int probDim       = parlist.sublist("Problem").get("Problem Dimension",2);
    if (probDim> 3 || probDim < 2) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> PDE-OPT/poisson/pde_poisson.hpp: Problem dimension is not 2 or 3!");
    }
    if (basisOrder> 2 || basisOrder < 1) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> PDE-OPT/poisson/pde_poisson.hpp: Basis order is not 1 or 2!");
    }
    if (probDim == 2) {
      if (basisOrder == 1) {
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
      else if (basisOrder == 2) {
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
    }
    else if (probDim == 3) {
      if (basisOrder == 1) {
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_HEX_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
      else if (basisOrder == 2) {
        basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_HEX_C2_FEM<Real, Intrepid::FieldContainer<Real>>>();
      }
    }
    basisPtrs_.clear();
    for (int i=0; i<probDim; ++i) {
      basisPtrs_.push_back(basisPtr_);  // Displacement component
    }

    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology(); // get cell type from basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                // create cubature factory
    cellCub_ = cubFactory.create(cellType, cubDegree);                // create default cubature

    shards::CellTopology bdryCellType = cellType.getCellTopologyData(probDim-1, 0);
    bdryCub_ = cubFactory.create(bdryCellType, bdryCubDegree);

    ROL::ParameterList list;
    list.set("Use Plain Stress", parlist.sublist("Problem").get("Use Plain Stress", true));
    std::vector<Real> ym = ROL::getArrayFromStringParameter<Real>(parlist.sublist("Problem"), "Young's Modulus");
    std::vector<Real> pr = ROL::getArrayFromStringParameter<Real>(parlist.sublist("Problem"), "Poisson Ratio");
    int T = ym.size();
    matTensor_.resize(T);
    for (int t=0; t<T; ++t) {
      list.set("Young's Modulus",ym[t]);
      list.set("Poisson Ratio",pr[t]);
      matTensor_[t] = ROL::makePtr<MaterialTensor<Real>>(list);
    }
    std::string example = parlist.sublist("Problem").get("Example","Default");
    load_    = ROL::makePtr<Load<Real>>(parlist.sublist("Problem"),example);
    traction_= ROL::makePtr<Traction<Real>>(parlist.sublist("Problem"),example);
    dirichlet_ = ROL::makePtr<Dirichlet<Real>>(parlist.sublist("Problem"),example);

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

    xlbnd_ = parlist.sublist("Geometry").get("X0",0.0);
    xubnd_ = xlbnd_ + parlist.sublist("Geometry").get("Width",2.0);
    ylbnd_ = parlist.sublist("Geometry").get("Y0",0.0);
    yubnd_ = ylbnd_ + parlist.sublist("Geometry").get("Height",1.0);
    M_     = parlist.sublist("Problem").get("Number of Horizontal Cells",10);
    N_     = parlist.sublist("Problem").get("Number of Vertical Cells",20);

    minDensity_ = static_cast<Real>(1e-4)/static_cast<Real>(T);
    maxDensity_ = static_cast<Real>(1);
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real>> & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    int matd = matTensor_[0]->getMatrixDim();
    int T    = matTensor_.size(); 
 
    // Initialize residuals.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> R(d);
    for (int i=0; i<d; ++i) {
      R[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f);
    }

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);

    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> UMat, rhoCUMat, tmp;
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> CUMat(T), gradDisp_eval(d);
    rhoCUMat = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
    tmp      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
    for (int i=0; i<d; ++i) {
      gradDisp_eval[i] =  ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      fe_->evaluateGradient(gradDisp_eval[i], U[i]);
    }

    matTensor_[0]->computeUmat(UMat, gradDisp_eval);
    for (int t=0; t<T; ++t) {
      CUMat[t] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
      matTensor_[t]->applyTensor(CUMat[t],UMat);
    }
    
    Real z(0);
    for (int i=0; i<M_; ++i) {
      for (int j=0; j<N_; ++j) {
        for (int t=0; t<T; ++t) {
          z = (*z_param)[i+M_*(j+N_*t)];
          Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*tmp, *weights_[i][j], *CUMat[t]);
          for (int k=0; k<c; ++k) {
            for (int l=0; l<p; ++l) {
              for (int m=0; m<matd; ++m) {
                (*rhoCUMat)(k,l,m) += (minDensity_ + (maxDensity_-minDensity_)*z) * (*tmp)(k,l,m);
              }
            }
          }
        }
      }
    }

    for (int i=0; i<d; ++i) {
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *rhoCUMat,                // C B U
                                                    *matTensor_[0]->BdetJ(i), // B'
                                                    Intrepid::COMP_CPP,
                                                    false);
    }

    // EVALUATE LOAD
    if (!load_->isNull()) {
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> load(d);
      for (int i=0; i<d; ++i) {
        load[i] =  ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      }
      load_->compute(load, fe_, PDE<Real>::getParameter(), static_cast<Real>(-1));
      for (int i=0; i<d; ++i) {
        Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                      *load[i],           // F
                                                      *fe_->NdetJ(),      // N
                                                      Intrepid::COMP_CPP,
                                                      true);
      }
    }

    // APPLY TRACTION CONDITIONS
    if (!traction_->isNull()) {
      traction_->apply(R, feBdry_, PDE<Real>::getParameter(), static_cast<Real>(-1));
    }

    // APPLY DIRICHLET CONDITIONS
    dirichlet_->applyResidual(R,U);

    // Combine the residuals.
    fieldHelper_->combineFieldCoeff(res, R);
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
    int matd = matTensor_[0]->getMatrixDim();
    int T    = matTensor_.size();
 
    // Initialize Jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J(d);
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j) {
        J[i].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f));
      }
    }

    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> tmp;
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> rhoCBMat(d);
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> CBMat(T);
    tmp = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p, matd);

    for (int t=0; t<T; ++t) {
      CBMat[t].resize(d);
      for (int i=0; i<d; ++i) {
        CBMat[t][i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p, matd);
        matTensor_[t]->applyTensor(CBMat[t][i],matTensor_[t]->B(i));
      }
    }
    
    Real z(0);
    for (int o=0; o<d; ++o) {
      rhoCBMat[o] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p, matd);
      for (int i=0; i<M_; ++i) {
        for (int j=0; j<N_; ++j) {
          for (int t=0; t<T; ++t) {
            z = (*z_param)[i+M_*(j+N_*t)];
            Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*tmp, *weights_[i][j], *CBMat[t][o]);
            for (int k=0; k<c; ++k) {
              for (int l=0; l<f; ++l) {
                for (int m=0; m<p; ++m) {
                  for (int n=0; n<matd; ++n) {
                    (*rhoCBMat[o])(k,l,m,n) += (minDensity_ + (maxDensity_-minDensity_)*z) * (*tmp)(k,l,m,n);
                  }
                }
              }
            }
          }
        }
      }
    }

    /*** Evaluate weak form of the Jacobian. ***/
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j) {
        Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][j],
                                                      *rhoCBMat[i],            // rho B
                                                      *matTensor_[0]->BdetJ(j), // B' C
                                                      Intrepid::COMP_CPP,
                                                      false);
      }
    }

    // APPLY DIRICHLET CONDITIONS
    dirichlet_->applyJacobian1(J);

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);

  }


  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_TopoOpt::Jacobian_2): Jacobian is zero.");
  }

  void Jacobian_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    int matd = matTensor_[0]->getMatrixDim();
    int T    = matTensor_.size(); 
 
    // Initialize residuals.
    int size = z_param->size();
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J(size);
    jac.resize(size);
    for (int i=0; i<size; ++i) {
      jac[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, numDofs_);
      J[i].resize(d);
      for (int j=0; j<d; ++j) {
        J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      }
    }

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);

    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> UMat;
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> CUMat(T), rhoCUMat(size), gradDisp_eval(d);
    for (int i=0; i<d; ++i) {
      gradDisp_eval[i] =  ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      fe_->evaluateGradient(gradDisp_eval[i], U[i]);
    }

    matTensor_[0]->computeUmat(UMat, gradDisp_eval);
    for (int t=0; t<T; ++t) {
      CUMat[t] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
      matTensor_[t]->applyTensor(CUMat[t],UMat);
    }
    
    for (int i=0; i<M_; ++i) {
      for (int j=0; j<N_; ++j) {
        for (int t=0; t<T; ++t) {
          rhoCUMat[i+M_*(j+N_*t)] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
          Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*rhoCUMat[i+M_*(j+N_*t)], *weights_[i][j], *CUMat[t]);
          Intrepid::RealSpaceTools<Real>::scale(*rhoCUMat[i+M_*(j+N_*t)],maxDensity_-minDensity_);
        }
      }
    }

    for (int i=0; i<size; ++i) {
      for (int j=0; j<d; ++j) {
        Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][j],
                                                      *rhoCUMat[i],             // C B U
                                                      *matTensor_[0]->BdetJ(j), // B'
                                                      Intrepid::COMP_CPP,
                                                      false);
      }
    }

    // Combine the jacobians.
    for (int i=0; i<size; ++i) {
      // APPLY DIRICHLET CONDITIONS
      dirichlet_->applyJacobian_3(J[i]);
      fieldHelper_->combineFieldCoeff(jac[i], J[i]);
    }
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_TopoOpt::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_TopoOpt::Hessian_12): Hessian is zero.");
  }

  void Hessian_13(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    int matd = matTensor_[0]->getMatrixDim();
    int T    = matTensor_.size(); 
 
    // Initialize residuals.
    int size = z_param->size();
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> H(size);
    hess.resize(size);
    for (int i=0; i<size; ++i) {
      hess[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, numDofs_);
      H[i].resize(d);
      for (int j=0; j<d; ++j) {
        H[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f);
      }
    }

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> L;
    fieldHelper_->splitFieldCoeff(L, l_coeff);
    dirichlet_->applyMultiplier(L);

    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real>> LMat;
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> CLMat(T), rhoCLMat(size), gradDisp_eval(d);
    for (int i=0; i<d; ++i) {
      gradDisp_eval[i] =  ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      fe_->evaluateGradient(gradDisp_eval[i], L[i]);
    }

    matTensor_[0]->computeUmat(LMat, gradDisp_eval);
    for (int t=0; t<T; ++t) {
      CLMat[t] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
      matTensor_[t]->applyTensor(CLMat[t],LMat);
    }
    
    for (int i=0; i<M_; ++i) {
      for (int j=0; j<N_; ++j) {
        for (int t=0; t<T; ++t) {
          rhoCLMat[i+M_*(j+N_*t)] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, matd);
          Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*rhoCLMat[i+M_*(j+N_*t)], *weights_[i][j], *CLMat[t]);
          Intrepid::RealSpaceTools<Real>::scale(*rhoCLMat[i+M_*(j+N_*t)],maxDensity_-minDensity_);
        }
      }
    }

    for (int i=0; i<size; ++i) {
      for (int j=0; j<d; ++j) {
        Intrepid::FunctionSpaceTools::integrate<Real>(*H[i][j],
                                                      *rhoCLMat[i],             // C B U
                                                      *matTensor_[0]->BdetJ(j), // B'
                                                      Intrepid::COMP_CPP,
                                                      false);
      }
    }

    // Combine the jacobians.
    for (int i=0; i<size; ++i) {
      // APPLY DIRICHLET CONDITIONS
      dirichlet_->applyJacobian_3(H[i]);
      fieldHelper_->combineFieldCoeff(hess[i], H[i]);
    }
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_TopoOpt::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_TopoOpt::Hessian_22): Hessian is zero.");
  }

  void Hessian_23(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_TopoOpt::Hessian_23): Hessian is zero.");
  }

  void Hessian_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    Hessian_13(hess,l_coeff,u_coeff,z_coeff,z_param);
  }

  void Hessian_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_TopoOpt::Hessian_32): Hessian is zero.");
  }

  void Hessian_33(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_TopoOpt::Hessian_33): Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    //throw Exception::NotImplemented(">>> (PDE_TopoOpt::RieszMap_1): Not implemented.");

    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int d = fe_->gradN()->dimension(3);
 
    // Initialize Jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J(d);
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j) {
        J[i].push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,f));
      }
    }

    for (int i=0; i<d; ++i) {
      *(J[i][i]) = *(fe_->stiffMat());
      Intrepid::RealSpaceTools<Real>::add(*(J[i][i]),*(fe_->massMat()));
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    throw Exception::NotImplemented(">>> (PDE_TopoOpt::RieszMap_2): Not implemented.");
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
    if (!traction_->isNull()) {
      traction_->setCellNodes(bdryCellNodes_,bdryCellLocIds_);
    }
    dirichlet_->setCellNodes(bdryCellNodes_,bdryCellLocIds_,fidx_);
    int T = matTensor_.size();
    for (int t=0; t<T; ++t) {
      matTensor_[t]->setFE(fe_);
    }
    // Construct boundary FE
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets> 0) {
      feBdry_.resize(numSideSets);
      for (int i = 0; i < numSideSets; ++i) {
        int numLocSides = bdryCellNodes[i].size();
        feBdry_[i].resize(numLocSides);
        for (int j = 0; j < numLocSides; ++j) {
          if (bdryCellNodes[i][j] != ROL::nullPtr) {
            feBdry_[i][j] = ROL::makePtr<FE<Real>>(bdryCellNodes[i][j],basisPtr_,bdryCub_,j);
          }
        }
      }
    }
    computeWeight();
  }

  void setFieldPattern(const std::vector<std::vector<int>> & fieldPattern) {
    fieldPattern_ = fieldPattern;
    fieldHelper_ = ROL::makePtr<FieldHelper<Real>>(numFields_, numDofs_, numFieldDofs_, fieldPattern_);
  }

  const ROL::Ptr<FE<Real>> getFE(void) const {
    return fe_;
  }

  const std::vector<std::vector<ROL::Ptr<FE<Real>>>> getBdryFE(void) const {
    return feBdry_;
  }

  const std::vector<std::vector<std::vector<int>>> getBdryCellLocIds(void) const {
    return bdryCellLocIds_;
  }

  const ROL::Ptr<FieldHelper<Real>> getFieldHelper(void) const {
    return fieldHelper_;
  }

  const ROL::Ptr<Load<Real>> getLoad(void) const {
    return load_;
  }

  const ROL::Ptr<Traction<Real>> getTraction(void) const {
    return traction_;
  }

  const std::vector<ROL::Ptr<MaterialTensor<Real>>> getMaterialTensor(void) const {
    return matTensor_;
  }

}; // PDE_TopoOpt

#endif
