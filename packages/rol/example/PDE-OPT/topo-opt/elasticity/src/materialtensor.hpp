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

#ifndef PDEOPT_MATERIALTENSOR_HPP
#define PDEOPT_MATERIALTENSOR_HPP

#include "../../../TOOLS/fe.hpp"

template <class Real>
class MaterialTensor {
private:
  ROL::Ptr<FE<Real> > fe_;

  // Problem parameters.
  Real youngsModulus_;
  Real poissonRatio_;
  bool isPlainStress_;
  Real minDensity_;
  Real maxDensity_;
  Real powerSIMP_;
  bool usePhaseField_;
  std::vector<std::vector<Real> > materialMat_;

  // Precomputed quantities.
  std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > BMat_;
  std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > BdetJMat_;
  std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > NMat_;
  std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > NdetJMat_;
  std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > CBdetJMat_;

  void computeMaterialTensor(const int d) {
    int matd = (d*(d+1))/2;
    std::vector<Real> tmpvec(matd);
    materialMat_.resize(matd, tmpvec);
    if ((d==2) && (isPlainStress_)) {
      Real one(1), half(0.5);
      Real factor1 = youngsModulus_ /(one-poissonRatio_*poissonRatio_);
      materialMat_[0][0] = factor1;
      materialMat_[0][1] = factor1 * poissonRatio_;
      materialMat_[1][0] = factor1 * poissonRatio_;
      materialMat_[1][1] = factor1;
      materialMat_[2][2] = factor1 * half * (one-poissonRatio_);
    }
    else if ((d==2) && (!isPlainStress_)) {
      Real one(1), two(2), half(0.5);
      Real factor2 = youngsModulus_ /(one+poissonRatio_)/(one-two*poissonRatio_);
      materialMat_[0][0] = factor2 * (one-poissonRatio_);
      materialMat_[0][1] = factor2 * poissonRatio_;
      materialMat_[1][0] = factor2 * poissonRatio_;
      materialMat_[1][1] = factor2 * (one-poissonRatio_);
      materialMat_[2][2] = factor2 * half * (one-two*poissonRatio_);
    }
    else {
      Real one(1), two(2), half(0.5);
      Real lam = youngsModulus_*poissonRatio_/(one+poissonRatio_)/(one-two*poissonRatio_);
      Real mu = half*youngsModulus_/(one+poissonRatio_);
      materialMat_[0][0] = lam + two*mu;
      materialMat_[0][1] = lam;
      materialMat_[0][2] = lam;
      materialMat_[1][0] = lam;
      materialMat_[1][1] = lam + two*mu;
      materialMat_[1][2] = lam;
      materialMat_[2][0] = lam;
      materialMat_[2][1] = lam;
      materialMat_[2][2] = lam + two*mu;
      materialMat_[3][3] = mu;
      materialMat_[4][4] = mu;
      materialMat_[5][5] = mu;
    }
  }

  void computeNBmats(void) {
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    int matd = materialMat_.size();

    for (int i=0; i<d; ++i) {
      BMat_.push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,p,matd));
      BdetJMat_.push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,p,matd));
      CBdetJMat_.push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,p,matd));
      NMat_.push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,p,d));
      NdetJMat_.push_back(ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,p,d));
    }

    if (d==2) {
      for (int i=0; i<c; ++i) {
        for (int j=0; j<f; ++j) {
      	  for (int k=0; k<p; ++k) {
            (*NMat_[0])(i, j, k, 0) = (*fe_->N())(i, j, k);
            (*NMat_[1])(i, j, k, 1) = (*fe_->N())(i, j, k);
            (*NdetJMat_[0])(i, j, k, 0) = (*fe_->NdetJ())(i, j, k);
            (*NdetJMat_[1])(i, j, k, 1) = (*fe_->NdetJ())(i, j, k);
                
            (*BMat_[0])(i, j, k, 0) = (*fe_->gradN())(i, j, k, 0);
            (*BMat_[1])(i, j, k, 1) = (*fe_->gradN())(i, j, k, 1);
            (*BMat_[0])(i, j, k, 2) = (*fe_->gradN())(i, j, k, 1);
            (*BMat_[1])(i, j, k, 2) = (*fe_->gradN())(i, j, k, 0);

            (*BdetJMat_[0])(i, j, k, 0) = (*fe_->gradNdetJ())(i, j, k, 0);
            (*BdetJMat_[1])(i, j, k, 1) = (*fe_->gradNdetJ())(i, j, k, 1);
            (*BdetJMat_[0])(i, j, k, 2) = (*fe_->gradNdetJ())(i, j, k, 1);
            (*BdetJMat_[1])(i, j, k, 2) = (*fe_->gradNdetJ())(i, j, k, 0);
          }
        }
      }
    }

    if(d==3) {
      for (int i=0; i<c; ++i) {
        for (int j=0; j<f; ++j) {
          for (int k=0; k<p; ++k) {
            (*NMat_[0])(i, j, k, 0) = (*fe_->N())(i, j, k);
            (*NMat_[1])(i, j, k, 1) = (*fe_->N())(i, j, k);
            (*NMat_[2])(i, j, k, 2) = (*fe_->N())(i, j, k);
            (*NdetJMat_[0])(i, j, k, 0) = (*fe_->NdetJ())(i, j, k);
            (*NdetJMat_[1])(i, j, k, 1) = (*fe_->NdetJ())(i, j, k);
            (*NdetJMat_[2])(i, j, k, 2) = (*fe_->NdetJ())(i, j, k);
            
            (*BMat_[0])(i, j, k, 0) = (*fe_->gradN())(i, j, k, 0);
            (*BMat_[1])(i, j, k, 1) = (*fe_->gradN())(i, j, k, 1);
            (*BMat_[2])(i, j, k, 2) = (*fe_->gradN())(i, j, k, 2);
            (*BMat_[1])(i, j, k, 3) = (*fe_->gradN())(i, j, k, 2);
            (*BMat_[2])(i, j, k, 3) = (*fe_->gradN())(i, j, k, 1);
            (*BMat_[0])(i, j, k, 4) = (*fe_->gradN())(i, j, k, 2);
            (*BMat_[2])(i, j, k, 4) = (*fe_->gradN())(i, j, k, 0);
            (*BMat_[0])(i, j, k, 5) = (*fe_->gradN())(i, j, k, 1);
            (*BMat_[1])(i, j, k, 5) = (*fe_->gradN())(i, j, k, 0);
                
            (*BdetJMat_[0])(i, j, k, 0) = (*fe_->gradNdetJ())(i, j, k, 0);
            (*BdetJMat_[1])(i, j, k, 1) = (*fe_->gradNdetJ())(i, j, k, 1);
            (*BdetJMat_[2])(i, j, k, 2) = (*fe_->gradNdetJ())(i, j, k, 2);
            (*BdetJMat_[1])(i, j, k, 3) = (*fe_->gradNdetJ())(i, j, k, 2);
            (*BdetJMat_[2])(i, j, k, 3) = (*fe_->gradNdetJ())(i, j, k, 1);
            (*BdetJMat_[0])(i, j, k, 4) = (*fe_->gradNdetJ())(i, j, k, 2);
            (*BdetJMat_[2])(i, j, k, 4) = (*fe_->gradNdetJ())(i, j, k, 0);
            (*BdetJMat_[0])(i, j, k, 5) = (*fe_->gradNdetJ())(i, j, k, 1);
            (*BdetJMat_[1])(i, j, k, 5) = (*fe_->gradNdetJ())(i, j, k, 0);
          }
        }
      }
    } 

    for (int i=0; i<c; ++i) {
      for (int j=0; j<f; ++j) {
        for (int k=0; k<p; ++k) {
          for (int l=0; l<d; ++l) {
            for (int m=0; m<matd; ++m) {
              for (int n=0; n<matd; ++n) {
                (*CBdetJMat_[l])(i,j,k,m) += materialMat_[m][n] * (*BdetJMat_[l])(i,j,k,n);
              }
            }
          }
        }
      }
    }

  }

public:
  MaterialTensor(Teuchos::ParameterList &parlist) {
    // Problem parameters.
    youngsModulus_  = parlist.get("Young's Modulus",     1.0);
    poissonRatio_   = parlist.get("Poisson Ratio",       0.3);
    isPlainStress_  = parlist.get("Use Plain Stress",    true);
    minDensity_     = parlist.get("Minimum Density",     1e-4);
    maxDensity_     = parlist.get("Maximum Density",     1.0);
    powerSIMP_      = parlist.get("SIMP Power",          3.0);
    usePhaseField_  = parlist.get("Use Phase Field",     false);
  }

  void setFE(ROL::Ptr<FE<Real> > &fe) {
    fe_ = fe;
    int d = fe_->gradN()->dimension(3);
    computeMaterialTensor(d);
    computeNBmats();
  }

  void applyTensor(ROL::Ptr<Intrepid::FieldContainer<Real> > & out, const ROL::Ptr<Intrepid::FieldContainer<Real> > & inData) const {
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int matd = materialMat_.size();

    for (int i=0; i<c; ++i) {
      for (int j=0; j<p; ++j) {
        for (int k=0; k<matd; ++k) {
          for (int l=0; l<matd; ++l) {
            (*out)(i,j,k) += materialMat_[k][l] * (*inData)(i,j,l);
          }
        }
      }
    }
  }

  void computeUmat(ROL::Ptr<Intrepid::FieldContainer<Real> > & UMat, std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & gradU) const {
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    int matd = materialMat_.size();

    UMat = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,matd);

    if (d==2) {
      for (int i=0; i<c; ++i) {
        for (int j=0; j<p; ++j) {
          (*UMat)(i, j, 0) = (*gradU[0])(i, j, 0);
          (*UMat)(i, j, 1) = (*gradU[1])(i, j, 1);
          (*UMat)(i, j, 2) = (*gradU[0])(i, j, 1) + (*gradU[1])(i, j, 0);
        }
      }
    }

    if(d==3) {
      for (int i=0; i<c; ++i) {
        for (int j=0; j<p; ++j) {
          (*UMat)(i, j, 0) = (*gradU[0])(i, j, 0);
          (*UMat)(i, j, 1) = (*gradU[1])(i, j, 1);
          (*UMat)(i, j, 2) = (*gradU[2])(i, j, 2);
          (*UMat)(i, j, 3) = (*gradU[1])(i, j, 2) + (*gradU[2])(i, j, 1);
          (*UMat)(i, j, 4) = (*gradU[0])(i, j, 2) + (*gradU[2])(i, j, 0);
          (*UMat)(i, j, 5) = (*gradU[0])(i, j, 1) + (*gradU[1])(i, j, 0);
        }
      }
    } 

  }

  void computeDensity(ROL::Ptr<Intrepid::FieldContainer<Real> > & rho,
                      const ROL::Ptr<Intrepid::FieldContainer<Real> > & Z,
                      const int deriv = 0) const {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);

    Real z(0);
    const Real half(0.5), one(1), diff = maxDensity_-minDensity_;
    for (int i=0; i<c; ++i) {
      for (int j=0; j<p; ++j) {
        z = (*Z)(i,j);
        if (!usePhaseField_) {
          if (deriv==0) {
            (*rho)(i,j) = minDensity_ + diff*std::pow(z, powerSIMP_);
          }
          else if (deriv==1) {
            (*rho)(i,j) = powerSIMP_*diff*std::pow(z, powerSIMP_-one);
          }
          else if (deriv==2) {
            (*rho)(i,j) = powerSIMP_*(powerSIMP_-one)*diff*std::pow(z, powerSIMP_-one);
          }
        }
        else {
          z += one;
          z *= half;
          if (deriv==0) {
            (*rho)(i,j) = minDensity_ + diff*std::pow(z, 2);
          }
          else if (deriv==1) { 
            (*rho)(i,j) = diff*z;
          }
          else if (deriv==2) {
            (*rho)(i,j) = half*diff;
          }
        }
      }
    }
  }

  const ROL::Ptr<Intrepid::FieldContainer<Real> > B(const int i) const {
    return BMat_[i];
  }

  const ROL::Ptr<Intrepid::FieldContainer<Real> > BdetJ(const int i) const {
    return BdetJMat_[i];
  }

  const ROL::Ptr<Intrepid::FieldContainer<Real> > N(const int i) const {
    return NMat_[i];
  }

  const ROL::Ptr<Intrepid::FieldContainer<Real> > NdetJ(const int i) const{
    return NdetJMat_[i];
  }

  const ROL::Ptr<Intrepid::FieldContainer<Real> > CBdetJ(const int i) const {
    return CBdetJMat_[i];
  }

  int getMatrixDim(void) const {
    return materialMat_.size();
  }
};

#endif
