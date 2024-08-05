// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  obj.hpp
    \brief Provides the interface for local (cell-based) objective function computations.
*/

#ifndef PDEOPT_MATERIALTENSOR_MULTIMAT_HPP
#define PDEOPT_MATERIALTENSOR_MULTIMAT_HPP

#include "../../../TOOLS/fe.hpp"

template <class Real>
class MultiMat_MaterialTensor {
private:
  ROL::Ptr<FE<Real>> fe_;

  // Problem parameters.
  Real youngsModulus_;
  Real poissonRatio_;
  bool isPlainStress_;
  Real minDensity_;
  Real maxDensity_;
  Real powerSIMP_;
  std::vector<std::vector<Real>> materialMat_;

  // Precomputed quantities.
  std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> BMat_;
  std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> BdetJMat_;
  std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> NMat_;
  std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> NdetJMat_;
  std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> CBdetJMat_;

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
  MultiMat_MaterialTensor(Teuchos::ParameterList &parlist) {
    // Problem parameters.
    youngsModulus_  = parlist.get("Young's Modulus",     1.0);
    poissonRatio_   = parlist.get("Poisson Ratio",       0.3);
    isPlainStress_  = parlist.get("Use Plain Stress",    true);
    minDensity_     = parlist.get("Minimum Density",     1e-4);
    maxDensity_     = parlist.get("Maximum Density",     1.0);
    powerSIMP_      = parlist.get("SIMP Power",          3.0);
  }

  void setFE(ROL::Ptr<FE<Real>> &fe) {
    fe_ = fe;
    int d = fe_->gradN()->dimension(3);
    computeMaterialTensor(d);
    computeNBmats();
  }

  void applyTensor(ROL::Ptr<Intrepid::FieldContainer<Real>> & out, const ROL::Ptr<Intrepid::FieldContainer<Real>> & inData) const {
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

  void computeUmat(ROL::Ptr<Intrepid::FieldContainer<Real>> & UMat, std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & gradU) const {
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

  void computeDensity(ROL::Ptr<Intrepid::FieldContainer<Real>> & rho,
                      const ROL::Ptr<Intrepid::FieldContainer<Real>> & Z,
                      const int deriv = 0) const {
    // Retrieve dimensions.
    int c = rho->dimension(0);
    int p = rho->dimension(1);

    Real z(0);
    const Real zero(0), one(1), two(2), diff = maxDensity_-minDensity_;
    for (int i=0; i<c; ++i) {
      for (int j=0; j<p; ++j) {
        z = (*Z)(i,j);
        if (deriv==0) {
          (*rho)(i,j) = minDensity_ + diff*(powerSIMP_==one ? z : std::pow(z, powerSIMP_));
        }
        else if (deriv==1) {
          (*rho)(i,j) = diff*(powerSIMP_==one ? one : powerSIMP_*std::pow(z, powerSIMP_-one));
        }
        else if (deriv==2) {
          (*rho)(i,j) = diff*(powerSIMP_==one ? zero : powerSIMP_*(powerSIMP_-one)*std::pow(z, powerSIMP_-two));
        }
      }
    }
  }

  const ROL::Ptr<Intrepid::FieldContainer<Real>> B(const int i) const {
    return BMat_[i];
  }

  const ROL::Ptr<Intrepid::FieldContainer<Real>> BdetJ(const int i) const {
    return BdetJMat_[i];
  }

  const ROL::Ptr<Intrepid::FieldContainer<Real>> N(const int i) const {
    return NMat_[i];
  }

  const ROL::Ptr<Intrepid::FieldContainer<Real>> NdetJ(const int i) const{
    return NdetJMat_[i];
  }

  const ROL::Ptr<Intrepid::FieldContainer<Real>> CBdetJ(const int i) const {
    return CBdetJMat_[i];
  }

  int getMatrixDim(void) const {
    return materialMat_.size();
  }
};

#endif
