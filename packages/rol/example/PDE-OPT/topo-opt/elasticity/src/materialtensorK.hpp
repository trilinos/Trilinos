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

#ifndef PDEOPT_MATERIALTENSORK_HPP
#define PDEOPT_MATERIALTENSORK_HPP

#include "../../../TOOLS/feK.hpp"

template <class Real, class DeviceType>
class MaterialTensor {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;

private:
  ROL::Ptr<FE<Real,DeviceType>> fe_;

  // Problem parameters.
  Real youngsModulus_;
  Real poissonRatio_;
  bool isPlainStress_;
  Real minDensity_;
  Real maxDensity_;
  Real powerSIMP_;
  bool usePhaseField_;
  bool useProjection_;
  Real beta_, eta_;
  std::vector<std::vector<Real>> materialMat_;

  // Precomputed quantities.
  std::vector<scalar_view> BMat_;
  std::vector<scalar_view> BdetJMat_;
  std::vector<scalar_view> NMat_;
  std::vector<scalar_view> NdetJMat_;
  std::vector<scalar_view> CBdetJMat_;

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
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);
    int matd = materialMat_.size();

    for (int i=0; i<d; ++i) {
      BMat_.push_back(scalar_view("BMat", c,f,p,matd));
      BdetJMat_.push_back(scalar_view("BdetJMat", c,f,p,matd));
      CBdetJMat_.push_back(scalar_view("CBdetJMat", c,f,p,matd));
      NMat_.push_back(scalar_view("NMat", c,f,p,d));
      NdetJMat_.push_back(scalar_view("NdetJMat", c,f,p,d));
    }

    if (d==2) {
      for (int i=0; i<c; ++i) {
        for (int j=0; j<f; ++j) {
      	  for (int k=0; k<p; ++k) {
            (NMat_[0])(i, j, k, 0) = (fe_->N())(i, j, k);
            (NMat_[1])(i, j, k, 1) = (fe_->N())(i, j, k);
            (NdetJMat_[0])(i, j, k, 0) = (fe_->NdetJ())(i, j, k);
            (NdetJMat_[1])(i, j, k, 1) = (fe_->NdetJ())(i, j, k);
                
            (BMat_[0])(i, j, k, 0) = (fe_->gradN())(i, j, k, 0);
            (BMat_[1])(i, j, k, 1) = (fe_->gradN())(i, j, k, 1);
            (BMat_[0])(i, j, k, 2) = (fe_->gradN())(i, j, k, 1);
            (BMat_[1])(i, j, k, 2) = (fe_->gradN())(i, j, k, 0);

            (BdetJMat_[0])(i, j, k, 0) = (fe_->gradNdetJ())(i, j, k, 0);
            (BdetJMat_[1])(i, j, k, 1) = (fe_->gradNdetJ())(i, j, k, 1);
            (BdetJMat_[0])(i, j, k, 2) = (fe_->gradNdetJ())(i, j, k, 1);
            (BdetJMat_[1])(i, j, k, 2) = (fe_->gradNdetJ())(i, j, k, 0);
          }
        }
      }
    }

    if(d==3) {
      for (int i=0; i<c; ++i) {
        for (int j=0; j<f; ++j) {
          for (int k=0; k<p; ++k) {
            (NMat_[0])(i, j, k, 0) = (fe_->N())(i, j, k);
            (NMat_[1])(i, j, k, 1) = (fe_->N())(i, j, k);
            (NMat_[2])(i, j, k, 2) = (fe_->N())(i, j, k);
            (NdetJMat_[0])(i, j, k, 0) = (fe_->NdetJ())(i, j, k);
            (NdetJMat_[1])(i, j, k, 1) = (fe_->NdetJ())(i, j, k);
            (NdetJMat_[2])(i, j, k, 2) = (fe_->NdetJ())(i, j, k);
            
            (BMat_[0])(i, j, k, 0) = (fe_->gradN())(i, j, k, 0);
            (BMat_[1])(i, j, k, 1) = (fe_->gradN())(i, j, k, 1);
            (BMat_[2])(i, j, k, 2) = (fe_->gradN())(i, j, k, 2);
            (BMat_[1])(i, j, k, 3) = (fe_->gradN())(i, j, k, 2);
            (BMat_[2])(i, j, k, 3) = (fe_->gradN())(i, j, k, 1);
            (BMat_[0])(i, j, k, 4) = (fe_->gradN())(i, j, k, 2);
            (BMat_[2])(i, j, k, 4) = (fe_->gradN())(i, j, k, 0);
            (BMat_[0])(i, j, k, 5) = (fe_->gradN())(i, j, k, 1);
            (BMat_[1])(i, j, k, 5) = (fe_->gradN())(i, j, k, 0);
                
            (BdetJMat_[0])(i, j, k, 0) = (fe_->gradNdetJ())(i, j, k, 0);
            (BdetJMat_[1])(i, j, k, 1) = (fe_->gradNdetJ())(i, j, k, 1);
            (BdetJMat_[2])(i, j, k, 2) = (fe_->gradNdetJ())(i, j, k, 2);
            (BdetJMat_[1])(i, j, k, 3) = (fe_->gradNdetJ())(i, j, k, 2);
            (BdetJMat_[2])(i, j, k, 3) = (fe_->gradNdetJ())(i, j, k, 1);
            (BdetJMat_[0])(i, j, k, 4) = (fe_->gradNdetJ())(i, j, k, 2);
            (BdetJMat_[2])(i, j, k, 4) = (fe_->gradNdetJ())(i, j, k, 0);
            (BdetJMat_[0])(i, j, k, 5) = (fe_->gradNdetJ())(i, j, k, 1);
            (BdetJMat_[1])(i, j, k, 5) = (fe_->gradNdetJ())(i, j, k, 0);
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
                (CBdetJMat_[l])(i,j,k,m) += materialMat_[m][n] * (BdetJMat_[l])(i,j,k,n);
              }
            }
          }
        }
      }
    }

  }

public:
  MaterialTensor(ROL::ParameterList &parlist) {
    // Problem parameters.
    youngsModulus_  = parlist.get("Young's Modulus",      1.0);
    poissonRatio_   = parlist.get("Poisson Ratio",        0.3);
    isPlainStress_  = parlist.get("Use Plain Stress",     true);
    minDensity_     = parlist.get("Minimum Density",      1e-4);
    maxDensity_     = parlist.get("Maximum Density",      1.0);
    powerSIMP_      = parlist.get("SIMP Power",           3.0);
    usePhaseField_  = parlist.get("Use Phase Field",      false);
    useProjection_  = parlist.get("Use Projection",       false);
    beta_           = parlist.get("Projection Parameter", 1.0);
    eta_            = parlist.get("Projection Threshold", 0.5);
  }

  void setFE(ROL::Ptr<FE<Real,DeviceType>> &fe) {
    fe_ = fe;
    int d = fe_->gradN().extent_int(3);
    computeMaterialTensor(d);
    computeNBmats();
  }

  void applyTensor(scalar_view & out, const scalar_view inData) const {
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    int matd = materialMat_.size();

    for (int i=0; i<c; ++i) {
      for (int j=0; j<p; ++j) {
        for (int k=0; k<matd; ++k) {
          for (int l=0; l<matd; ++l) {
            out(i,j,k) += materialMat_[k][l] * inData(i,j,l);
          }
        }
      }
    }
  }

  void computeUmat(scalar_view & UMat, std::vector<scalar_view> & gradU) const {
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);
    int matd = materialMat_.size();

    UMat = scalar_view("Umat", c,p,matd);

    if (d==2) {
      for (int i=0; i<c; ++i) {
        for (int j=0; j<p; ++j) {
          UMat(i, j, 0) = (gradU[0])(i, j, 0);
          UMat(i, j, 1) = (gradU[1])(i, j, 1);
          UMat(i, j, 2) = (gradU[0])(i, j, 1) + (gradU[1])(i, j, 0);
        }
      }
    }

    if(d==3) {
      for (int i=0; i<c; ++i) {
        for (int j=0; j<p; ++j) {
          UMat(i, j, 0) = (gradU[0])(i, j, 0);
          UMat(i, j, 1) = (gradU[1])(i, j, 1);
          UMat(i, j, 2) = (gradU[2])(i, j, 2);
          UMat(i, j, 3) = (gradU[1])(i, j, 2) + (gradU[2])(i, j, 1);
          UMat(i, j, 4) = (gradU[0])(i, j, 2) + (gradU[2])(i, j, 0);
          UMat(i, j, 5) = (gradU[0])(i, j, 1) + (gradU[1])(i, j, 0);
        }
      }
    } 

  }

  void computeDensity(scalar_view & rho,
                      const scalar_view & Z,
                      const int deriv = 0) const {
    // Retrieve dimensions.
    int c = rho.extent_int(0);
    int p = rho.extent_int(1);

    Real z(0);
    const Real half(0.5), one(1), two(2), diff = maxDensity_-minDensity_;
    for (int i=0; i<c; ++i) {
      for (int j=0; j<p; ++j) {
        z = Z(i,j);
        if (!usePhaseField_) {
          if (!useProjection_) {
            if (deriv==0) {
              rho(i,j) = minDensity_ + diff*std::pow(z, powerSIMP_);
            }
            else if (deriv==1) {
              rho(i,j) = powerSIMP_*diff*std::pow(z, powerSIMP_-one);
            }
            else if (deriv==2) {
              rho(i,j) = powerSIMP_*(powerSIMP_-one)*diff*std::pow(z, powerSIMP_-two);
            }
          }
          else {
            if (deriv==0) {
              Real pz = project(z,0);
              rho(i,j) = minDensity_ + diff*std::pow(pz, powerSIMP_);
            }
            else if (deriv==1) {
              Real pz = project(z,0), pz1 = project(z,1);
              rho(i,j) = powerSIMP_*diff*std::pow(pz, powerSIMP_-one)*pz1;
            }
            else if (deriv==2) {
              Real pz = project(z,0), pz1 = project(z,1), pz2 = project(z,2);
              rho(i,j) = powerSIMP_*(powerSIMP_-one)*diff*std::pow(pz, powerSIMP_-two)*pz1*pz1
                         +powerSIMP_*diff*std::pow(pz, powerSIMP_-one)*pz2;
            }
          }
        }
        else {
          z += one;
          z *= half;
          if (deriv==0) {
            rho(i,j) = minDensity_ + diff*std::pow(z, 2);
          }
          else if (deriv==1) { 
            rho(i,j) = diff*z;
          }
          else if (deriv==2) {
            rho(i,j) = half*diff;
          }
        }
      }
    }
  }

  const scalar_view B(const int i) const {
    return BMat_[i];
  }

  const scalar_view BdetJ(const int i) const {
    return BdetJMat_[i];
  }

  const scalar_view N(const int i) const {
    return NMat_[i];
  }

  const scalar_view NdetJ(const int i) const{
    return NdetJMat_[i];
  }

  const scalar_view CBdetJ(const int i) const {
    return CBdetJMat_[i];
  }

  int getMatrixDim(void) const {
    return materialMat_.size();
  }

private:
  Real project(Real rho, int deriv) const {
    const Real one(1), tbe(std::tanh(beta_*eta_)), tbe1(std::tanh(beta_*(one-eta_)));
    if (deriv==1) {
      const Real tber(std::tanh(beta_*(rho-eta_)));
      return (one-tber*tber)/(tbe+tbe1)*beta_;
    }
    else if (deriv==2) {
      const Real tber(std::tanh(beta_*(rho-eta_))), two(2);
      return -two*tber*(one-tber*tber)/(tbe+tbe1)*beta_*beta_;
    }
    return (tbe+std::tanh(beta_*(rho-eta_)))/(tbe+tbe1);
  }
};

#endif
