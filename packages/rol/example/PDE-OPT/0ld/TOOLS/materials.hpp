// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  materials.hpp
*/

#ifndef MATERIALS_HPP
#define MATERIALS_HPP

#include "Intrepid_FieldContainer.hpp"

/** \class  Materials
*/

template <class Real>
class Material {
private:

  int dim_;
  int tensorMatSize_;

  bool planeStrain_;
  Real modulus_;
  Real poissonRatio_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > materialTensor_;

public:
  
  Material() {}
  virtual ~Material() {}

  virtual void InitializeMaterial(const int dim, const bool planeStrain,
                                  const Real modulus, const Real poissonRatio) {
    TEUCHOS_TEST_FOR_EXCEPTION(dim > 3 || dim < 1, std::invalid_argument,
      ">>> ERROR (InitializeMaterial): dim less than one or greater than three!");
    dim_          = dim;
    planeStrain_  = planeStrain;
    modulus_      = modulus;
    poissonRatio_ = poissonRatio;
    if(dim_==1) {
      tensorMatSize_ = 1;
    }
    else if(dim_==2) {
      tensorMatSize_ = 3;
    }
    else {
      tensorMatSize_ = 6;
    }
	
    ComputeMaterialTensor();
  }

  void ComputeMaterialTensor(void) {
    materialTensor_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(tensorMatSize_, tensorMatSize_);
    materialTensor_ -> initialize(0.0);
    if(dim_==1) {
      (*materialTensor_)(0, 0)=modulus_;
    }
    else if(dim_==2 && !planeStrain_) {
      Real one(1), half(0.5);
      Real factor1 = modulus_ /(one-poissonRatio_*poissonRatio_);
      (*materialTensor_)(0, 0) = factor1;
      (*materialTensor_)(0, 1) = factor1 * poissonRatio_;
      (*materialTensor_)(1, 0) = factor1 * poissonRatio_;
      (*materialTensor_)(1, 1) = factor1;
      (*materialTensor_)(2, 2) = factor1 * half * (one-poissonRatio_);
    }		
    else if(dim_==2 && planeStrain_) {
      Real one(1), two(2), half(0.5);
      Real factor2 = modulus_ /(one+poissonRatio_)/(one-two*poissonRatio_);
      (*materialTensor_)(0, 0) = factor2 * (one-poissonRatio_);
      (*materialTensor_)(0, 1) = factor2 * poissonRatio_;
      (*materialTensor_)(1, 0) = factor2 * poissonRatio_;
      (*materialTensor_)(1, 1) = factor2 * (one-poissonRatio_);
      (*materialTensor_)(2, 2) = factor2 * half * (one-two*poissonRatio_);
    }
    else {
      Real one(1), two(2), half(0.5);
      Real lam = modulus_*poissonRatio_/(one+poissonRatio_)/(one-two*poissonRatio_);
      Real mu = half*modulus_/(one+poissonRatio_);
      (*materialTensor_)(0, 0) = lam + two*mu;
      (*materialTensor_)(0, 1) = lam;
      (*materialTensor_)(0, 2) = lam;
      (*materialTensor_)(1, 0) = lam;
      (*materialTensor_)(1, 1) = lam + two*mu;
      (*materialTensor_)(1, 2) = lam;
      (*materialTensor_)(2, 0) = lam;
      (*materialTensor_)(2, 1) = lam;
      (*materialTensor_)(2, 2) = lam + two*mu;
      (*materialTensor_)(3, 3) = mu;
      (*materialTensor_)(4, 4) = mu;
      (*materialTensor_)(5, 5) = mu;
    }				
  }

  const ROL::Ptr<Intrepid::FieldContainer<Real> > GetMaterialTensor(void) const {
    return materialTensor_;
  }
  
  int GetMaterialTensorDim(void) const {
    return tensorMatSize_;
  }
  
  Real GetYoungsModulus(void) const {
    return modulus_;
  }
  
  Real GetPoissonRatio(void) const {
    return poissonRatio_;
  }
  
  Real GetBulkModulus(void) const { }
  
  Real GetShearModulus(void) const {
    Real half(0.5), one(1);
    return half*modulus_/(one+poissonRatio_); 
  }
  
  void CheckMaterialTensor(void) {
    std::cout<< *materialTensor_ <<std::endl;
  }

};


// new class for topology optimization
template<class Real>
class Material_SIMP : public Material <Real> {
private:

  Real density_;
  int powerP_;
  Real minDensity_;

public:

  virtual void InitializeSIMP(const int dim, const bool planeStrain,
                              const Real modulus, const Real poissonRatio,
                              const Real density, const int powerP, const Real minDensity) {
    TEUCHOS_TEST_FOR_EXCEPTION(powerP < 1, std::invalid_argument,
      ">>> ERROR (InitializeSIMP): SIMP power is less than one!");
    Material<Real>::InitializeMaterial(dim, planeStrain, modulus, poissonRatio);
    density_    = density;
    powerP_     = powerP;
    minDensity_ = minDensity;
  }
  
  Real getSIMPScaleFactor(void) const {
    TEUCHOS_TEST_FOR_EXCEPTION(powerP_ < 1, std::invalid_argument,
      ">>> ERROR (getSIMPScaleFactor): SIMP power is less than one!");
    Real one(1);
    return minDensity_ + (one-minDensity_)*std::pow(density_, powerP_);
  }
  
  Real getSIMPFirstDerivativeScaleFactor(void) const {
    TEUCHOS_TEST_FOR_EXCEPTION(powerP_ < 1, std::invalid_argument,
      ">>> ERROR (getSIMPFirstDerivativeScaleFactor): SIMP power is less than one!");
    Real one(1);
    return (one-minDensity_) * powerP_ * (std::pow(density_, powerP_-1));
  }
  
  Real getSIMPSecondDerivativeScaleFactor(void) const {
    TEUCHOS_TEST_FOR_EXCEPTION(powerP_ < 1, std::invalid_argument,
      ">>> ERROR (getSIMPSecondDerivativeScaleFactor): SIMP power is less than one!");
    Real one(1), scale(0);
    if ( powerP_ > 1 ) {
      scale = (one-minDensity_) * powerP_ * (powerP_-1) * (std::pow(density_, powerP_-2));
    }
    return scale;
  }
  
  void setDensity(const Real dens) {
    density_ = dens;
  }
  
  // the following function is not in use
  ROL::Ptr<Intrepid::FieldContainer<Real> > computeScaledMaterialTensor(const Real scale) const {
    Real zero(0);
    int TMS = Material<Real>::getMaterialTensorDim();
    ROL::Ptr<Intrepid::FieldContainer<Real> > scaledMaterialTensor
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(TMS, TMS);
    scaledMaterialTensor -> initialize(zero);
  	
    for(int i=0; i<TMS; i++) {
      for(int j=0; j<TMS; j++) {
        (*scaledMaterialTensor)(i, j) = scale * (*Material<Real>::getMaterialTensor())(i, j);
      }
    }
    return scaledMaterialTensor;
  }

};



#endif

