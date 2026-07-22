// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PERMEABILITY_HPP
#define PERMEABILITY_HPP

#include <vector>
#include "ROL_ParameterList.hpp"
#include "ROL_Ptr.hpp"
#include "Intrepid_FieldContainer.hpp"

template<class Real>
class Permeability {
private:
  Real dRadius_, viscosity_, minPerm_, maxPerm_;
  int type_;

  Real value(const Real z) const {
    if (type_==1) {
      const Real a(minPerm_/viscosity_), b(std::log(maxPerm_/minPerm_));
      return a * std::exp(b * z);
    }
    else {
      const Real a(minPerm_/viscosity_), b(maxPerm_/viscosity_);
      return a + z * (b - a);
    }
  }

  Real deriv1(const Real z) const {
    if (type_==1) {
      const Real a(minPerm_/viscosity_), b(std::log(maxPerm_/minPerm_));
      return a * b * std::exp(b * z);
    }
    else {
      const Real a(minPerm_/viscosity_), b(maxPerm_/viscosity_);
      return b - a;
    }
  }

  Real deriv2(const Real z) const {
    if (type_==1) {
      const Real a(minPerm_/viscosity_), b(std::log(maxPerm_/minPerm_));
      return a * b * b * std::exp(b * z);
    }
    else {
      return static_cast<Real>(0);
    }
  }

public:
  Permeability(ROL::ParameterList &list) {
    bool useDarcy = list.sublist("Problem").get("Use Darcy Flow",true);
    dRadius_   = list.sublist("Problem").get("Diffuser Radius",5.0);        // mm
    viscosity_ = list.sublist("Problem").get("Dynamic Viscosity", 0.84e-8); // kg/mm-s
    minPerm_   = list.sublist("Problem").get("Minimum Permeability", 3e-7); // mm^2
    maxPerm_   = list.sublist("Problem").get("Maximum Permeability", 3e-6); // mm^2
    if (useDarcy) dRadius_ += static_cast<Real>(10);
    type_      = list.sublist("Problem").get("Parametrization Type", 0);
    // type_ = 0 ... K = z kmin + (1-z) kmax
    // type_ = 1 ... K = exp((1-z) log(kmin) + z log(kmax))
  }

  void compute(ROL::Ptr<Intrepid::FieldContainer<Real>> &alpha,
         const ROL::Ptr<Intrepid::FieldContainer<Real>> &z,
         const ROL::Ptr<const Intrepid::FieldContainer<Real>> &pts,
         const int deriv=0) const {
    const Real tol = std::sqrt(ROL::ROL_EPSILON<Real>());
    const Real zero(0), one(1);
    const int c = pts->dimension(0);
    const int p = pts->dimension(1);
    const int d = pts->dimension(2);
    Real weight(0), norm(0);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        norm = zero;
        for (int k = 0; k < d; ++k) norm += (*pts)(i,j,k)*(*pts)(i,j,k);
        weight = (std::sqrt(norm) <= dRadius_ + tol ? one : zero);
        // Compute spatially dependent viscosity
        if (deriv==1)      (*alpha)(i,j) = weight * deriv1((*z)(i,j));
        else if (deriv==2) (*alpha)(i,j) = weight * deriv2((*z)(i,j));
        else               (*alpha)(i,j) = weight * value((*z)(i,j));
      }
    }
  }
};

#endif
