// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.hpp
    \brief Contains Rosenbrock's objective function definition for ArrayFire arrays.
    \addtogroup examples_group
*/


// Whether or not to use the exact Hessian-times-a-vector
#ifndef USE_HESSVEC 
#define USE_HESSVEC 0
#endif

#ifndef ROL_ROSENBROCK_AF_HPP
#define ROL_ROSENBROCK_AF_HPP

#include "ROL_Objective.hpp"
#include "ROL_ArrayFireVector.hpp"

/** \brief Rosenbrock's function with ArrayFire data structures.
*/
template< class Real, class Element = Real>
class Objective_Rosenbrock_AF : public ROL::Objective<Real> {

private:
  Element alpha_;

public:

  Objective_Rosenbrock_AF(Element alpha = 100.0) : alpha_(alpha) {}


  Real value(const ROL::Vector<Real> &x, Real &tol) {

    const ROL::ArrayFireVector<Real, Element> &ex = dynamic_cast<const ROL::ArrayFireVector<Real, Element>&>(x);
    ROL::Ptr<const af::array> afVector = ex.getVector();

    Element one(1);

    af::array val = af::sum(alpha_ * af::pow(af::pow((*afVector)(af::seq(0, af::end, 2)), 2) - (*afVector)(af::seq(1, af::end, 2)), 2) +
                            af::pow((*afVector)(af::seq(0, af::end, 2)) - one, 2), -1);

    return(static_cast<Real>(val.scalar<Element>()));
  }	


  void gradient(ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real &tol) {

    const ROL::ArrayFireVector<Real, Element> &ex = dynamic_cast<const ROL::ArrayFireVector<Real, Element>&>(x);
    ROL::ArrayFireVector<Real, Element> &eg = dynamic_cast<ROL::ArrayFireVector<Real, Element>&>(g);
    ROL::Ptr<const af::array> afx = ex.getVector();
    ROL::Ptr<af::array> afg = eg.getVector();

    Element one(1), two(2), four(4);

    // Vectorized implementation.
    (*afg)(af::seq(0, af::end, 2), 2) =  four * alpha_ * (*afx)(af::seq(0, af::end, 2)) *
                                                         ( af::pow((*afx)(af::seq(0, af::end, 2)), 2) - (*afx)(af::seq(1, af::end, 2)) ) +
                                         two * ( (*afx)(af::seq(0, af::end, 2)) - one );
    (*afg)(af::seq(1, af::end, 2), 2) = -two * alpha_ *  ( af::pow((*afx)(af::seq(0, af::end, 2)), 2) - (*afx)(af::seq(1, af::end, 2)) );

    /*
    // Componentwise implementation.
    dim_t n = afx->dims(0);
    for (dim_t i = 0; i<n/2; i++) {
      (*afg)(2 * i) = four*alpha_*(af::pow((*afx)(2 * i), 2) - (*afx)(2 * i + 1))*(*afx)(2 * i) + two*((*afx)(2 * i) - one);
      (*afg)(2 * i + 1) = -two*alpha_*(af::pow((*afx)(2 * i), 2) - (*afx)(2 * i + 1));
    }
    */
  }

#if USE_HESSVEC
  /* needs reimplementing
		void hessVec(Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {

			
			ROL::Ptr<const vector> xp = getVector<XPrim>(x);
			ROL::Ptr<const vector> vp = getVector<XPrim>(v);
			ROL::Ptr<vector> hvp = getVector<XDual>(hv);

			uint n = xp->size();
			for (uint i = 0; i<n / 2; i++) {
				Real h11 = 4.0*alpha_*(3.0*pow((*xp)[2 * i], 2) - (*xp)[2 * i + 1]) + 2.0;
				Real h12 = -4.0*alpha_*(*xp)[2 * i];
				Real h22 = 2.0*alpha_;

				(*hvp)[2 * i] = h11*(*vp)[2 * i] + h12*(*vp)[2 * i + 1];
				(*hvp)[2 * i + 1] = h12*(*vp)[2 * i] + h22*(*vp)[2 * i + 1];
			}
		}
  */
#endif

};


#endif //ROL_ROSENBROCK_AF_HPP
