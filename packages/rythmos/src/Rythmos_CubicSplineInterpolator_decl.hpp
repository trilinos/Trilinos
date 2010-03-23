//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_CUBIC_SPLINE_INTERPOLATOR_DECL_H
#define Rythmos_CUBIC_SPLINE_INTERPOLATOR_DECL_H

#include "Rythmos_InterpolatorBase.hpp"
#include "Rythmos_Types.hpp"


namespace Rythmos {

using Teuchos::Ptr;

template<class Scalar>
class CubicSplineCoeff {
  public:
    Array<Scalar> t;
    Array<RCP<Thyra::VectorBase<Scalar> > > a;
    Array<RCP<Thyra::VectorBase<Scalar> > > b;
    Array<RCP<Thyra::VectorBase<Scalar> > > c;
    Array<RCP<Thyra::VectorBase<Scalar> > > d;
};


/** \brief Concrete implemenation of <tt>InterpolatorBase</tt> that implements
 * cubic spline interpolation
 */
template<class Scalar>
class CubicSplineInterpolator : virtual public InterpolatorBase<Scalar>
{
public:

  /** \brief . */
  ~CubicSplineInterpolator() {};

  /** \brief . */
  CubicSplineInterpolator();

  /** \brief . */
  bool supportsCloning() const;

  /** \brief . */
  RCP<InterpolatorBase<Scalar> > cloneInterpolator() const;

  /** \brief . */
  void setNodes(
    const RCP<const typename DataStore<Scalar>::DataStoreVector_t> & nodes
    );

  /** \brief . */
  void interpolate(
    const Array<Scalar> &t_values,
    typename DataStore<Scalar>::DataStoreVector_t *data_out
    ) const;

  /** \brief . */
  int order() const; 

  /** \brief . */
  std::string description() const;

  /** \brief . */
  void describe(
    FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);

  /** \brief . */
  RCP<ParameterList> getNonconstParameterList();

  /** \brief . */
  RCP<ParameterList> unsetParameterList();

  /** \brief. */
  RCP<const Teuchos::ParameterList> getValidParameters() const;

private:


  RCP<const typename DataStore<Scalar>::DataStoreVector_t> nodes_;
#ifdef RYTHMOS_DEBUG
  RCP<typename DataStore<Scalar>::DataStoreVector_t> nodes_copy_;
#endif // RYTHMOS_DEBUG

  mutable CubicSplineCoeff<Scalar> splineCoeff_;
  mutable bool splineCoeffComputed_;
  bool nodesSet_;

  RCP<ParameterList> parameterList_;

};

// non-member constructor
template<class Scalar>
RCP<CubicSplineInterpolator<Scalar> > cubicSplineInterpolator();

// Non-member helper function
// Algorithm from "CRC Standard Mathematical Tables and Formulae" by Daniel
// Zwillinger, 31st Edition, pg 738, natural cubic splines
// input:  n, {x_0, x_1, ..., x_n}
//         a_0 = f(x_0), a_1 = f(x_1), ... a_n = f(x_n)
// output: {a_j,b_j,c_j,d_j}, j=0..n-1
template<class Scalar>
void computeCubicSplineCoeff(
    const typename DataStore<Scalar>::DataStoreVector_t & data,
    const Ptr<CubicSplineCoeff<Scalar> > & coeffPtr
    );

template<class Scalar>
void validateCubicSplineCoeff(const CubicSplineCoeff<Scalar>& coeff);

// Non-member helper function
// Evaluate cubic spline 
template<class Scalar>
void evaluateCubicSpline(
    const CubicSplineCoeff<Scalar>& coeff,
    Teuchos::Ordinal j, 
    const Scalar& t,
    const Ptr<Thyra::VectorBase<Scalar> >& S,
    const Ptr<Thyra::VectorBase<Scalar> >& Sp = Teuchos::null, 
    const Ptr<Thyra::VectorBase<Scalar> >& Spp = Teuchos::null
    );


} // namespace Rythmos


#endif // Rythmos_CUBIC_SPLINE_INTERPOLATOR_DECL_H
