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
#ifndef Rythmos_CONVERGENCETEST_HELPERS_H
#define Rythmos_CONVERGENCETEST_HELPERS_H

#include "Rythmos_Types.hpp"
#include "Rythmos_StepperBase.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_as.hpp"

namespace Rythmos {

using Teuchos::as;

template<class Scalar>
class LinearRegression
{
  public:
    LinearRegression(); 
    void setData(Array<Scalar>& x, Array<Scalar>& y);
    Scalar getSlope() const; 
    Scalar getYIntercept() const;
  private:
    // Private functions
    void compute_();
    void validateXYData_(Array<Scalar>& x, Array<Scalar>& y);

    // Private data
    Array<Scalar> x_;
    Array<Scalar> y_;
    Scalar slope_;
    Scalar yIntercept_;
    bool isInitialized_;
};

// Member functions:

template<class Scalar>
LinearRegression<Scalar>::LinearRegression() 
{
  isInitialized_ = false;
}

template<class Scalar>
void LinearRegression<Scalar>::setData(Array<Scalar>& x, Array<Scalar>& y)
{
  validateXYData_(x,y);
  x_ = x; // copy x data
  y_ = y; // copy y data
  isInitialized_ = true;
  compute_();
}

template<class Scalar>
void LinearRegression<Scalar>::validateXYData_(Array<Scalar>& x, Array<Scalar>& y)
{
  TEST_FOR_EXCEPT(x.size() != y.size());
  TEST_FOR_EXCEPT(x.size() < 2);
  int N = as<int>(x.size());
  // There must be at least two unique x values
  Scalar alpha = x[0];
  int numUnique = 1;
  for (int i=1; i<N ; ++i) {
    if (x[i] != alpha) {
      numUnique++;
    }
  }
  TEST_FOR_EXCEPT(numUnique==1);
}

template<class Scalar>
Scalar LinearRegression<Scalar>::getSlope() const
{
  TEST_FOR_EXCEPT(!isInitialized_);
  return slope_;
}

template<class Scalar>
Scalar LinearRegression<Scalar>::getYIntercept() const
{
  TEST_FOR_EXCEPT(!isInitialized_);
  return yIntercept_;
}

template<class Scalar>
void LinearRegression<Scalar>::compute_() 
{
  TEST_FOR_EXCEPT(!isInitialized_);
  typedef Teuchos::ScalarTraits<Scalar> ST;

  int N = Teuchos::as<int>(x_.size());

  Scalar sum1 = ST::zero();
  Scalar sum2 = ST::zero();
  for (int i=0 ; i<N ; ++i) {
    sum1 += x_[i]*y_[i];
    sum2 += x_[i]*x_[i];
  }
  sum1 *= Scalar(-2*ST::one());
  sum2 *= Scalar(-2*ST::one());

  Scalar sum3 = ST::zero();
  Scalar sum4 = ST::zero();
  for (int i=0 ; i<N ; ++i) {
    for (int j=0 ; j<N ; ++j) {
      sum3 += x_[i]*y_[j];
      sum4 += x_[i]*x_[j];
    }
  }
  sum3 *= Scalar(2*ST::one()/Scalar(N));
  sum4 *= Scalar(2*ST::one()/Scalar(N));

  slope_ = ( sum3 + sum1 ) / ( sum4 + sum2 );

  yIntercept_ = ST::zero();
  for (int i=0 ; i<N ; ++i ) {
    yIntercept_ += y_[i]-slope_*x_[i];
  }
  yIntercept_ *= Scalar(ST::one()/Scalar(N));
}

// Nonmember helper functions:
template<class Scalar>
Scalar computeLinearRegressionSlope(Array<Scalar>& x, Array<Scalar>& y) 
{
  LinearRegression<Scalar> lr;
  lr.setData(x,y);
  return(lr.getSlope());
}

template<class Scalar>
RCP<LinearRegression<Scalar> > linearRegression()
{
  RCP<LinearRegression<Scalar> > lr = rcp(new LinearRegression<Scalar>());
  return lr;
}

template<class Scalar>
class StepperFactoryBase
{
  public:
    virtual RCP<StepperBase<Scalar> > create() const =0;
};

double computeOrderByLocalErrorConvergenceStudy(const StepperFactoryBase<double>& stepperFactory);

double computeOrderByGlobalErrorConvergenceStudy(const StepperFactoryBase<double>& stepperFactory);

} // namespace Rythmos

#endif // Rythmos_CONVERGENCETEST_HELPERS_H

