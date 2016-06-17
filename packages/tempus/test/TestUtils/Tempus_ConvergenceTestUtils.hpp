#ifndef TEMPUS_CONVERGENCE_TEST_UTILS_HPP
#define TEMPUS_CONVERGENCE_TEST_UTILS_HPP

#include "Teuchos_Array.hpp"
#include "Teuchos_as.hpp"

namespace Tempus_Test {

/** \brief Linear regression class.
 *  Copied and modified from Rythmos.
 */
template<class Scalar>
class LinearRegression
{
  public:
    LinearRegression();
    void setData(Teuchos::Array<Scalar>& x, Teuchos::Array<Scalar>& y);
    Scalar getSlope() const;
    Scalar getYIntercept() const;
  private:
    // Private functions
    void compute_();
    void validateXYData_(Teuchos::Array<Scalar>& x, Teuchos::Array<Scalar>& y);

    // Private data
    Teuchos::Array<Scalar> x_;
    Teuchos::Array<Scalar> y_;
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
void LinearRegression<Scalar>::setData(
  Teuchos::Array<Scalar>& x, Teuchos::Array<Scalar>& y)
{
  validateXYData_(x,y);
  x_ = x; // copy x data
  y_ = y; // copy y data
  isInitialized_ = true;
  compute_();
}

template<class Scalar>
void LinearRegression<Scalar>::validateXYData_(
  Teuchos::Array<Scalar>& x, Teuchos::Array<Scalar>& y)
{
  TEUCHOS_TEST_FOR_EXCEPT(x.size() != y.size());
  TEUCHOS_TEST_FOR_EXCEPT(x.size() < 2);
  int N = Teuchos::as<int>(x.size());
  // There must be at least two unique x values
  Scalar alpha = x[0];
  int numUnique = 1;
  for (int i=1; i<N ; ++i) {
    if (x[i] != alpha) {
      numUnique++;
    }
  }
  TEUCHOS_TEST_FOR_EXCEPT(numUnique==1);
}

template<class Scalar>
Scalar LinearRegression<Scalar>::getSlope() const
{
  TEUCHOS_TEST_FOR_EXCEPT(!isInitialized_);
  return slope_;
}

template<class Scalar>
Scalar LinearRegression<Scalar>::getYIntercept() const
{
  TEUCHOS_TEST_FOR_EXCEPT(!isInitialized_);
  return yIntercept_;
}

template<class Scalar>
void LinearRegression<Scalar>::compute_()
{
  TEUCHOS_TEST_FOR_EXCEPT(!isInitialized_);
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
Scalar computeLinearRegression(
  Teuchos::Array<Scalar>& x, Teuchos::Array<Scalar>& y)
{
  LinearRegression<Scalar> lr;
  lr.setData(x,y);
  return(lr.getSlope());
}

template<class Scalar>
void computeLinearRegression(
  Teuchos::Array<Scalar>& x, Teuchos::Array<Scalar>& y,
  Scalar & slope, Scalar & yIntercept)
{
  LinearRegression<Scalar> lr;
  lr.setData(x,y);
  slope = lr.getSlope();
  yIntercept = lr.getYIntercept();
  return;
}

template<class Scalar>
Scalar computeLinearRegressionLogLog(
  Teuchos::Array<Scalar>& x, Teuchos::Array<Scalar>& y)
{
  TEUCHOS_TEST_FOR_EXCEPT(x.size() != y.size());
  int N = Teuchos::as<int>(x.size());
  Teuchos::Array<Scalar> xlog;
  Teuchos::Array<Scalar> ylog;

  for (int i=0 ; i<N ; ++i) {
    xlog.push_back(log(x[i]));
    ylog.push_back(log(y[i]));
  }

  LinearRegression<Scalar> lr;
  lr.setData(xlog,ylog);
  return(lr.getSlope());
}

template<class Scalar>
RCP<LinearRegression<Scalar> > linearRegression()
{
  RCP<LinearRegression<Scalar> > lr =
    Teuchos::rcp(new LinearRegression<Scalar>());
  return lr;
}

} // namespace Tempus_Test

#endif // TEMPUS_CONVERGENCE_TEST_UTILS_HPP

