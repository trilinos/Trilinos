//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_CONVERGENCE_TEST_UTILS_HPP
#define TEMPUS_CONVERGENCE_TEST_UTILS_HPP

#include <fstream>

#include "Teuchos_as.hpp"

#include "Tempus_NumericalUtils.hpp"
#include "Tempus_String_Utilities.hpp"
#include "Tempus_Stepper.hpp"

namespace Tempus_Test {

/** \brief Linear regression class.
 *  Copied and modified from Rythmos.
 */
template <class Scalar>
class LinearRegression {
 public:
  LinearRegression();
  void setData(std::vector<Scalar>& x, std::vector<Scalar>& y);
  Scalar getSlope() const;
  Scalar getYIntercept() const;

 private:
  // Private functions
  void compute_();
  void validateXYData_(std::vector<Scalar>& x, std::vector<Scalar>& y);

  // Private data
  std::vector<Scalar> x_;
  std::vector<Scalar> y_;
  Scalar slope_;
  Scalar yIntercept_;
  bool isInitialized_;
};

// Member functions:

template <class Scalar>
LinearRegression<Scalar>::LinearRegression()
{
  isInitialized_ = false;
}

template <class Scalar>
void LinearRegression<Scalar>::setData(std::vector<Scalar>& x,
                                       std::vector<Scalar>& y)
{
  validateXYData_(x, y);
  x_             = x;  // copy x data
  y_             = y;  // copy y data
  isInitialized_ = true;
  compute_();
}

template <class Scalar>
void LinearRegression<Scalar>::validateXYData_(std::vector<Scalar>& x,
                                               std::vector<Scalar>& y)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      x.size() != y.size(), std::logic_error,
      "x and y data are note the same size for linear regression\n");
  TEUCHOS_TEST_FOR_EXCEPTION(x.size() < 2, std::logic_error,
                             "Not enough data points for linear regression!\n");
  int N = Teuchos::as<int>(x.size());
  // There must be at least two unique x values
  Scalar alpha  = x[0];
  int numUnique = 1;
  for (int i = 1; i < N; ++i) {
    if (x[i] != alpha) {
      numUnique++;
    }
  }
  TEUCHOS_TEST_FOR_EXCEPT(numUnique == 1);
}

template <class Scalar>
Scalar LinearRegression<Scalar>::getSlope() const
{
  TEUCHOS_TEST_FOR_EXCEPT(!isInitialized_);
  return slope_;
}

template <class Scalar>
Scalar LinearRegression<Scalar>::getYIntercept() const
{
  TEUCHOS_TEST_FOR_EXCEPT(!isInitialized_);
  return yIntercept_;
}

template <class Scalar>
void LinearRegression<Scalar>::compute_()
{
  TEUCHOS_TEST_FOR_EXCEPT(!isInitialized_);
  typedef Teuchos::ScalarTraits<Scalar> ST;

  int N = Teuchos::as<int>(x_.size());

  Scalar sum1 = ST::zero();
  Scalar sum2 = ST::zero();
  for (int i = 0; i < N; ++i) {
    sum1 += x_[i] * y_[i];
    sum2 += x_[i] * x_[i];
  }
  sum1 *= Scalar(-2 * ST::one());
  sum2 *= Scalar(-2 * ST::one());

  Scalar sum3 = ST::zero();
  Scalar sum4 = ST::zero();
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      sum3 += x_[i] * y_[j];
      sum4 += x_[i] * x_[j];
    }
  }
  sum3 *= Scalar(2 * ST::one() / Scalar(N));
  sum4 *= Scalar(2 * ST::one() / Scalar(N));

  slope_ = (sum3 + sum1) / (sum4 + sum2);

  yIntercept_ = ST::zero();
  for (int i = 0; i < N; ++i) {
    yIntercept_ += y_[i] - slope_ * x_[i];
  }
  yIntercept_ *= Scalar(ST::one() / Scalar(N));
}

// Nonmember helper functions:
template <class Scalar>
Scalar computeLinearRegression(std::vector<Scalar>& x, std::vector<Scalar>& y)
{
  LinearRegression<Scalar> lr;
  lr.setData(x, y);
  return (lr.getSlope());
}

template <class Scalar>
void computeLinearRegression(std::vector<Scalar>& x, std::vector<Scalar>& y,
                             Scalar& slope, Scalar& yIntercept)
{
  LinearRegression<Scalar> lr;
  lr.setData(x, y);
  slope      = lr.getSlope();
  yIntercept = lr.getYIntercept();
  return;
}

template <class Scalar>
Scalar computeLinearRegressionLogLog(std::vector<Scalar>& x,
                                     std::vector<Scalar>& y)
{
  TEUCHOS_TEST_FOR_EXCEPT(x.size() != y.size());
  int N = Teuchos::as<int>(x.size());
  std::vector<Scalar> xlog;
  std::vector<Scalar> ylog;

  for (int i = 0; i < N; ++i) {
    if (!(Tempus::approxZero(x[i]) || Tempus::approxZero(y[i]))) {
      xlog.push_back(log(x[i]));
      ylog.push_back(log(y[i]));
    }
  }

  LinearRegression<Scalar> lr;
  lr.setData(xlog, ylog);
  return (lr.getSlope());
}

template <class Scalar>
Teuchos::RCP<LinearRegression<Scalar>> linearRegression()
{
  Teuchos::RCP<LinearRegression<Scalar>> lr =
      Teuchos::rcp(new LinearRegression<Scalar>());
  return lr;
}

template <class Scalar>
void writeOrderError(
    const std::string filename, Teuchos::RCP<Tempus::Stepper<Scalar>> stepper,
    std::vector<Scalar>& StepSize,
    std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar>>>& solutions,
    std::vector<Scalar>& xErrorNorm, Scalar& xSlope,
    std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar>>>& solutionsDot,
    std::vector<Scalar>& xDotErrorNorm, Scalar& xDotSlope,
    std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar>>>& solutionsDotDot,
    std::vector<Scalar>& xDotDotErrorNorm, Scalar& xDotDotSlope,
    Teuchos::FancyOStream& out)
{
  if (solutions.empty()) {
    out << "No solutions to write out!\n"
        << std::endl;
    return;
  }
  std::size_t lastElem = solutions.size() - 1;  // Last element is the ref soln.
  std::vector<Scalar> dt;

  auto ref_solution = solutions[lastElem];
  for (std::size_t i = 0; i < lastElem; ++i) {
    dt.push_back(StepSize[i]);
    auto tmp = solutions[i];
    Thyra::Vp_StV(tmp.ptr(), -1.0, *ref_solution);
    Scalar L2norm = Thyra::norm_2(*tmp);
    xErrorNorm.push_back(L2norm);
  }
  xSlope = computeLinearRegressionLogLog<Scalar>(dt, xErrorNorm);

  if (!solutionsDot.empty()) {
    auto ref_solutionDot = solutionsDot[lastElem];
    for (std::size_t i = 0; i < lastElem; ++i) {
      auto tmp = solutionsDot[i];
      Thyra::Vp_StV(tmp.ptr(), -1.0, *ref_solutionDot);
      Scalar L2norm = Thyra::norm_2(*tmp);
      xDotErrorNorm.push_back(L2norm);
    }
    xDotSlope = computeLinearRegressionLogLog<Scalar>(dt, xDotErrorNorm);
  }

  if (!solutionsDotDot.empty()) {
    auto ref_solutionDotDot = solutionsDotDot[solutions.size() - 1];
    for (std::size_t i = 0; i < lastElem; ++i) {
      auto tmp = solutionsDotDot[i];
      Thyra::Vp_StV(tmp.ptr(), -1.0, *ref_solutionDotDot);
      Scalar L2norm = Thyra::norm_2(*tmp);
      xDotDotErrorNorm.push_back(L2norm);
    }
    xDotDotSlope = computeLinearRegressionLogLog<Scalar>(dt, xDotDotErrorNorm);
  }

  Scalar order = stepper->getOrder();
  out << "  Stepper = " << stepper->description() << std::endl;
  out << "  =================================" << std::endl;
  out << "  Expected Order = " << order << std::endl;
  out << "  x        Order = " << xSlope << std::endl;
  if (!solutionsDot.empty())
    out << "  xDot     Order = " << xDotSlope << std::endl;
  if (!solutionsDotDot.empty())
    out << "  xDotDot  Order = " << xDotDotSlope << std::endl;
  out << "  =================================" << std::endl;

  std::ofstream ftmp(filename);
  Scalar factor = 0.0;
  for (std::size_t n = 0; n < lastElem; n++) {
    factor = 0.8 * (pow(dt[n] / dt[0], order));
    ftmp << dt[n] << "   " << xErrorNorm[n] << "   " << factor * xErrorNorm[0];
    if (xDotErrorNorm.size() == lastElem)
      ftmp << "   " << xDotErrorNorm[n] << "   " << factor * xDotErrorNorm[0];
    if (xDotDotErrorNorm.size() == lastElem)
      ftmp << "   " << xDotDotErrorNorm[n] << "   "
           << factor * xDotDotErrorNorm[0];
    ftmp << std::endl;
  }
  ftmp.close();
}

template <class Scalar>
void writeOrderError(
    const std::string filename, Teuchos::RCP<Tempus::Stepper<Scalar>> stepper,
    std::vector<Scalar>& StepSize,
    std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar>>>& solutions,
    std::vector<Scalar>& xErrorNorm, Scalar& xSlope,
    std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar>>>& solutionsDot,
    std::vector<Scalar>& xDotErrorNorm, Scalar& xDotSlope,
    Teuchos::FancyOStream& out)
{
  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar>>> solutionsDotDot;
  std::vector<Scalar> xDotDotErrorNorm;
  Scalar xDotDotSlope = Scalar(0.0);
  writeOrderError(filename, stepper, StepSize, solutions, xErrorNorm, xSlope,
                  solutionsDot, xDotErrorNorm, xDotSlope, solutionsDotDot,
                  xDotDotErrorNorm, xDotDotSlope, out);
}

template <class Scalar>
void writeOrderError(
    const std::string filename, Teuchos::RCP<Tempus::Stepper<Scalar>> stepper,
    std::vector<Scalar>& StepSize,
    std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar>>>& solutions,
    std::vector<Scalar>& xErrorNorm, Scalar& xSlope, Teuchos::FancyOStream& out)
{
  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar>>> solutionsDot;
  std::vector<Scalar> xDotErrorNorm;
  Scalar xDotSlope = Scalar(0.0);
  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar>>> solutionsDotDot;
  std::vector<Scalar> xDotDotErrorNorm;
  Scalar xDotDotSlope = Scalar(0.0);
  writeOrderError(filename, stepper, StepSize, solutions, xErrorNorm, xSlope,
                  solutionsDot, xDotErrorNorm, xDotSlope, solutionsDotDot,
                  xDotDotErrorNorm, xDotDotSlope, out);
}

template <class Scalar>
void writeSolution(
    const std::string filename,
    Teuchos::RCP<const Tempus::SolutionHistory<Scalar>> solutionHistory)
{
  std::ofstream ftmp(filename);
  Teuchos::RCP<const Thyra::VectorBase<Scalar>> x;
  Teuchos::RCP<const Thyra::VectorBase<Scalar>> xDot;
  Teuchos::RCP<const Thyra::VectorBase<Scalar>> xDotDot;
  for (int i = 0; i < solutionHistory->getNumStates(); i++) {
    Teuchos::RCP<const Tempus::SolutionState<Scalar>> solutionState =
        (*solutionHistory)[i];
    Scalar time = solutionState->getTime();
    x           = solutionState->getX();
    xDot        = solutionState->getXDot();
    xDotDot     = solutionState->getXDotDot();
    int J       = x->space()->dim();

    ftmp << time;
    for (int j = 0; j < J; j++) ftmp << "   " << get_ele(*x, j);
    if (xDot != Teuchos::null)
      for (int j = 0; j < J; j++) ftmp << "   " << get_ele(*xDot, j);
    if (xDotDot != Teuchos::null)
      for (int j = 0; j < J; j++) ftmp << "   " << get_ele(*xDotDot, j);
    ftmp << std::endl;
  }
  ftmp.close();
}

template <class Scalar>
void writeSolution(
    const std::string filename,
    Teuchos::RCP<Tempus::SolutionHistory<Scalar>> solutionHistory)
{
  Teuchos::RCP<const Tempus::SolutionHistory<Scalar>> sh(solutionHistory);
  writeSolution(filename, sh);
}

}  // namespace Tempus_Test

#endif  // TEMPUS_CONVERGENCE_TEST_UTILS_HPP
