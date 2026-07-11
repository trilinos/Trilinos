
#include <iostream>
#include <cmath>

#include "Teuchos_RCP.hpp"
#include "Teuchos_as.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"

#include "Tempus_SolutionState.hpp"

/** \file Tutorial_Regression_Tester.hpp
 *
 *  \brief Regression utilities for the Tempus progression tutorials.
 *
 *  This file defines \ref tutorialRegressionTest, which performs the common
 *  final-solution regression comparison used by the van der Pol tutorial
 *  examples.
 */

/** \brief Regression check for the Tempus tutorial van der Pol examples.
 *
 *  This helper compares a computed solution against a fixed "gold" reference
 *  solution for the van der Pol tutorial problem at final time \f$t = 2.0\f$
 *  using 2000 Forward Euler steps.
 *
 *  The purpose of this check is not to validate the physical model against
 *  an analytic solution, since no closed-form analytic solution is available
 *  for the van der Pol problem. Rather, it is a regression check intended to
 *  ensure that the tutorial examples continue to reproduce the same numerical
 *  result to approximately 8 digits of relative accuracy for the fixed
 *  timestep size \f$\Delta t = 0.001\f$.
 *
 *  The current gold values are:
 *  \f[
 *    x_{\mathrm{gold}} =
 *    \left[
 *      \begin{array}{c}
 *        -1.59496108218721311 \\
 *         0.96359412806611255
 *      \end{array}
 *    \right]
 *  \f]
 *
 *  In practice, this tutorial example often reproduces the gold values much
 *  more accurately than the required 8 digits, but the acceptance threshold is
 *  intentionally looser in order to keep the regression check robust.
 *
 *  Replacing this regression test with a comparison to an analytic solution
 *  would answer a different question: it would primarily measure temporal
 *  discretization error at this timestep size. That error would be much larger
 *  than the regression tolerance used here, and would make it harder to detect
 *  simple coding mistakes that this tutorial regression check is meant to
 *  catch, such as taking one too many timesteps or changing the update logic.
 *
 *  \section tempus_tutorial_regression_update How to update the gold values
 *
 *  If the gold values ever need to be updated, use the Tempus Forward Euler
 *  verification test, `test/ForwardEuler/Tempus_ForwardEulerTest.cpp`,
 *  which performs a temporal refinement study to verify first-order accuracy.
 *  A typical update workflow is:
 *
 *  - rerun a temporal refinement study for the same problem to verify that the
 *    Forward Euler implementation still exhibits the expected first-order
 *    convergence behavior and that no first-order bugs have been introduced
 *  - once satisfied that the implementation is correct, record the final
 *    solution for the tutorial configuration:
 *    - van der Pol problem
 *    - final time \f$t = 2.0\f$
 *    - 2000 Forward Euler steps
 *  - replace the hardcoded gold values in this helper with the new reference
 *    values
 *
 *  \param[in] x_n Computed final solution vector.
 *  \param[in] out Output stream used for reporting.
 *  \param[in] relTol Relative error tolerance. The default is `1.0e-8`.
 *
 *  \return `true` if the regression comparison passes; otherwise `false`.
 */
inline bool tutorialRegressionTest(
    const Teuchos::RCP<Thyra::VectorBase<double> >& x_n,
    std::ostream& out = std::cout,
    const double relTol = Teuchos::as<double>(1.0e-8))
{
  out << std::scientific;
  auto x_regress = x_n->clone_v();
  {
    Thyra::DetachedVectorView<double> x_regress_view(*x_regress);
    x_regress_view[0] = Teuchos::as<double>(-1.59496108218721311);
    x_regress_view[1] = Teuchos::as<double>( 0.96359412806611255);
  }

  auto x_error = x_n->clone_v();
  Thyra::V_VmV(x_error.ptr(), *x_n, *x_regress);

  const double x_L2norm_error   = Thyra::norm_2(*x_error);
  const double x_L2norm_regress = Thyra::norm_2(*x_regress);
  const double relError = x_L2norm_error / x_L2norm_regress;

  out << "Relative L2 Norm of the error (regression) = "
      << relError << std::endl;

  if (x_L2norm_error > relTol * x_L2norm_regress) {
    out << "FAILED regression constraint!" << std::endl;
    return false;
  }

  return true;
}

/** \brief Overload for \ref Tempus::SolutionState.
 *
 *  Extracts the solution vector from the supplied \ref Tempus::SolutionState
 *  and applies \ref tutorialRegressionTest().
 *
 *  \param[in] solState Final solution state.
 *  \param[in] out Output stream used for reporting.
 *  \param[in] relTol Relative error tolerance. The default is `1.0e-8`.
 *
 *  \return `true` if the regression comparison passes, else `false`.
 */
inline bool tutorialRegressionTest(
    const Teuchos::RCP<Tempus::SolutionState<double> >& solState,
    std::ostream& out = std::cout,
    const double relTol = Teuchos::as<double>(1.0e-8))
{
  return tutorialRegressionTest(solState->getX(), out, relTol);
}

/** \brief Overload for the raw-array representation used in Example 0.
 *
 *  This overload converts the raw array into a \ref Thyra::VectorBase and
 *  then applies the common regression check implemented by the Thyra-vector
 *  overload of \ref tutorialRegressionTest().
 *
 *  For the vander Pol problem tutorial, N is expected to be an array of
 *  length 2.
 *
 *  \param[in] x_n Computed final solution stored as an array.
 *  \param[in] out Output stream used for reporting.
 *  \param[in] relTol Relative error tolerance. The default is `1.0e-8`.
 *
 *  \return `true` if the regression comparison passes; otherwise `false`.
 */
template<int N>
inline bool tutorialRegressionTest(
    const double (&x_n)[N],
    std::ostream& out = std::cout,
    const double relTol = Teuchos::as<double>(1.0e-8))
{
  Teuchos::RCP<const Thyra::VectorSpaceBase<double> > xSpace =
    Thyra::defaultSpmdVectorSpace<double>(N);
  Teuchos::RCP<Thyra::VectorBase<double> > x_vec = Thyra::createMember(xSpace);
  {
    Thyra::DetachedVectorView<double> x_view(*x_vec);
    for (int i = 0; i < N; ++i)
      x_view[i] = x_n[i];
  }

  return tutorialRegressionTest(x_vec, out, relTol);
}
