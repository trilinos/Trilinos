// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_LINESEARCH_NONLINEARCG_H
#define NOX_LINESEARCH_NONLINEARCG_H

#include "NOX_LineSearch_Generic.H" // base class
#include "Teuchos_RCP.hpp"

// Forward declarations
namespace NOX {
  class Utils;
  namespace StatusTest {
    class FiniteValue;
  }
}

namespace NOX {
namespace LineSearch {

//! Use NonlinearCG linesearch
/*!

  This is a simple linesearch intended to be used with
  NOX::Direction::NonlinearCG, which provides search direction \f$
  d \f$, in computing an update to the current solution vector \f$
  x_{new} = x_{old} + \lambda d \f$.  It is designed to compute a step
  length \f$ \lambda \f$ consistent with the exact linesearch of %Linear
  CG for linear problems, and it avoids use of matrices by employing a
  directional derivative (details below).  The step length, \f$ \lambda
  \f$ is computed from a single evaluation of,

  \f[
    \lambda = - \frac{F(x_{old})^T d}{d^T J(x_{old})d}
  \f]

  where \f$ J \f$ is the n x n Jacobian matrix.  Explicit construction of
  \f$ J \f$ is avoided by performing the product \f$ Jd \f$ using a
  directional derivative (cf NOX::Epetra::MatrixFree):

  \f[
  J(x_{old})d \approx \frac{F(x_{old} + \delta d) - F(x_{old})}{\delta}
  \f]

  where \f$ \delta = 10^{-6} (10^{-6} + ||x_{old}|| / ||d||) \f$ .


  <b> Derivation / Theory: </b>

  This linesearch is derived by attempting to achieve in a single step,
  the following minimization:

  \f[
    \min_\lambda \phi(\lambda)\equiv\phi (x_{old}+ \lambda d)
  \f]

  where \f$ \phi \f$ is a merit function chosen (but never explicitly
  given) so that an equivalence to %Linear CG holds, ie \f$ \nabla\phi(x)
  \leftrightarrow F(x) \f$.  The minimization above can now be cast as
  an equation:

  \f[
    \phi ' (\lambda) = \nabla\phi (x_{old}+ \lambda d)^T d =
    F(x_{old}+ \lambda d)^T d = 0~~.
  \f]

  An approximate solution to this equation can be obtained from a
  second-order expansion of \f[ \phi(\lambda) \f],

  \f[
    \phi(\lambda)\approx\phi (0) + \phi ' (0)\lambda + \phi '' (0)
    \frac{\lambda^2}{2}
  \f]

  from which it immediately follows

  \f[
    \%lambda_{min} \approx - \frac{\phi ' (0)}{\phi '' (0)} =
    - \frac{F(x_{old})^T d}{d^T J(x_{old})d}
  \f]


  <b> Input Parameters </b>

  The %NonlinearCG linesearch is selected using:

  "Line Search"

  <li> "Method" = "%NonlinearCG" [required]

  Currently, no adjustable parameters exist for this linesarch.

<b>References</b>

  <ul>

  This linesearch is adapted from ideas presented in Section 14.2 of:

  <li>Jonathan Richard Shewchuk,
  <A HREF="http://www-2.cs.cmu.edu/~jrs/jrspapers.html"/> "An
  Introduction to the Conjugate Gradient Method Without the Agonizing
  Pain</A>," 1994.</li> Though presented within the context of nonlinear
  optimization, the connection to solving nonlinear equation systems is
  made via the equivalence \f$ f'(x) \leftrightarrow F(x) \f$.

  \author Russ Hooper, Org. 9233, Sandia National Labs

*/

class NonlinearCG : public Generic {

public:

  //! Constructor
  NonlinearCG(const Teuchos::RCP<NOX::GlobalData>& gd,
          Teuchos::ParameterList& params);

  //! Destructor
  ~NonlinearCG();

  // derived
  bool reset(const Teuchos::RCP<NOX::GlobalData>& gd,
         Teuchos::ParameterList& params);

  // derived
  bool compute(NOX::Abstract::Group& newgrp, double& step,
           const NOX::Abstract::Vector& dir,
           const NOX::Solver::Generic& s);

private:

  //! Method for computing directional derivatives numerically
  NOX::Abstract::Vector& computeDirectionalDerivative(
                            const Abstract::Vector& dir,
                            const Abstract::Group& grp);

  //! Printing utilities
  Teuchos::RCP<NOX::Utils> utils;

  //! Temporary Vector pointer used to compute directional derivatives
  Teuchos::RCP<NOX::Abstract::Vector> vecPtr;

  //! Temporary Group pointer used to compute directional derivatives
  Teuchos::RCP<NOX::Abstract::Group> grpPtr;

  //! Utililty object to test for NaN or Inf values for computed steps
  Teuchos::RCP<NOX::StatusTest::FiniteValue> finiteValueTester;

};
} // namespace LineSearch
} // namespace NOX
#endif
