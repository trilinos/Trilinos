// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_LINESEARCH_UTILS_SLOPE_H
#define NOX_LINESEARCH_UTILS_SLOPE_H

#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "Teuchos_RCP.hpp"

// Forward declaration
namespace NOX {
  class GlobalData;
  namespace Abstract {
    class Vector;
    class Group;
  }
}

namespace NOX {

namespace LineSearch {

namespace Utils {

//! %Common line search utilites for computing the slope of a function.
/*!
  This class provides routines for computing the slope of a give function.  There are two methods, one that uses a Jacobian and the other that estimates the action of the Jacobian by directional derivatives.
*/

class Slope {

public:
  //! Constructor (requires call to reset() to initialize object).
  Slope();

  //! Constructor
  Slope(const Teuchos::RCP<NOX::GlobalData>& gd);

  //! Destructor
  virtual ~Slope();

  //! Reset method.
  void reset(const Teuchos::RCP<NOX::GlobalData>& gd);

  //! Compute the inner product of the given direction and the gradient associated with the given group.
  /*!
    Calculates and returns
    \f[
    \zeta = d^T \nabla f(x).
    \f]

    Here \f$d\f$ represents the input parameter \c dir and \f$\nabla
    f(x)\f$ is the gradient associated with the given group.
  */
  double computeSlope(const NOX::Abstract::Vector& dir,
              const NOX::Abstract::Group& grp);

  //! This is a variant of the computeSlope() method above optimized to work with out having to compute an explicit Jacobian.
  /*!
    Calculates and returns
    \f[
    \zeta = d^T \nabla f(x) = d^TJ^TF
    \f]

    Here \f$d\f$ represents the input parameter \c dir \f$\nabla
    f(x)\f$ is the gradient associated with the given group (for nonlinear solves this equates to \f$ J^TF \f$ where \f$ J \f$ is the Jacobian and \f$ F \f$ is the original nonlinear function).

    We can rewrite this equation as:

    \f[ d^TJ^TF = F^TJd \f]

    which allows us to use directional derivatives to estimate \f$ J^TF \f$:

    \f[ F^TJd = F^T \frac{F(x + \eta d) - F(x)}{\eta} \f]

    This may allow for faster computations of the slope if the Jacobian is expensive to evaluate.

where \f$\eta\f$ is a scalar perturbation calculated by:

\f[ \eta = \lambda * (\lambda + \frac{\| x\|}{\| d\|} ) \f]

\f$ \lambda \f$ is a constant fixed at 1.0e-6.

  */
  double computeSlopeWithOutJac(const NOX::Abstract::Vector& dir,
                const NOX::Abstract::Group& grp);

private:

  //! Printing object
  NOX::Utils utils;

  //! A vector that *may* be used in computeSlope()
  Teuchos::RCP<NOX::Abstract::Vector> vecPtr;

  //! An extra group that will only be allocated if computeSlopeWithOutJac() is called.
  Teuchos::RCP<NOX::Abstract::Group> grpPtr;

};
} // namespace Utils
} // namespace LineSearch
} // namespace NOX

#endif
