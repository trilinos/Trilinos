// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_LINESEARCH_GENERIC_H
#define NOX_LINESEARCH_GENERIC_H

#include "Teuchos_RCP.hpp"

// Forward Declarations
namespace Teuchos {
  class ParameterList;
}
namespace NOX {
  class GlobalData;
  namespace Abstract {
    class Vector;
    class Group;
  }
  namespace Solver {
    class Generic;
  }
}

namespace NOX {

//! NOX linesearches
namespace LineSearch {

/*!  \brief Base class line search interface

  Every line search should respect the following Parameter:

  - "Max Iters" - maximum number of iterations (i.e., RHS computations)
*/
class Generic {

public:

  //! Default constructor
  Generic() {};

  //! Destructor
  virtual ~Generic() {};

  //! Perform a line search
  /*!

  On input:
  @param grp The initial solution vector, \f$x_{\rm old}\f$.
  @param dir A vector of directions to be used in the line search, \f$d\f$.
  @param s The nonlinear solver.

  On output:
  @param step The distance the direction was scaled, \f$ \lambda \f$.
  @param grp The \c grp is updated with a new solution, \f$ x_{\rm new} \f$, resulting from the linesearch.  Normally, for a single direction line search, this is computed as:

    \f[
    x_{\rm new} = x_{\rm old} + \lambda d.
    \f]

  Ideally, \f$ \|F(x_{\rm new})\| < \|F(x_{\rm old})\| \f$ (e.g the final direction is a descent direction).

  Note that the \c dir object is a std::vector. For typical line
  searches as described in the above equation, this vector is of size
  one.  We have used a std::vector to allow for special cases of
  multi-directional line searches such as the Bader/Schnabel
  curvillinear line search.

  Return value is true for a successful line search computation.

   */
  virtual bool compute(NOX::Abstract::Group& grp,
               double& step,
               const NOX::Abstract::Vector& dir,
               const NOX::Solver::Generic& s) = 0;

};
} // namespace LineSearch
} // namespace NOX

#endif
