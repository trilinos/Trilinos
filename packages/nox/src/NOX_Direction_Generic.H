// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_DIRECTION_GENERIC_H
#define NOX_DIRECTION_GENERIC_H

#include "Teuchos_RCP.hpp"

// Forward declarations
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
    class LineSearchBased;
  }
}

namespace NOX {

//! Search direction strategies
/*!
  The Generic class spells out the interface. The Manager class can
  instantiate any other concrete implementation. The remaining classes
  are concrete implementations.
*/
namespace Direction {

//! %Generic direction interface
/*!
  %Generic interface for calculating a search direction, \f$d\f$, to
  be used in updating the iterate.
 */
class Generic {

public:

  //! Constructor.
  /*!
    Constructors of derived objects should look like reset().
   */
  Generic() {};

  //! Destructor
  virtual ~Generic() {};

  //! Reset direction based on possibly new parameters.
  virtual bool reset(const Teuchos::RCP<NOX::GlobalData>& gd,
             Teuchos::ParameterList& params) = 0;

  /*!
    \brief Compute the direction vector, \c dir, for a specific method given
    the current group, \c grp.

    The \c grp is not const so that we can compute the F vector, the
    Jacobian matrix, the %Newton vector, and so on.

    Const access to the solver is used for getting additional
    information such as the past solution, the iteration number, and
    so on.
  */
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
               const NOX::Solver::Generic& solver) = 0;

  /*!

    \brief Same as compute(NOX::Abstract::Vector&,
    NOX::Abstract::Group&, const NOX::Solver::Generic&)

    Enables direct support for line search based solvers for the
    purpose of efficiency since the LineSearchBased object has a
    getStep() function that some directions require.

    If it is not redefined in the derived class, it will just call the
    compute with the NOX::Solver::Generic argument.
  */
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
               const NOX::Solver::LineSearchBased& solver);
};
} // namespace Direction
} // namespace NOX

#endif
