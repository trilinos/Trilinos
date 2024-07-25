// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_DIRECTION_MODIFIEDNEWTON_H
#define NOX_DIRECTION_MODIFIEDNEWTON_H

#include "NOX_Common.H"

#ifdef WITH_PRERELEASE

#include "NOX_Direction_Generic.H" // base class
#include "Teuchos_ParameterList.hpp"       // class data element
#include "Teuchos_RCP.hpp"

// Forward Declarations
namespace NOX {
  class Utils;
  class GlobalData;
}


namespace NOX {
namespace Direction {

//! %Modified Newton direction computation
/*!
  Computes the Newton direction by solving
  \f[ Js = -f. \f]
  where \f[ J \f] may be an old Jacobian.

  To use this direction, specify that the "Method" is "Modified-Newton" in the
  "Direction" sublist of the parameters that are passed to the solver
  (see NOX::Direction::Manager for more information on choosing the
  search direction).  The default maximum age of the Jacobian is 10, and
  can be set using the "Max Age of Jacobian" parameter in the "Modified-Newton"
  sublist.
  Max Age of Jacobian  = 0 => Secant method (never update Jacobian)
                       = 1 => Newton's method (always update Jacobian)
                      >= 2 => Modified Newton's method

  \b References

  <ol>

  <li>C. T. Kelley, Iterative methods for linear and nonlinear equations, volume 16 of Frontiers in Applied Mathematics, SIAM, Philadelphia PA, 1995.

  </ol>
 */
class ModifiedNewton : public Generic {

public:

  //! Constructor
  ModifiedNewton(const Teuchos::RCP<NOX::GlobalData>& gd,
         Teuchos::ParameterList& params);

  //! Destructor
  virtual ~ModifiedNewton();

  // derived
  virtual bool reset(const Teuchos::RCP<NOX::GlobalData>& gd,
             Teuchos::ParameterList& params);

  // derived
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
               const NOX::Solver::Generic& solver);

  // derived
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
               const NOX::Solver::LineSearchBased& solver);

protected:

private:

  //! Return true if we can recover from a bad Newton solve
  /*!
    In other words, sometimes we'll use the calculated Newton
    direction even if it doesn't satisfy the convergence tolerance.

    If doRescue is true,
    \f[
    \|J s - (-F) \| / \|F\| < 1,
    \f]
    and we don't encounter any other errors, return true.
  */
  bool rescueBadNewtonSolve(const NOX::Abstract::Group& soln) const;

  //! Print an error message
  void throwError(const std::string& functionName, const std::string& errorMsg);

private:

  //! Global data pointer.  Keep this so the parameter list remains valid.
  Teuchos::RCP<NOX::GlobalData> globalDataPtr;

  //! Printing Utilities
  Teuchos::RCP<NOX::Utils> utils;

  //! "Direction" sublist with parameters for the direction vector
  /*!
    \note This is pointer rather than a reference to allow for the
    reset function.
   */
  Teuchos::ParameterList* paramsPtr;

  //! Determined based on "Rescue Bad Newton Solve". Used in rescueBadNewtonSolve
  bool doRescue;

  //! Group containing the old Jacobian
  Teuchos::RCP<NOX::Abstract::Group> oldJacobianGrpPtr;

  //! Integer counting how old the Jacobian is
  int ageOfJacobian;

  //! Integer containing the maximum age of the Jacobian
  int maxAgeOfJacobian;

};
} // namespace Direction
} // namespace NOX
#endif // NOX_DIRECTION_MODIFIEDNEWTON_H
#endif // WITH_PRERELEASE
