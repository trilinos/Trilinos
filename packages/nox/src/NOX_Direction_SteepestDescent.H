// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef NOX_DIRECTION_STEEPESTDESCENT_H
#define NOX_DIRECTION_STEEPESTDESCENT_H

#include "NOX_Direction_Generic.H" // base class

// Forward declaration
namespace NOX {
class Utils;
}

namespace NOX {
namespace Direction {

//! %Calculates the steepest descent direction.
/*! 

Calculates the direction
\f[
s = - \gamma \nabla f(x) = -\gamma J(x)^T F(x)
\f]

  To use this direction, specify that the "Method" is "Steepest Descent" in the
  "Direction" sublist of the parameters that are passed to the solver
  (see NOX::Direction::Manager for more information on choosing the
  search direction).

  The following options may be specified in the "Steepest Descent"
  sublist of the "Direction" sublist of the solver parameters.

<ul>
<li>"Scaling Type" can be either of:

    <ul>
    <li>"2-Norm" - \f$ \gamma = \displaystyle\frac{1}{\|s\|_2} \f$ [default]

    <li>"Quadratic Model Min" - \f$ \gamma = \displaystyle\frac{\|s\|_2}{s^T J^T J s} \f$

    <li>"F 2-Norm" - \f$ \gamma = \displaystyle\frac{1}{\|F(x)\|_2} \f$ [default]

    <li>"None" - No scaling

    </ul>
</ul>

 */
class SteepestDescent : public Generic {

public:

  //! Scaling types
  enum ScalingType {
    //! Scale by the 2-norm of the gradient
    TwoNorm, 
    //! Scale by the step that minimizes the 1-D quadratic
    QuadMin, 
    //! Scale by the 2-norm of F
    FunctionTwoNorm, 
    //! No scaling at all
    None};

  //! Constructor
  SteepestDescent(const NOX::Utils& u, NOX::Parameter::List& params);

  //! Destructor
  ~SteepestDescent();

  // derived
  bool reset(NOX::Parameter::List& params);

  // derived
  bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp, 
	       const NOX::Solver::Generic& solver);

private:

  //! Print error message and throw error
  void throwError(const string& functionName, const string& errorMsg);

private:

  //! Printing Utils
  const NOX::Utils& utils;

  //! Temporary vector used in Cauchy scaling factor calculation
  NOX::Abstract::Vector* tmpVecPtr;

  //! Type of scaling of the Cauchy direction
  ScalingType scaleType;

};
} // namespace Direction
} // namespace NOX
#endif
