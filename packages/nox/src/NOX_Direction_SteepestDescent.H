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

namespace NOX {
namespace Direction {

//! %Calculates the steepest descent direction.
/*! 

Calculates the direction
\f[
d = - \gamma \nabla f(x) = -\gamma J(x)^T F(x)
\f]
where we are solving 
\f$ F(x) = 0 \f$ 
and 
\f$ f(x) = \frac{1}{2} \|F(x)\|_2 \f$ and \f$ J(x) = \nabla F(x) \f$.
The following parameters are valid for this solver:

<ul>
<li>"Scaling Type" can be either of:
    <ul>
    <li>"2-Norm" - \f$ \gamma = \displaystyle\frac{1}{\|d\|_2} \f$
    <li>"Quadratic Model Min" - \f$ \gamma = \displaystyle\frac{\|d\|_2}{d^T J^T J d} \f$
    <li>"None" - No scaling
    </ul>
</ul>

 */
class SteepestDescent : public Generic {

public:

  //! Scaling types
  enum ScalingType {TwoNorm, QuadMin, None};

  //! Constructor
  SteepestDescent(Parameter::List& params);

  //! Destructor
  ~SteepestDescent();

  // derived
  bool reset(Parameter::List& params);

  // derived
  bool compute(Abstract::Vector& dir, Abstract::Group& grp, 
		  const Solver::Generic& solver);

private:

  //! Temporary vector used in Cauchy scaling factor calculation
  NOX::Abstract::Vector* tmpVecPtr;

  //! Type of scaling of the Cauchy direction
  ScalingType scaleType;

};
} // namespace Direction
} // namespace NOX
#endif
