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

#ifndef NOX_DIRECTION_NEWTON_H
#define NOX_DIRECTION_NEWTON_H

#include "NOX_Direction_Generic.H" // base class
#include "NOX_Parameter_List.H"	   // class data element 

namespace NOX {
namespace Direction {

//! %Newton direction computation
/*! 
  - "Forcing Term Method" - Method to compute the forcing term, i.e.,
    the tolerance for the linear solver. Defaults to ""
    (nothing). Choices are "Type 1" and "Type 2".

  - "Forcing Term Initial Tolerance" - Initial linear
    solver tolerance. Defaults to 1.0e-1.

  - "Forcing Term Minimum Tolerance" - Minimum acceptable linear
    solver tolerance. Defaults to 1.0e-6.

  - "Forcing Term Maximum Tolerance" = Maximum acceptable linear
    solver tolerance. Default to 0.01.

  - "Forcing Term Alpha" - Used for the "Type 2" forcing term
    calcuation. Defaults to 1.5.

  - "Forcing Term Gamma" - Used for the "Type 2" forcing term
    calcuation. Defaults to 0.9.
 */  
class Newton : public Generic {

public:

  //! Constructor
  Newton(Parameter::List& params);

  //! Destructor
  virtual ~Newton();

  // derived
  virtual bool reset(Parameter::List& params);

  // derived
  virtual bool compute(Abstract::Vector& dir, Abstract::Group& grp, 
			  const Solver::Generic& solver);

protected:

  //! Called each iteration to reset the forcing term (ie, the
  //! convergence tolerance for the linear solver).
  virtual bool resetForcingTerm(const Abstract::Group& soln, const Abstract::Group& oldsoln, int niter);


private:

  //! "Direction" sublist with parameters for the direction vector
  Parameter::List* paramsptr;

  //! Vector containing the predicted RHS
  Abstract::Vector* predf;

  //! Vector containing the computed direction
  Abstract::Vector* stepdir;


};
} // namespace Direction
} // namespace NOX
#endif
