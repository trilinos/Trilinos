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

#ifndef NOX_DIRECTION_MANAGER_H
#define NOX_DIRECTION_MANAGER_H

#include "NOX_Direction_Generic.H" // base class
#include "NOX_Common.H"            // class data element (string)

namespace NOX {
namespace Direction {

//! %Manager for all direction computations
/*! 
   The parameters passed to the constructor or reset should specify
   the "Method", as described below, as well as any additional
   parameters for that particular method.
 
   <ul>
   <li> "Method" - Name of the direction method. Valid choices are
   <ul> 
   <li> "%Newton" (NOX::Direction::Newton) [Default]
   <li> "Steepest Descent" (NOX::Direction::SteepestDescent)
   <li> "NonlinearCG" (NOX::Direction::NonlinearCG)
   </ul>
   </ul>
 */  
class Manager : public Generic {

public:

  //! Empty Constructor
  Manager();

  //! Constructor
  Manager(NOX::Parameter::List& params);

  //! Destructor
  ~Manager();

  // derived
  bool reset(NOX::Parameter::List& params);

  // derived
  bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp, 
	       const NOX::Solver::Generic& solver);

private:

  //! Method name
  string method;

  //! Pointer to actual direction calculation method
  NOX::Direction::Generic* ptr;

};
} // namespace Direction
} // namespace NOX

#endif
