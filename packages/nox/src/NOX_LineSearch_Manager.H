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

#ifndef NOX_LINESEARCH_MANAGER_H
#define NOX_LINESEARCH_MANAGER_H

#include "NOX_LineSearch_Generic.H" // base class
#include "NOX_Common.H"		    // class data element (string)


namespace NOX {
namespace LineSearch {

//! %Manager for all different line searches
/*! 
   The parameters must specify the type of line search as well as all
   the corresponding parameters for that line search.
 
   <ul>
   <li> "Method" - Name of the line search. Valid choices are
   <ul> 
   <li> "Full Step" (NOX::LineSearch::FullStep)
   <li> "Backtrack" (NOX::LineSearch::Backtrack)
   <li> "%Polynomial" (NOX::LineSearch::Polynomial)
   <li> "Quadratic" (NOX::LineSearch::Quadratic)
   <li> "More'-Thuente" (NOX::LineSearch::MoreThuente)
   <li> "Secant" (NOX::LineSearch::Secant)
   </ul>
   </ul>
 */  
class Manager : public Generic {

public:

  //! Constructor
  Manager(Parameter::List& params);

  //! Destructor
  ~Manager();

  // derived
  bool reset(Parameter::List& params);

  // derived
  bool compute(Abstract::Group& newgrp, double& step, 
	       const Abstract::Vector& dir,
	       const Solver::Generic& s) ;

private:
  //! Method name
  string method;

  //! Pointer to method
  LineSearch::Generic* ptr;

};
} // namespace LineSearch
} // namespace NOX
#endif
