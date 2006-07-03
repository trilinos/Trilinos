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
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef NOX_DIRECTION_MANAGER_H
#define NOX_DIRECTION_MANAGER_H

#include "NOX_Direction_Generic.H" // base class
#include "NOX_Common.H"            // class data element (string)
#include "Teuchos_RefCountPtr.hpp" 

// Forward declaration
namespace NOX {
  class Utils;
  class GlobalData;
}

namespace NOX {
namespace Direction {

//! %Manager for all direction computations
/*! 
  <b>Parameters</b>

  "Direction":
 
   <ul>
   <li> "Method" - Name of the direction method. Valid choices are
   <ul> 
   <li> "Newton" (NOX::Direction::Newton) [Default]
   <li> "Steepest Descent" (NOX::Direction::SteepestDescent)
   <li> "User Defined" - see below
   </ul>

   <li> "User Defined Constructor" - see below.

   </ul>

<b>Using a User-Defined %Direction</b>

The user has the option of passing in a user-defined direction as
follows.  

<ol>

<li> Create a direction that derives from
NOX::Direction::Generic. For example, the direction \c Foo might be
defined as shown below.

\code
class Foo : public NOX::Direction::Generic {
// Insert class definition here
}
\endcode

<li> Define a class that derives from
NOX::Parameter::DirectionConstructor. The easiest way to do this is to
use the templated version of this class, called
NOX::Parameter::DirectionConstructorT.

\code
NOX::Parameter::DirectionConstructorT<Foo> foo;
\endcode

<li> Create the appropriate entries in the parameter list, as follows.

\code
params.sublist("Direction").set("Method", "User Defined");
params.sublist("Direction").set("User Defined Constructor", foo);
\endcode

</ol>


 */  
class Manager : public Generic {

public:

  //! Empty Constructor
  Manager(const Teuchos::RefCountPtr<NOX::GlobalData>& gd);

  //! Constructor
  Manager(const Teuchos::RefCountPtr<NOX::GlobalData>& gd,
	  Teuchos::ParameterList& params);

  //! Destructor
  ~Manager();

  // derived
  bool reset(const Teuchos::RefCountPtr<NOX::GlobalData>& gd,
	     Teuchos::ParameterList& params);

  // derived
  bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp, 
	       const NOX::Solver::Generic& solver);

  // derived
  bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp, 
	       const NOX::Solver::LineSearchBased& solver);

private:

  /*!
    \brief Prints a warning message. 
    
    \param name - Name of the subroutine where the error occurs
    \param warning - The warning
  */
  void printWarning(const string& name, const string& warning);

private:

  //! Printing Utils.
  Teuchos::RefCountPtr<NOX::Utils> utils;

  //! Method name.
  string method;

  //! Pointer to actual direction calculation method.
  Teuchos::RefCountPtr<NOX::Direction::Generic> ptr;

};
} // namespace Direction
} // namespace NOX

#endif
