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

#ifndef NOX_LINESEARCH_MANAGER_H
#define NOX_LINESEARCH_MANAGER_H

#include "NOX_LineSearch_Generic.H"     // base class
#include "NOX_Common.H"		        // class data element (string)
#include "Teuchos_RefCountPtr.hpp"      // class data element

// Forward declarations
namespace NOX {
  class Utils;
}

namespace NOX {
namespace LineSearch {

//! %Manager for all different line searches
/*! 
  <b>Parameters</b>
  
  "LineSearch":
   <ul>
   <li> "Method" - Name of the line search. Valid choices are
   <ul> 
   <li> "Full Step" (NOX::LineSearch::FullStep)
   <li> "Backtrack" (NOX::LineSearch::Backtrack)
   <li> "Polynomial" (NOX::LineSearch::Polynomial)
   <li> "More'-Thuente" (NOX::LineSearch::MoreThuente)
   <li> "User Defined" - see below

   </ul>

   <li> "User Defined Constructor" - see below

   </ul>

<b>Using a User-Defined Line Search</b>

The user has the option of passing in a user-defined line search as
follows.  

<ol>

<li> Create a line search that derives from
NOX::LineSearch::Generic. For example, the direction \c Foo might be
defined as shown below.

\code
class Foo : public NOX::LineSearch::Generic {
// Insert class definition here
}
\endcode

<li> Define a class that derives from
NOX::Parameter::LineSearchConstructor. The easiest way to do this is to
use the templated version of this class, called
NOX::Parameter::LineSearchConstructorT.

\code
NOX::Parameter::LineSearchConstructorT<Foo> foo;
\endcode

<li> Create the appropriate entries in the parameter list, as follows.

\code
params.sublist("Line Search").set("Method", "User Defined");
params.sublist("Line Search").set("User Defined Constructor", foo);
\endcode

</ol>

 */  
class Manager : public Generic {

public:

  //! Empty constructor, reset must be called if this ctor is used.
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
  bool compute(NOX::Abstract::Group& newgrp, double& step, 
	       const NOX::Abstract::Vector& dir,
	       const NOX::Solver::Generic& s) ;

private:

  /*!
    \brief Prints a warning message. 
    
    \param name - Name of the subroutine where the error occurs
    \param warning - The warning message 
  */
  void printWarning(const string& name, const string& warning);

private:
  //! Printing utilities
  Teuchos::RefCountPtr<NOX::Utils> utils;

  //! Method name
  string method;

  //! Pointer to method
  Teuchos::RefCountPtr<NOX::LineSearch::Generic> ptr;

};
} // namespace LineSearch
} // namespace NOX
#endif
