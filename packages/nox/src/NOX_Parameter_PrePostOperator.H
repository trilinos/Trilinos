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

#ifndef NOX_PARAMETER_PREPOSTOPERATOR_H
#define NOX_PARAMETER_PREPOSTOPERATOR_H

#include "NOX_Common.H"  // for NOX_Config.h
#include "NOX_Parameter_Arbitrary.H" // base class

// Forward Declarations
namespace NOX {
  namespace Solver {
    class Generic;
  }
}

namespace NOX {
namespace Parameter {

/*! 

  \brief %NOX's pure virtual class to allow users to insert pre and
  post operations into nox's solvers (before and after the
  NOX::Solver::Generic::iterate() and NOX::Solver::Generic::solve()
  methods).

  This class is a member of the namespace NOX::Parameter and is
  derived from the NOX:Parameter::Arbitrary class.

  The user should implement their own concrete implementation of this
  class and register it in the "Solver Options" sublist.

  
  To create and use a user defined pre/post operators:

  <ol>

  <li> Create a pre/post operator that derives from
  NOX::Parameter::PrePostOperator. For example, the pre/post oeprator \c
  Foo might be defined as shown below.

  \code
  class Foo : public NOX::Parameter::PrePostOperator {
  // Insert class definition here
  }
  \endcode

  <li> Create the appropriate entries in the parameter list, as follows.

  \code
  Foo foo();
  params.sublist("Line Search Based Sovler").setParameter("User Defined Pre/Post Operator", foo);
  \endcode

  </ol>

  \author Tammy Kolda (SNL 8950), Roger Pawlowski (SNL 9233)
*/

class PrePostOperator : public NOX::Parameter::Arbitrary {

public:
  
  //! %Abstract %Vector constructor (does nothing)
  PrePostOperator();

  //! Copy constructor (doesnothing)
  PrePostOperator(const NOX::Parameter::PrePostOperator& source);

  //! %Abstract %Vector destructor (does nothing)
  virtual ~PrePostOperator();

  //! User defined method that will be executed at the start of a call to NOX::Solver::Generic::iterate().
  virtual void runPreIterate(const NOX::Solver::Generic& solver);

  //! User defined method that will be executed at the end of a call to NOX::Solver::Generic::iterate().
  virtual void runPostIterate(const NOX::Solver::Generic& solver);

  //! User defined method that will be executed at the start of a call to NOX::Solver::Generic::solve().
  virtual void runPreSolve(const NOX::Solver::Generic& solver);

  //! User defined method that will be executed at the end of a call to NOX::Solver::Generic::solve().
  virtual void runPostSolve(const NOX::Solver::Generic& solver);

}; // class PrePostOperator
} // namespace Parameter
} // namespace NOX

#endif
