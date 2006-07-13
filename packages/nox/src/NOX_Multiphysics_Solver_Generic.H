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

#ifndef NOX_MULTIPHYSICS_SOLVER_GENERIC_H
#define NOX_MULTIPHYSICS_SOLVER_GENERIC_H

#include "NOX_StatusTest_Generic.H"  // base class
#include "NOX_Multiphysics_DataExchange_Interface.H"
#include "NOX_Solver_Manager.H"  // for NOX::StatusTest::StatusType

// Forward declarations
namespace NOX {
  namespace Parameter {
    class List;
  }
  namespace Abstract {
    class Group;
  }
}
#include "Teuchos_RefCountPtr.hpp"

namespace NOX {
namespace Multiphysics {

/*! 
  \brief %NOX nonlinear solvers namespace.

  The Generic class spells out the interface. The Manager class can
  instantiate any other concrete implementation. The remaining classes
  are concrete implementations.
*/
namespace Solver {

/*! 
  \brief %Abstract nonlinear solver method interface.

  Defines the type of access methods into the iterative nonlinear solvers. 

  <ul> 
  <li>Instantiate or reset() the solver.
  <li>Find the solution via solve() or perform a single iterations via iterate().
  <li>Get information about the current solver state via
      getSolutionGroup(), getPreviousSolutionGroup(),
      getNumIterations(), and getParameterList() --- particularily
      useful for NOX::StatusTest methods.
  <li>Get the current status of the solver via getStatus().
  </ul>
 */
class Generic : virtual public NOX::Solver::Generic {

public:

  //! Constructor (does nothing) 
  Generic() {};

  //! Destructor (does nothing)
  virtual ~Generic() {};

  /*! 
    \brief Reset the nonlinear solver for a new solve

    \param grp On input, contains the intial guess. On output,
    contains the solution. 

    \param tests Status tests to check for convergence or
    failure. These tests will be modified by the solver.

    \param params List of parameters. These parameters will be
    modified by the solver.

    \note All the objects passed to reset() will be modified. 

    \note The group object will be cloned via
    NOX::Abstract::Group::clone(), and the vectors within will also be
    individually cloned via NOX::Abstract::Vector::clone().

    \warning If the contents of grp, tests, or params are modified by
    the calling program after calling reset(), then the behavior of
    iterate() and solve() are completely undefined. To remedy this,
    call reset() again with the modified objects.
  */
  virtual bool 
  reset(const Teuchos::RefCountPtr<vector<NOX::Solver::Manager*> >& solvers, 
	const Teuchos::RefCountPtr<NOX::Multiphysics::DataExchange::Interface>& interface, 
	const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& tests, 
	const Teuchos::RefCountPtr<Teuchos::ParameterList>& params) = 0;

  //! reset methods inherited from NOX::Solver::Generic and needed here to avoid hiding this overloaded virtual method
  virtual bool 
  reset(const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp, 
	const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& tests, 
	const Teuchos::RefCountPtr<Teuchos::ParameterList>& params) = 0;

  virtual bool 
  reset(const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp, 
	const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& tests) = 0;

};
} // namespace Solver
} // namespace Multiphysics
} // namespace NOX


#endif
