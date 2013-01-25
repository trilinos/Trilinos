//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_MULTIPHYSICS_SOLVER_GENERIC_H
#define NOX_MULTIPHYSICS_SOLVER_GENERIC_H

#include "NOX_StatusTest_Generic.H"  // base class
#include "NOX_Multiphysics_DataExchange_Interface.H"
#include "NOX_Solver_Generic.H"  

// Forward declarations
namespace NOX {
  namespace Parameter {
    class List;
  }
  namespace Abstract {
    class Group;
  }
}
#include "Teuchos_RCP.hpp"
#include <vector>

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
  reset(const Teuchos::RCP<std::vector<Teuchos::RCP<NOX::Solver::Generic> > >& solvers, 
	const Teuchos::RCP<NOX::Multiphysics::DataExchange::Interface>& interface, 
	const Teuchos::RCP<NOX::StatusTest::Generic>& tests, 
	const Teuchos::RCP<Teuchos::ParameterList>& params) = 0;

  //! reset methods inherited from NOX::Solver::Generic and needed here to avoid hiding this overloaded virtual method

  virtual void 
  reset(const NOX::Abstract::Vector& initialGuess) = 0;

  virtual void 
  reset(const NOX::Abstract::Vector& initialGuess, 
	const Teuchos::RCP<NOX::StatusTest::Generic>& tests) = 0;

};
} // namespace Solver
} // namespace Multiphysics
} // namespace NOX


#endif
