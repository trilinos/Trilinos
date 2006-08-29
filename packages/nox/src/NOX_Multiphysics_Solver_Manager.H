//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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

#ifndef _NOX_MULTIPHYSICS_SOLVER_MANAGER_H
#define _NOX_MULTIPHYSICS_SOLVER_MANAGER_H

#include "NOX_Solver_Manager.H"	// base class
#include "NOX_Multiphysics_Solver_Generic.H"	// base class
#include "NOX_Common.H"		// class data element (string)
#include "NOX_Utils.H"          // class data element

namespace NOX {
namespace Multiphysics {
namespace Solver{

/*! 

  \brief %Manager class to control the instantiation of the objects
  derived from the NOX::Solver::Generic object.

  <B>Parameters</B>
  
  The following entries may be specified in the parameter list.

   <ul>
   <li> "Nonlinear %Solver" - Name of the solver method. Valid choices are
   <ul> 
   <li> "Line Search Based" (NOX::Solver::LineSearchBased) [Default]
   <li> "Trust Region Based" (NOX::Solver::TrustRegionBased)
   </ul>
   </ul>
   
   \deprecated The "Nonlinear %Solver" choices "Newton" and "Line
   Search" are deprecated and revert to "Line Search Based". Likewise,
   the choice "Trust Region" is deprecated and reverts to "Trust
   Region Based".

  \author Russell Hooper (SNL 1416)
*/

class Manager : public NOX::Solver::Manager {

public:

  //! Empty constructor - reset called later to really construct it
  Manager();

  /*!
    \brief Constructor

    See reset() for a full description.
  */
  Manager(const Teuchos::RefCountPtr<vector<Teuchos::RefCountPtr<NOX::Solver::Manager> > >& solvers, 
          const Teuchos::RefCountPtr<NOX::Multiphysics::DataExchange::Interface>& i, 
	  const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& t, 
	  const Teuchos::RefCountPtr<Teuchos::ParameterList>& p);

  /*!
    \brief Constructor

    See reset() for a full description.
  */
  Manager(const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp, 
	  const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& t, 
	  const Teuchos::RefCountPtr<Teuchos::ParameterList>& p);

  //! Destructor.
  virtual ~Manager();

  virtual bool reset(const Teuchos::RefCountPtr<vector<Teuchos::RefCountPtr<NOX::Solver::Manager> > >& solvers, 
                     const Teuchos::RefCountPtr<NOX::Multiphysics::DataExchange::Interface>& i, 
	             const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& tests, 
	             const Teuchos::RefCountPtr<Teuchos::ParameterList>& params);
  
  virtual bool reset(const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp, 
	             const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& tests, 
	             const Teuchos::RefCountPtr<Teuchos::ParameterList>& params);
  
  virtual bool reset(const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp,
	             const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& tests);

  virtual bool reset(const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp);

  virtual NOX::StatusTest::StatusType getStatus();

  virtual NOX::StatusTest::StatusType step();

  virtual NOX::StatusTest::StatusType solve();

  virtual const NOX::Abstract::Group& getSolutionGroup() const;
  
  virtual const NOX::Abstract::Group& getPreviousSolutionGroup() const;

  virtual int getNumIterations() const;

  virtual const Teuchos::ParameterList& getList() const;

private:

  //! Print a warning message that oldName is deprecated and newName should be used instead
  virtual void deprecated(const string& oldName, const string& newName) const;

  /*! 
    \brief Check that ptr is non-null. If it is null, print a
    message that prints fname and throw an error. 
  */
  virtual void checkNullPtr(const string& fname) const;

private:

  //! Printing utilities
  NOX::Utils utils;

  //! Name of method being used
  string method;


protected:

  //! Pointer to the coupling solver object
  NOX::Multiphysics::Solver::Generic* cplPtr;


};
}
}
}
#endif
