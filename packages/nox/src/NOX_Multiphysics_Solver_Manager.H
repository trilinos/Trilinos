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

#ifndef _NOX_MULTIPHYSICS_SOLVER_MANAGER_H
#define _NOX_MULTIPHYSICS_SOLVER_MANAGER_H

#include "NOX_Solver_Generic.H"	// base class
#include "NOX_Multiphysics_Solver_Generic.H"	// base class
#include "NOX_Common.H"		// class data element (std::string)
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

class Manager : public NOX::Solver::Generic {

public:

  //! Empty constructor - reset called later to really construct it
  Manager();

  /*!
    \brief Constructor

    See reset() for a full description.
  */
  Manager(const Teuchos::RCP<std::vector<Teuchos::RCP<NOX::Solver::Generic> > >& solvers, 
          const Teuchos::RCP<NOX::Multiphysics::DataExchange::Interface>& i, 
	  const Teuchos::RCP<NOX::StatusTest::Generic>& t, 
	  const Teuchos::RCP<Teuchos::ParameterList>& p);

  /*!
    \brief Constructor

    See reset() for a full description.
  */
  Manager(const Teuchos::RCP<NOX::Abstract::Group>& grp, 
	  const Teuchos::RCP<NOX::StatusTest::Generic>& t, 
	  const Teuchos::RCP<Teuchos::ParameterList>& p);

  //! Destructor.
  virtual ~Manager();

  virtual bool reset(const Teuchos::RCP<std::vector<Teuchos::RCP<NOX::Solver::Generic> > >& solvers, 
                     const Teuchos::RCP<NOX::Multiphysics::DataExchange::Interface>& i, 
	             const Teuchos::RCP<NOX::StatusTest::Generic>& tests, 
	             const Teuchos::RCP<Teuchos::ParameterList>& params);
  
  virtual void reset(const NOX::Abstract::Vector& initialGuess,
	             const Teuchos::RCP<NOX::StatusTest::Generic>& tests);
  virtual void reset(const NOX::Abstract::Vector& initialGuess);
  virtual NOX::StatusTest::StatusType getStatus();
  virtual NOX::StatusTest::StatusType step();
  virtual NOX::StatusTest::StatusType solve();
  virtual const NOX::Abstract::Group& getSolutionGroup() const;
  virtual const NOX::Abstract::Group& getPreviousSolutionGroup() const;
  virtual int getNumIterations() const;
  virtual const Teuchos::ParameterList& getList() const;

  virtual Teuchos::RCP< const NOX::Abstract::Group > getSolutionGroupPtr() const;
  virtual Teuchos::RCP< const NOX::Abstract::Group > getPreviousSolutionGroupPtr() const;
  virtual Teuchos::RCP< const Teuchos::ParameterList > getListPtr() const;

private:

  //! Print a warning message that oldName is deprecated and newName should be used instead
  virtual void deprecated(const std::string& oldName, const std::string& newName) const;

  /*! 
    \brief Check that ptr is non-null. If it is null, print a
    message that prints fname and throw an error. 
  */
  virtual void checkNullPtr(const std::string& fname) const;

private:

  //! Printing utilities
  NOX::Utils utils;

  //! Name of method being used
  std::string method;


protected:

  //! Pointer to the coupling solver object
  NOX::Multiphysics::Solver::Generic* cplPtr;


};
}
}
}
#endif
