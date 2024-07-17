// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _NOX_MULTIPHYSICS_SOLVER_MANAGER_H
#define _NOX_MULTIPHYSICS_SOLVER_MANAGER_H

#include "NOX_Solver_Generic.H"    // base class
#include "NOX_Multiphysics_Solver_Generic.H"    // base class
#include "NOX_Common.H"        // class data element (std::string)
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
  virtual void reset();
  virtual NOX::StatusTest::StatusType getStatus() const;
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
