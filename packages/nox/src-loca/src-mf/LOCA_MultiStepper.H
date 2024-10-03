// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_STEPPERMULTI_H
#define LOCA_STEPPERMULTI_H

#include "Teuchos_RCP.hpp"

#include "LOCA_Abstract_Iterator.H"    // for status
#include "Teuchos_ParameterList.hpp"           // class data element
#include "NOX_Solver_Generic.H"           // class data element
#include "NOX_StatusTest_Generic.H"    // for enum StatusType
#include "NOX_Common.H"                // for <string> class

#include "MFLOCA.H"

// Forward Declarations
namespace LOCA {
  class GlobalData;
  namespace MultiContinuation {
    class ExtendedGroup;
    class AbstractGroup;
  }
}

namespace LOCA {

  /*!
    \brief Implementation of LOCA::Abstract::Iterator for computing points
    along a continuation curve.
  */
  /*!
    The %Stepper class implements the pure virtual methods of the
    LOCA::Abstract::Iterator for iteratively computing points along a
    continuation curve.
  */

  class MultiStepper {

  public:

    //! Constructor
    MultiStepper(const Teuchos::RCP<LOCA::GlobalData>& global_data,
         const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& initialGuess,
           const Teuchos::RCP<NOX::StatusTest::Generic>& t,
           const Teuchos::RCP<Teuchos::ParameterList>& p);

    //! Destructor
    virtual ~MultiStepper();

    //! Reset the Stepper to start a new continuation run
    virtual bool reset(
      const Teuchos::RCP<LOCA::GlobalData>& global_data,
      const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& initialGuess,
      const Teuchos::RCP<NOX::StatusTest::Generic>& t,
      const Teuchos::RCP<Teuchos::ParameterList>& p);

    //! Run continuation
    virtual LOCA::Abstract::Iterator::IteratorStatus run();

    //! Return the current solution group.
    virtual Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup>
    getSolutionGroup();

    //! Return the output parameters from the stepper algorithm.
    virtual Teuchos::RCP<const Teuchos::ParameterList>
    getList() const;

  protected:

    //! Print to the screen the initialization information
    virtual void printInitializationInfo();

    //! Parses parameter list to get continuation parameter data
    virtual void getConParamData();

    //! Builds a constrained group is specified through parameter lists
    Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
    buildConstrainedGroup(
      const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp);

  private:

    //! Private to prohibit copying
    MultiStepper(const MultiStepper&);

    //! Private to prohibit copying
    MultiStepper& operator = (const MultiStepper&);

  protected:

    //! Global data
    Teuchos::RCP<LOCA::GlobalData> globalData;

    //! Parsed parameter list
    Teuchos::RCP<LOCA::Parameter::SublistParser> parsedParams;

    //! Predictor strategy
    Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy> predictor;

    //! Current continuation group
    Teuchos::RCP<LOCA::MultiContinuation::AbstractStrategy> curGroupPtr;

    //! Pointer to bifurcation group
    Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> bifGroupPtr;

    //! Pointer to the status test.
    Teuchos::RCP<NOX::StatusTest::Generic> statusTestPtr;

    //! Pointer to parameter list passed in constructor/reset method
    Teuchos::RCP<Teuchos::ParameterList> paramListPtr;

    //! Pointer to stepper parameters
    Teuchos::RCP<Teuchos::ParameterList> stepperList;

    //! Solver (a.k.a corrector)
    Teuchos::RCP<NOX::Solver::Generic> solverPtr;

    //! Stores continuation parameter vector
    LOCA::ParameterVector paramVec;

    //! Stores vector of continuation parameter IDs
    std::vector<int> conParamIDVec;

    //! Stores list of continuation parameter data
    std::list<ParamData> conParamData;
  };
} // end namespace LOCA

#endif

