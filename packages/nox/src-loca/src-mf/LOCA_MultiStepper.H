// $Id$ 
// $Source$ 

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

#ifndef LOCA_STEPPERMULTI_H
#define LOCA_STEPPERMULTI_H

#include "Teuchos_RefCountPtr.hpp"

#include "LOCA_Abstract_Iterator.H"    // for status
#include "Teuchos_ParameterList.hpp"	       // class data element
#include "NOX_Solver_Manager.H"	       // class data element
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
    MultiStepper(const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		 const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& initialGuess,
	       const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& t,
	       const Teuchos::RefCountPtr<Teuchos::ParameterList>& p);
    
    //! Destructor
    virtual ~MultiStepper();

    //! Reset the Stepper to start a new continuation run
    virtual bool reset(
	  const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	  const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& initialGuess,
	  const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& t,
	  const Teuchos::RefCountPtr<Teuchos::ParameterList>& p);

    //! Run continuation
    virtual LOCA::Abstract::Iterator::IteratorStatus run();

    //! Return the current solution group.
    virtual Teuchos::RefCountPtr<const LOCA::MultiContinuation::AbstractGroup> 
    getSolutionGroup();

    //! Return the output parameters from the stepper algorithm. 
    virtual Teuchos::RefCountPtr<const Teuchos::ParameterList> 
    getList() const;

  protected:

    //! Print to the screen the initialization information
    virtual void printInitializationInfo();

    //! Parses parameter list to get continuation parameter data
    virtual void getConParamData();

    //! Builds a constrained group is specified through parameter lists
    Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>
    buildConstrainedGroup(
      const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& grp);

  private:

    //! Private to prohibit copying
    MultiStepper(const MultiStepper&);

    //! Private to prohibit copying
    MultiStepper& operator = (const MultiStepper&);

  protected:

    //! Global data
    Teuchos::RefCountPtr<LOCA::GlobalData> globalData;

    //! Parsed parameter list
    Teuchos::RefCountPtr<LOCA::Parameter::SublistParser> parsedParams;

    //! Predictor strategy
    Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy> predictor;

    //! Current continuation group
    Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractStrategy> curGroupPtr;

    //! Pointer to bifurcation group
    Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup> bifGroupPtr;

    //! Pointer to the status test.
    Teuchos::RefCountPtr<NOX::StatusTest::Generic> statusTestPtr;

    //! Pointer to parameter list passed in constructor/reset method
    Teuchos::RefCountPtr<Teuchos::ParameterList> paramListPtr;

    //! Pointer to stepper parameters
    Teuchos::RefCountPtr<Teuchos::ParameterList> stepperList;

    //! Solver (a.k.a corrector)
    Teuchos::RefCountPtr<NOX::Solver::Manager> solverPtr;

    //! Stores continuation parameter vector
    LOCA::ParameterVector paramVec;

    //! Stores vector of continuation parameter IDs
    vector<int> conParamIDVec;

    //! Stores list of continuation parameter data
    list<ParamData> conParamData;
  };
} // end namespace LOCA

#endif

