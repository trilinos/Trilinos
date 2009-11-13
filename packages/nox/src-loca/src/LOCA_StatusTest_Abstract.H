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

#ifndef LOCA_STATUSTEST_ABSTRACT_H
#define LOCA_STATUSTEST_ABSTRACT_H

// FIXME Replace by LOCA equivalent
#include "NOX_Common.H"         // for ostream

namespace LOCA {

  // forward declaration of LOCA::Stepper
  class Stepper;

//! %Status checkers.
namespace StatusTest {

/*! 

  \brief %Generic status test to check for convergence or failure of
  the LOCA continuation.

*/

  //! Enumerated type for status of the iterator
  enum StatusType {
    Finished = 1,       //!< The iterator is finished
    Failed = 0,         //!< The iterator has failed
    NotFinished = -1,   //!< The iterator is not finished
    Unevaluated = -2    //!< The iterator status has not been evaluated
  };

  //! Type of check that should be done by checkStatus
  enum CheckType
    {
      //! Evaluate every test and subtest
      Complete,
      //! Evaluate all tests and subtests that are necessary to determine the status
      Minimal,
      //! Evaluation is entirely optional
      None
    };

class Abstract {

public:

  //! Constructor.
  Abstract() {};

  //! Destructor.
  virtual ~Abstract() {};

  //! %Test the stopping criterion
  /*!  
    The test can (and should, if possible) be skipped if 
    checkType is NOX::StatusType::None.  If the test is skipped, then
    the status should be set to NOX::StatusTest::Unevaluated.
  */
  virtual LOCA::StatusTest::StatusType
  checkStatus(const LOCA::Stepper& stepper,
		    LOCA::StatusTest::CheckType checkType) = 0;

  //! Return the result of the most recent checkStatus call
  virtual LOCA::StatusTest::StatusType
  getStatus() const = 0;

  //! Output formatted description of stopping test to output stream.
  virtual ostream&
  print(ostream& stream, int indent = 0) const = 0;

};

} // namespace StatusTest
} // namespace LOCA

ostream&
operator<<(ostream& os, LOCA::StatusTest::StatusType status);


#endif

