// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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

#ifndef LOCA_STATUSTEST_ABSTRACT_H
#define LOCA_STATUSTEST_ABSTRACT_H

// FIXME Replace by LOCA equivalent
#include "NOX_Common.H"         // for std::ostream
#include "LOCA_Abstract_Iterator.H"
#include <ostream>

namespace LOCA {


  // forward declaration of LOCA::Stepper
  // Modify interfaces using LOCA::Stepper in favor of its LOCA::Abstract::Iterator parent (GAH)
//  class Stepper;

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
    checkType is LOCA::StatusType::None.  If the test is skipped, then
    the status should be set to LOCA::StatusTest::Unevaluated.
  */
  // Modify interfaces using LOCA::Stepper in favor of its LOCA::Abstract::Iterator parent (GAH)
  virtual LOCA::StatusTest::StatusType
//  checkStatus(const LOCA::Stepper& stepper,
  checkStatus(const LOCA::Abstract::Iterator& stepper,
		    LOCA::StatusTest::CheckType checkType) = 0;

  //! Return the result of the most recent checkStatus call
  virtual LOCA::StatusTest::StatusType
  getStatus() const = 0;

  //! Output formatted description of stopping test to output stream.
  virtual std::ostream&
  print(std::ostream& stream, int indent = 0) const = 0;

};

std::ostream&
operator<<(std::ostream& os, LOCA::StatusTest::StatusType status);

} // namespace StatusTest
} // namespace LOCA


#endif
