
//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
//

#ifndef BELOS_STATUS_TEST_HPP
#define BELOS_STATUS_TEST_HPP

/*!
  \file BelosStatusTest.hpp
  \brief Pure virtual base class for defining the status testing capabilities of Belos.
*/

#include "BelosTypes.hpp"
#include "BelosIteration.hpp"
#include "BelosConfigDefs.hpp"
#include "Teuchos_Describable.hpp"

  /*! 
    \class Belos::StatusTest
    \brief A pure virtual class for defining the status tests for the Belos iterative solvers
    
    Belos::StatusTest is an interface that can be implemented to create convergence tests for
    all Belos solvers.  Almost any kind of test can be expressed using this mechanism, 
    including composite tests (see Belos::StatusTestCombo). 
  */

namespace Belos {

  //! @name StatusTest Exceptions
  //@{

  /** \brief Exception thrown to signal error in a status test during Belos::StatusTest::checkStatus().
   */
  class StatusTestError : public BelosError 
  {public: StatusTestError(const std::string& what_arg) : BelosError(what_arg) {}};

  class StatusTestNaNError : public StatusTestError
  {public: StatusTestNaNError(const std::string& what_arg) : StatusTestError(what_arg) {}};

  //@}

template <class ScalarType, class MV, class OP>
class StatusTest : public Teuchos::Describable {

 public:
   //! @name Constructors/destructors
  //@{ 

  //! Constructor
  StatusTest() {};

  //! Destructor
  virtual ~StatusTest() {};
  //@}

  //! @name Status methods
  //@{ 
  //! Check convergence status: Unconverged, Converged, Failed.
  /*! This method checks to see if the convergence criteria are met.  The calling routine may pass in the
    current native residual std::vector (the one naturally produced as part of the iterative method) or a
    pre-computed estimate of the two-norm of the current residual, or both or neither.  The calling routine
    should also indicate if the solution of the linear problem has been updated to be compatible with
    the residual.  Some methods, such as GMRES do not update the solution at each iteration.

    \return Belos::StatusType: Unconverged, Converged or Failed.
  */
  virtual StatusType checkStatus( Iteration<ScalarType,MV,OP>* iSolver ) = 0;

  //! Return the result of the most recent CheckStatus call.
  virtual StatusType getStatus() const = 0;
  //@}

  //! @name Reset methods
  //@{ 
  //! Informs the convergence test that it should reset its internal configuration to the initialized state.
  /*! This is necessary for the case when the status test is being reused by another solver or for another
    linear problem.  The status test may have information that pertains to a particular linear system.  The
    internal information will be reset back to the initialized state.  The user specified information that
    the convergence test uses will remain.
  */
  virtual void reset() = 0;
  //@}

  //! @name Print methods
  //@{ 

  //! Output formatted description of stopping test to output stream.
  virtual void print(std::ostream& os, int indent = 0) const = 0;
 
  //! Output the result of the most recent CheckStatus call.
  virtual void printStatus(std::ostream& os, StatusType type) const {
    os << std::left << std::setw(13) << std::setfill('.');
    switch (type) {
    case  Passed:
      os << "Passed";
      break;
    case  Failed:
      os << "Failed";
      break;
    case  Undefined:  
    default:
      os << "**";
      break;
    }
    os << std::left << std::setfill(' ');
    return;
  };
  //@}

};

} // end of Belos namespace

#endif /* BELOS_STATUS_TEST_HPP */
