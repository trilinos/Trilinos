// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_ERRORPOLLING_H
#define TEUCHOS_ERRORPOLLING_H

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_TestForException.hpp"

/*! \defgroup ErrorPolling_grp Utility code for synchronizing std::exception detection across processors. 
*/
//@{

namespace Teuchos
{
  class MPIComm;

  /** \brief ErrorPolling provides utilities for establishing agreement
   * between processors on whether an std::exception has been detected on any one
   * processor.
   *
   * The two functions must be used in a coordinated way. The simplest use
   * case is to embed a call to reportFailure() whenever an std::exception is
   * detected at the top-level try/catch block, and then to do a call to
   * pollForFailures() whenever it is desired to check for off-processor
   * errors before proceeding. The macro

    \code
    TEUCHOS_TEST_FOR_FAILURE(comm);
    \endcode  

   * calls pollForFailures() and throws an std::exception if the return value is
   * true.
   *
   * Polling is a collective operation (an MPI_Reduce) and so incurs some
   * performance overhead. It can be disabled with a call to 
   * \code
   * Teuchos::ErrorPolling::disable();
   * \endcode 
   * IMPORTANT: all processors must agree on whether collective error checking
   * is enabled or disabled. If there are inconsistent states, the reduction
   * operations in pollForFailures() will hang because some processors cannot be 
   * contacted. 
   */
  class TEUCHOS_LIB_DLL_EXPORT ErrorPolling
  {
  public:
    /** Call this function upon catching an std::exception in order to
     * inform other processors of the error. This function will do an
     * AllReduce in conjunction with calls to either this function or
     * its partner, pollForFailures(), on the other processors. This
     * procedure has the effect of communicating to the other
     * processors that an std::exception has been detected on this one. */
    static void reportFailure(const MPIComm& comm);
    
    /** Call this function after std::exception-free completion of a
     * try/catch block. This function will do an AllReduce in
     * conjunction with calls to either this function or its partner,
     * reportFailure(), on the other processors. If a failure has been
     * reported by another processor, the call to pollForFailures()
     * will return true and an std::exception can be thrown. */
    static bool pollForFailures(const MPIComm& comm);
    
    /** Activate error polling */
    static void enable() {isActive()=true;}

    /** Disable error polling */
    static void disable() {isActive()=false;}

  private:
    /** Set or check whether error polling is active */
    static bool& isActive() {static bool rtn = true; return rtn;}
  };

  /** 
   * This macro polls all processors in the given communicator to find
   * out whether an error has been reported by a call to 
   * ErrorPolling::reportFailure(comm).
   * 
   * @param comm [in] The communicator on which polling will be done
   */
#define TEUCHOS_POLL_FOR_FAILURES(comm)                                  \
  TEST_FOR_EXCEPTION(Teuchos::ErrorPolling::pollForFailures(comm), \
                     std::runtime_error,                                     \
                     "off-processor error detected by proc=" << (comm).getRank());
}

//@}

#endif
