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

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_MPIComm.hpp"
#include "Teuchos_ErrorPolling.hpp"
#include "Teuchos_Version.hpp"

using namespace Teuchos;
using std::string;

/* \example Test of polling for exceptions on other processors */

int main( int argc, char* argv[] )
{
  /* return value */
  int state=0;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  try
    {

      MPIComm comm = MPIComm::world();

     
      /*----- Demonstrate detection of an off-processor error  -------- */
      
      try
        {
          /* Try some code that will fail on one of the processors */
          try
            {
              /* Generate an std::exception on proc 1 */
              TEST_FOR_EXCEPTION(comm.getRank()==1, std::runtime_error,
                                 "std::exception [expected] detected on proc="
                                 << comm.getRank());
              /* On all other procs, do some calculation */
              double x=0;
              for (int i=0; i<100; i++) x += i;

            }
          catch(std::exception& ex1)
            {
              /* If we catch an std::exception, report the failure to the other 
               * processors. This call to reportFailure() must be
               * paired with a call to pollForFailures() in the 
               * branch that did not detect an std::exception.
               */
              ErrorPolling::reportFailure(comm);
              TEUCHOS_TRACE(ex1);
            }

          /* 
           * Here we poll for the state of other processors. If all processors
           * report OK, pollForFailures() will return zero and an
           * std::exception will not be thrown. If another
           * processor has called reportFailure(), then pollForFailures()
           * will return a nonzero number and an std::exception will be thrown.
           */
          TEST_FOR_EXCEPTION(ErrorPolling::pollForFailures(comm),
                             std::runtime_error, 
                             "off-processor error [expected] detected "
                             "on proc=" << comm.getRank());



          /* Do a collective operation. In the present example,
           * this code should never be reached
           * because all processors should have detected either a local
           * std::exception or a remote std::exception. */
          std::cerr << "this is bad! Processor=" << comm.getRank() 
               << "should not have reached this point" << std::endl;

          /* report the bad news to the testharness 
           * using the return value... */
          state = 1;

          /* Throw an std::exception. This is not a drill!!! */
          TEST_FOR_EXCEPTION(state, std::runtime_error,
                             "std::exception [UNEXPECTED!!!] detected in test "
                             "of polling on processor=" << comm.getRank());

          /* This collective operation would fail if executed here, because
           * one of the processors has thrown an std::exception and never
           * reached this point. Good thing we've polled for errors! */
          int x=comm.getRank();
          int sum;
          comm.allReduce( (void*) &x, (void*) &sum, 1, MPIComm::INT,
                          MPIComm::SUM);
          std::cerr << "sum=" << sum << std::endl;
        }
      catch(std::exception& ex)
        {
          std::cerr << ex.what() << std::endl;
        }

      std::cerr << "p=" << MPIComm::world().getRank() 
           << ": std::exception polling successful" << std::endl;


      /*-- Demonstrate safe pass-through when no off-proc error happens --- */

      try
        {
          /* Try some code that will not fail on any processors */
          try
            {
              /* On all procs, do some foolproof calculation */
              double x=0;
              for (int i=0; i<100; i++) x += i;

            }
          catch(std::exception& ex1)
            {
              /* If we catch an std::exception, report the failure to the other 
               * processors. This call to reportFailure() must be
               * paired with a call to pollForFailures() in the 
               * branch that did not detect an std::exception.
               */
              ErrorPolling::reportFailure(comm);
              TEUCHOS_TRACE(ex1);
            }

          /* 
           * Here we poll for the state of other processors. If all processors
           * report OK, pollForFailures() will return zero and an
           * std::exception will not be thrown. If another
           * processor has called reportFailure(), then pollForFailures()
           * will return a nonzero number and an std::exception will be thrown.
           */
          TEST_FOR_EXCEPTION(ErrorPolling::pollForFailures(comm),
                             std::runtime_error, 
                             "off-processor error [UNEXPECTED!!!] detected "
                             "on proc=" << comm.getRank());



          /* 
           * Do a collective operation. In the present example,
           * this code will be reached on all processors because
           * no std::exception has been thrown by any processor.
           */
          std::cerr << "Processor=" << comm.getRank() 
               << "ready to do collective operation" << std::endl;

          /* 
           * This collective operation is safe because we have polled
           * all processors and known that everyone is still up and running.
           */
          int x=comm.getRank();
          int sum;
          comm.allReduce( (void*) &x, (void*) &sum, 1, MPIComm::INT,
                          MPIComm::SUM);
          if (comm.getRank()==0) std::cerr << "sum=" << sum << std::endl;
        }
      catch(std::exception& ex)
        {
          std::cerr << "std::exception [UNEXPECTED!!!] detected" << std::endl;
          std::cerr << ex.what() << std::endl;
          state = 1;
        }
    }
  catch(std::exception& e)
    {
      std::cerr << e.what() << std::endl;
      state = 1;
    }

  return state;

}
