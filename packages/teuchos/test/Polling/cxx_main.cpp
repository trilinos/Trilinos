// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
