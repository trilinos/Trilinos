/*
// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
// 
// ***********************************************************************
// 
// @HEADER

*/

#include <iostream>
#include <fstream>
#include <string>
#include <sys/types.h>
#include <unistd.h>


#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RCP.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
   #include "mpi.h"
#else
   #include "Epetra_SerialComm.h"
#endif

#include "Test_Utils.hpp"
#include "src/tDiagonalPreconditionerFactory.hpp"
#include "src/tLU2x2PreconditionerFactory.hpp"
#include "src/tLSCStablePreconditionerFactory.hpp"
#include "src/tSIMPLEPreconditionerFactory.hpp"
#include "src/tLSCStabilized.hpp"
#include "src/tJacobi2x2PreconditionerFactory.hpp"
#include "src/tBlockJacobiPreconditionerFactory.hpp"
#include "src/tBlockUpperTriInverseOp.hpp"
#include "src/tBlockLowerTriInverseOp.hpp"
#include "src/tLSCIntegrationTest.hpp"
#include "src/tLSCHIntegrationTest.hpp"
#include "src/tGraphLaplacian.hpp"
#include "src/tParallelInverse.hpp"
#include "src/tExplicitOps.hpp"
#include "src/tLumping.hpp"
#include "src/tAbsRowSum.hpp"
#include "src/tNeumannSeries.hpp"
#include "src/tPCDStrategy.hpp"
#include "src/Epetra/tEpetraOperatorWrapper.hpp"
#include "src/Epetra/tStridedEpetraOperator.hpp"
#include "src/Epetra/tInterlacedEpetra.hpp"
#include "src/Epetra/tBlockingEpetra.hpp"
#include "src/Epetra/tBlockedEpetraOperator.hpp"
#include "src/Epetra/tEpetraThyraConverter.hpp"

void gdbIn()
{
    int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("PID %d on %s ready for attach\n", getpid(), hostname);
    fflush(stdout);
    while (0 == i)
        sleep(5);
}

int main(int argc,char * argv[])
{
   // calls MPI_Init and MPI_Finalize
   Teuchos::GlobalMPISession mpiSession(&argc,&argv);

   // build MPI/Serial communicator
   #ifdef HAVE_MPI
      Epetra_MpiComm Comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm Comm;
   #endif

   Teko::Test::UnitTest::SetComm(Teuchos::rcpFromRef(Comm));

   Teuchos::CommandLineProcessor clp;

   int verbosity = 1;
   std::string faillog = "failure.log";
   bool isfast = false;

   clp.setOption("verb",&verbosity,"How verbose is the output? 1 is normal 10 is a lot.");
   clp.setOption("log",&faillog,"File for failure information to go to (also high verbosity text)");
   clp.setOption("fast","notfast",&isfast,"Run only fast tests");
   clp.parse(argc,argv);

   Teuchos::RCP<Teuchos::FancyOStream> termout = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
   Teuchos::RCP<Teuchos::FancyOStream> failout;
   std::ofstream failure;

   if(faillog=="stdout") {
      failout = termout;
   }
   else {
      failure.open(faillog.c_str());
      failout = Teuchos::getFancyOStream(Teuchos::rcpFromRef(failure));
   }

   termout->setOutputToRootOnly(0);
   failout->setOutputToRootOnly(0);

   // gdbIn();

   Teko_ADD_UNIT_TEST(Teko::Test::tSIMPLEPreconditionerFactory,SIMPLEPreconditionerFactory);
   Teko_ADD_UNIT_TEST(Teko::Test::tDiagonalPreconditionerFactory,DiagonalPreconditionerFactory);
   Teko_ADD_UNIT_TEST(Teko::Test::tLU2x2PreconditionerFactory,LU2x2PreconditionerFactory);
   Teko_ADD_UNIT_TEST(Teko::Test::tLSCStablePreconditionerFactory,LSCStablePreconditionerFactory);
   Teko_ADD_UNIT_TEST(Teko::Test::tLSCStabilized,LSCStabilized);
   Teko_ADD_UNIT_TEST(Teko::Test::tJacobi2x2PreconditionerFactory,Jacobi2x2PreconditionerFactory);
   Teko_ADD_UNIT_TEST(Teko::Test::tBlockJacobiPreconditionerFactory,BlockJacobiPreconditionerFactory);
   Teko_ADD_UNIT_TEST(Teko::Test::tBlockUpperTriInverseOp,BlockUpperTriInverseOp);
   Teko_ADD_UNIT_TEST(Teko::Test::tBlockLowerTriInverseOp,BlockLowerTriInverseOp);
   Teko_ADD_UNIT_TEST(Teko::Test::tEpetraOperatorWrapper,EpetraOperatorWrapper);
   Teko_ADD_UNIT_TEST(Teko::Test::tInterlacedEpetra,InterlacedEpetra);
   Teko_ADD_UNIT_TEST(Teko::Test::tBlockingEpetra,BlockingEpetra);
   Teko_ADD_UNIT_TEST(Teko::Test::tEpetraThyraConverter,EpetraThyraConverter);
   Teko_ADD_UNIT_TEST(Teko::Test::tGraphLaplacian,tGraphLaplacian);
   Teko_ADD_UNIT_TEST(Teko::Test::tParallelInverse,tParallelInverse);
   Teko_ADD_UNIT_TEST(Teko::Test::tExplicitOps,tExplicitOps);
   Teko_ADD_UNIT_TEST(Teko::Test::tLSCHIntegrationTest,LSCHIntegrationTest);
   Teko_ADD_UNIT_TEST(Teko::Test::tLumping,Lumping);
   Teko_ADD_UNIT_TEST(Teko::Test::tAbsRowSum,AbsRowSum);
   Teko_ADD_UNIT_TEST(Teko::Test::tNeumannSeries,NeumannSeries);
   Teko_ADD_UNIT_TEST(Teko::Test::tPCDStrategy,PCDStrategy);
   if(not isfast) {
      Teko_ADD_UNIT_TEST(Teko::Test::tLSCIntegrationTest,LSCIntegrationTest);
      Teko_ADD_UNIT_TEST(Teko::Test::tStridedEpetraOperator,tStridedEpetraOperator);
      Teko_ADD_UNIT_TEST(Teko::Test::tBlockedEpetraOperator,tBlockedEpetraOperator);
   }

   bool status = Teko::Test::UnitTest::RunTests(verbosity,*termout,*failout);


   if(not status)
      *termout << "Teko tests failed" << std::endl; 
   return status ? 0 : -1;
}
