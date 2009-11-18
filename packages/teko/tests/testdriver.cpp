#include <iostream>
#include <fstream>
#include <string>

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_MPISession.hpp"
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
#include "src/Epetra/tEpetraOperatorWrapper.hpp"
#include "src/Epetra/tStridedEpetraOperator.hpp"
#include "src/Epetra/tInterlacedEpetra.hpp"
#include "src/Epetra/tBlockingEpetra.hpp"
#include "src/Epetra/tBlockedEpetraOperator.hpp"
#include "src/Epetra/tEpetraThyraConverter.hpp"


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

   PB::Test::UnitTest::SetComm(Teuchos::rcpFromRef(Comm));

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

   PB_ADD_UNIT_TEST(PB::Test::tLU2x2PreconditionerFactory,LU2x2PreconditionerFactory);
   PB_ADD_UNIT_TEST(PB::Test::tLSCStablePreconditionerFactory,LSCStablePreconditionerFactory);
   PB_ADD_UNIT_TEST(PB::Test::tLSCStabilized,LSCStabilized);
   PB_ADD_UNIT_TEST(PB::Test::tJacobi2x2PreconditionerFactory,Jacobi2x2PreconditionerFactory);
   PB_ADD_UNIT_TEST(PB::Test::tBlockJacobiPreconditionerFactory,BlockJacobiPreconditionerFactory);
   PB_ADD_UNIT_TEST(PB::Test::tBlockUpperTriInverseOp,BlockUpperTriInverseOp);
   PB_ADD_UNIT_TEST(PB::Test::tBlockLowerTriInverseOp,BlockLowerTriInverseOp);
   PB_ADD_UNIT_TEST(PB::Test::tEpetraOperatorWrapper,EpetraOperatorWrapper);
   PB_ADD_UNIT_TEST(PB::Test::tInterlacedEpetra,InterlacedEpetra);
   PB_ADD_UNIT_TEST(PB::Test::tBlockingEpetra,BlockingEpetra);
   PB_ADD_UNIT_TEST(PB::Test::tEpetraThyraConverter,EpetraThyraConverter);
   PB_ADD_UNIT_TEST(PB::Test::tGraphLaplacian,tGraphLaplacian);
   PB_ADD_UNIT_TEST(PB::Test::tParallelInverse,tParallelInverse);
   PB_ADD_UNIT_TEST(PB::Test::tExplicitOps,tExplicitOps);
   PB_ADD_UNIT_TEST(PB::Test::tSIMPLEPreconditionerFactory,SIMPLEPreconditionerFactory);
   PB_ADD_UNIT_TEST(PB::Test::tLSCHIntegrationTest,LSCHIntegrationTest);
   PB_ADD_UNIT_TEST(PB::Test::tLumping,Lumping);
   PB_ADD_UNIT_TEST(PB::Test::tAbsRowSum,AbsRowSum);
   if(not isfast) {
      PB_ADD_UNIT_TEST(PB::Test::tLSCIntegrationTest,LSCIntegrationTest);
      PB_ADD_UNIT_TEST(PB::Test::tStridedEpetraOperator,tStridedEpetraOperator);
      PB_ADD_UNIT_TEST(PB::Test::tBlockedEpetraOperator,tBlockedEpetraOperator);
   }

   bool status = PB::Test::UnitTest::RunTests(verbosity,*termout,*failout);


   if(not status)
      *termout << "Teko tests failed" << std::endl; 
   return status ? 0 : -1;
}
