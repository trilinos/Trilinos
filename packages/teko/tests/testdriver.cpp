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
#include "src/tJacobi2x2PreconditionerFactory.hpp"
#include "src/tBlockJacobiPreconditionerFactory.hpp"
#include "src/tBlockUpperTriInverseOp.hpp"
#include "src/Epetra/tEpetraLSCIntegrationTest.hpp"
#include "src/Epetra/tEpetraOperatorWrapper.hpp"
#include "src/Epetra/tStridedEpetraOperator.hpp"
#include "src/Epetra/tInterlacedEpetra.hpp"
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

#if 1
   PB_ADD_UNIT_TEST(PB::Test::tLU2x2PreconditionerFactory,LU2x2PreconditionerFactory);
   PB_ADD_UNIT_TEST(PB::Test::tLSCStablePreconditionerFactory,LSCStablePreconditionerFactory);
   PB_ADD_UNIT_TEST(PB::Test::tJacobi2x2PreconditionerFactory,Jacobi2x2PreconditionerFactory);
   PB_ADD_UNIT_TEST(PB::Test::tBlockJacobiPreconditionerFactory,BlockJacobiPreconditionerFactory);
   PB_ADD_UNIT_TEST(PB::Test::tBlockUpperTriInverseOp,BlockUpperTriInverseOp);
   PB_ADD_UNIT_TEST(PB::Test::tEpetraOperatorWrapper,EpetraOperatorWrapper);
   PB_ADD_UNIT_TEST(PB::Test::tInterlacedEpetra,InterlacedEpetra);
   PB_ADD_UNIT_TEST(PB::Test::tEpetraThyraConverter,EpetraThyraConverter);
#endif
   if(not isfast) {
      PB_ADD_UNIT_TEST(PB::Test::tStridedEpetraOperator,tStridedEpetraOperator);
      PB_ADD_UNIT_TEST(PB::Test::tEpetraLSCIntegrationTest,EpetraLSCIntegrationTest);
   }

   bool status = PB::Test::UnitTest::RunTests(verbosity,*termout,*failout);

   return status ? 0 : -1;
}
