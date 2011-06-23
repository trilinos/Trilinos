#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_ParameterList.hpp>

#include <string>
#include <iostream>
#include <vector>
#include <set>

#include "Panzer_UniqueGlobalIndexer_Utilities.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_GeometricAggFieldPattern.hpp"

#include "UnitTest_UniqueGlobalIndexer.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
   #include "mpi.h"
#else
   #include "Epetra_SerialComm.h"
#endif

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {

TEUCHOS_UNIT_TEST(tUniqueGlobalIndexer_Utilities,GhostedFieldVector)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      RCP<Epetra_Comm> eComm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else
      RCP<Epetra_Comm> eComm = rcp(new Epetra_SerialComm());
   #endif

   Teuchos::FancyOStream fout(Teuchos::rcpFromRef(std::cout));
   fout.setOutputToRootOnly(-1); 
   fout.setShowProcRank(true); 

   int myRank = eComm->MyPID();
   int numProcs = eComm->NumProc();

   TEUCHOS_ASSERT(numProcs==2);

   RCP<panzer::UniqueGlobalIndexer<short,int> > globalIndexer 
         = rcp(new panzer::unit_test::UniqueGlobalIndexer(myRank,numProcs));

   std::vector<int> ghostedFields;
   std::vector<int> sharedIndices;
   globalIndexer->getOwnedAndSharedIndices(sharedIndices);
   panzer::buildGhostedFieldVector(*globalIndexer,ghostedFields);

   TEST_EQUALITY(ghostedFields.size(),sharedIndices.size());
   TEST_COMPARE(*std::min_element(ghostedFields.begin(),ghostedFields.end()),>,-1);

   std::stringstream ss;
   ss << "Field Numbers = ";
   for(std::size_t i=0;i<ghostedFields.size();i++) 
      ss << sharedIndices[i] << ":" << ghostedFields[i] << " ";
   out << std::endl;
   out << ss.str() << std::endl;

   if(myRank==0) {
      TEST_EQUALITY(ghostedFields[0], 0);
      TEST_EQUALITY(ghostedFields[1], 1);
      TEST_EQUALITY(ghostedFields[2], 0);
      TEST_EQUALITY(ghostedFields[3], 1);
      TEST_EQUALITY(ghostedFields[4], 0);
      TEST_EQUALITY(ghostedFields[5], 1);
      TEST_EQUALITY(ghostedFields[6], 0);
      TEST_EQUALITY(ghostedFields[7], 1);
      TEST_EQUALITY(ghostedFields[8], 0);
      TEST_EQUALITY(ghostedFields[9], 1);
      TEST_EQUALITY(ghostedFields[10],0);
      TEST_EQUALITY(ghostedFields[11],0);
      TEST_EQUALITY(ghostedFields[12],0);
      TEST_EQUALITY(ghostedFields[13],1);
   }
   else if(myRank==1) {
      TEST_EQUALITY(ghostedFields[0], 0);
      TEST_EQUALITY(ghostedFields[1], 1);
      TEST_EQUALITY(ghostedFields[2], 0);
      TEST_EQUALITY(ghostedFields[3], 1);
      TEST_EQUALITY(ghostedFields[4], 0);
      TEST_EQUALITY(ghostedFields[5], 1);
      TEST_EQUALITY(ghostedFields[6], 0);
      TEST_EQUALITY(ghostedFields[7], 1);
      TEST_EQUALITY(ghostedFields[8], 0);
      TEST_EQUALITY(ghostedFields[9], 0);
      TEST_EQUALITY(ghostedFields[10],0);
      TEST_EQUALITY(ghostedFields[11],0);
   }
   else 
      TEUCHOS_ASSERT(false);
}

}
