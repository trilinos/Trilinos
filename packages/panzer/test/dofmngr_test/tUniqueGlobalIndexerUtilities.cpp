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

void fillFieldContainer(int fieldNum,const std::string & blockId,
                        const panzer::UniqueGlobalIndexer<short,int> & ugi,
                        Intrepid::FieldContainer<int> & data)
{
   data.resize(1,4);

   const std::vector<short> & elements = ugi.getElementBlock(blockId);
   const std::vector<int> & fieldOffsets = ugi.getGIDFieldOffsets(blockId,fieldNum);
   std::vector<int> gids;
   for(std::size_t e=0;e<elements.size();e++) {
      ugi.getElementGIDs(elements[e],gids);
      for(std::size_t f=0;f<fieldOffsets.size();f++)
         data(e,f) = gids[fieldOffsets[f]];
   }
}

TEUCHOS_UNIT_TEST(tUniqueGlobalIndexer_Utilities,updateGhostedDataVector)
{
   typedef Intrepid::FieldContainer<int> IntFieldContainer;

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      RCP<Epetra_Comm> eComm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else
      RCP<Epetra_Comm> eComm = rcp(new Epetra_SerialComm());
   #endif

   int myRank = eComm->MyPID();
   int numProcs = eComm->NumProc();

   TEUCHOS_ASSERT(numProcs==2);

   RCP<panzer::UniqueGlobalIndexer<short,int> > globalIndexer 
         = rcp(new panzer::unit_test::UniqueGlobalIndexer(myRank,numProcs));

   int uFieldNum = globalIndexer->getFieldNum("U");
   int tFieldNum = globalIndexer->getFieldNum("T");

   Teuchos::RCP<Tpetra::Vector<int,std::size_t,int> > reducedFieldVector 
         = panzer::buildGhostedFieldReducedVector(*globalIndexer);

   Tpetra::Vector<int,std::size_t,int> reducedUDataVector(getFieldMap(uFieldNum,*reducedFieldVector));
   Tpetra::Vector<int,std::size_t,int> reducedTDataVector(getFieldMap(tFieldNum,*reducedFieldVector));

   TEST_EQUALITY(reducedUDataVector.getLocalLength(),8);
   TEST_EQUALITY(reducedTDataVector.getLocalLength(),4);

   IntFieldContainer dataU_b0, dataU_b1; 
   fillFieldContainer(uFieldNum,"block_0",*globalIndexer,dataU_b0);
   fillFieldContainer(uFieldNum,"block_1",*globalIndexer,dataU_b1);

   IntFieldContainer dataT_b0; 
   fillFieldContainer(tFieldNum,"block_0",*globalIndexer,dataT_b0);

   updateGhostedDataReducedVector("U","block_0",*globalIndexer,dataU_b0,reducedUDataVector);
   updateGhostedDataReducedVector("U","block_1",*globalIndexer,dataU_b1,reducedUDataVector);

   updateGhostedDataReducedVector("T","block_0",*globalIndexer,dataT_b0,reducedTDataVector);

   std::vector<int> ghostedFields_u(reducedUDataVector.getLocalLength());
   std::vector<int> ghostedFields_t(reducedTDataVector.getLocalLength());

   reducedUDataVector.get1dCopy(Teuchos::arrayViewFromVector(ghostedFields_u));
   reducedTDataVector.get1dCopy(Teuchos::arrayViewFromVector(ghostedFields_t));
   
   if(myRank==0) {
      TEST_EQUALITY(ghostedFields_u[0], 0);
      TEST_EQUALITY(ghostedFields_u[1], 2);
      TEST_EQUALITY(ghostedFields_u[2], 4);
      TEST_EQUALITY(ghostedFields_u[3], 6);
      TEST_EQUALITY(ghostedFields_u[4], 8);
      TEST_EQUALITY(ghostedFields_u[5],12);
      TEST_EQUALITY(ghostedFields_u[6],13);
      TEST_EQUALITY(ghostedFields_u[7],10);

      TEST_EQUALITY(ghostedFields_t[0], 1);
      TEST_EQUALITY(ghostedFields_t[1], 3);
      TEST_EQUALITY(ghostedFields_t[2], 5);
      TEST_EQUALITY(ghostedFields_t[3], 7);
   }
   else if(myRank==1) {
      TEST_EQUALITY(ghostedFields_u[0], 2);
      TEST_EQUALITY(ghostedFields_u[1], 8);
      TEST_EQUALITY(ghostedFields_u[2],10);
      TEST_EQUALITY(ghostedFields_u[3], 4);
      TEST_EQUALITY(ghostedFields_u[4],12);
      TEST_EQUALITY(ghostedFields_u[5],14);
      TEST_EQUALITY(ghostedFields_u[6],15);
      TEST_EQUALITY(ghostedFields_u[7],13);

      TEST_EQUALITY(ghostedFields_t[0], 3);
      TEST_EQUALITY(ghostedFields_t[1], 9);
      TEST_EQUALITY(ghostedFields_t[2],11);
      TEST_EQUALITY(ghostedFields_t[3], 5);
   }
   else 
      TEUCHOS_ASSERT(false);
}

}
