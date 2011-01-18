#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <string>
#include <iostream>

#include "Panzer_EpetraLinearObjFactory.hpp"
#include "UnitTest_UniqueGlobalIndexer.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"

#include "Thyra_EpetraThyraWrappers.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {

RCP<Epetra_MultiVector> getEpetraMultiVector(RCP<Thyra::MultiVectorBase<double> > & vec,const Epetra_Map & eMap)
{
   return Thyra::get_Epetra_MultiVector(eMap,vec);
}

TEUCHOS_UNIT_TEST(tEpetraLinearObjFactory, vector_constr)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_SerialComm());
   #endif

   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();

   RCP<panzer::UniqueGlobalIndexer<short,int> > indexer 
         = rcp(new unit_test::UniqueGlobalIndexer(myRank,numProc));
 
   // setup factory
   panzer::EpetraLinearObjFactory<short> la_factory(eComm,indexer);

   // build vectors from factory
   RCP<Thyra::MultiVectorBase<double> > ghostedVec = la_factory.getGhostedVector();
   RCP<Thyra::MultiVectorBase<double> > vec = la_factory.getVector();
   
   // convert to epetra vectors
   RCP<Epetra_Map> eMap = la_factory.getMap();
   RCP<Epetra_Map> eGhostedMap = la_factory.getGhostedMap();
   RCP<Epetra_MultiVector> ptrEVec = getEpetraMultiVector(vec,*eMap);
   RCP<Epetra_MultiVector> ptrEGhostedVec = getEpetraMultiVector(ghostedVec,*eGhostedMap);
   Epetra_Vector & eVec = *(*ptrEVec)(0);
   Epetra_Vector & eGhostedVec = *(*ptrEGhostedVec)(0);

   // check sizes of global epetra vectors
   TEST_EQUALITY(eVec.NumVectors(),1);
   TEST_EQUALITY(eVec.GlobalLength(),12);
   TEST_EQUALITY(eVec.MyLength(),6);

   TEST_EQUALITY(eGhostedVec.NumVectors(),1);
   TEST_EQUALITY(eGhostedVec.GlobalLength(),16);
   TEST_EQUALITY(eGhostedVec.MyLength(),8);

   // fill epetra vectors
   for(int i=0;i<eVec.MyLength();i++) 
      eVec[i] = double(eMap->GID(i))+0.1;
   for(int i=0;i<eGhostedVec.MyLength();i++) 
      eGhostedVec[i] = double(eGhostedMap->GID(i))+0.1;

   // run parallel assembly for global vector
   eVec.PutScalar(-10000.0);
   la_factory.ghostToGlobalVector(ghostedVec,vec);

   // check global vector 
   {
      for(int i=0;i<eVec.MyLength();i++) {
         int gid = eMap->GID(i); 
   
         if(gid==2 || gid==3 || gid==4 || gid==5) {
            TEST_FLOATING_EQUALITY(eVec[i],2.0*(double(gid)+0.1),1e-14);
         }
         else {
            TEST_FLOATING_EQUALITY(eVec[i],double(gid)+0.1,1e-14);
         }
      }
   }

   // construct ghosted vector
   eGhostedVec.PutScalar(-10000.0);
   la_factory.globalToGhostVector(vec,ghostedVec);

   // check ghosted vector 
   {
      eVec.PutScalar(0.0);
      for(int i=0;i<eGhostedVec.MyLength();i++) {
         int gid = eGhostedMap->GID(i); 
   
         if(gid==2 || gid==3 || gid==4 || gid==5) {
            TEST_FLOATING_EQUALITY(eGhostedVec[i],2.0*(double(gid)+0.1),1e-14);
         }
         else {
            TEST_FLOATING_EQUALITY(eGhostedVec[i],double(gid)+0.1,1e-14);
         }
      }
   }
}

}
