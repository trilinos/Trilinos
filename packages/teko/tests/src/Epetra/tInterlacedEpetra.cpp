// Thyra testing tools
#include "Thyra_TestingTools.hpp"
#include "Thyra_LinearOpTester.hpp"

// Thyra includes
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraOperatorWrapper.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"

// Teuchos includes
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

// TriUtils includes
#include "Trilinos_Util_CrsMatrixGallery.h"

#include "tInterlacedEpetra.hpp"

#include "Epetra/PB_InterlacedEpetra.hpp"


namespace PB {
namespace Test {

using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Thyra::VectorBase;
using Thyra::LinearOpBase;
using Thyra::createMember;
using Thyra::LinearOpTester;

void tInterlacedEpetra::initializeTest() 
{
}

int tInterlacedEpetra::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tInterlacedEpetra";

   status = test_buildSubMaps_num(verbosity,failstrm);
   allTests &= status;
   PB_TEST_MSG(stdstrm,1,"   \"buildSubMaps_num\" ... PASSED","   \"buildSubMaps_num\" ... FAILED");
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_buildSubMaps_vec(verbosity,failstrm);
   allTests &= status;
   PB_TEST_MSG(stdstrm,1,"   \"buildSubMaps_vec\" ... PASSED","   \"buildSubMaps_vec\" ... FAILED");
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      PB_TEST_MSG(failstrm,0,"tInterlacedEpetra...PASSED","tInterlacedEpetra...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      PB_TEST_MSG(failstrm,0,"...PASSED","tInterlacedEpetra...FAILED");
   }

   return failcount;
}

bool tInterlacedEpetra::test_buildSubMaps_num(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   const Epetra_Comm & comm = *GetComm();

   try {
      std::vector<std::pair<int,RCP<Epetra_Map> > > subMaps;
      int globals = 10;
      int numVars = 3;

      // build a set of submaps: this should fail
      PB::Epetra::buildSubMaps(globals,numVars,comm,subMaps); 

      TEST_ASSERT(false,
            "\n   tInerlacedEpetra::test_buildSubMaps_num: " << toString(status) << ": "
            "buildSubMaps(int,vector<pair<int,RCP<Epetra_Map> > >) did not throw "
            "with incorrect parameters");
   } 
   catch(...) {
      TEST_MSG(
            "\n   tInerlacedEpetra::test_buildSubMaps_num: "
         << "correctly threw an exception on incorrect parameters");
   }

   try {
      std::vector<std::pair<int,RCP<Epetra_Map> > > subMaps;
      int globals = 9;
      int numVars = 3;

      // build a set of submaps: this should fail
      PB::Epetra::buildSubMaps(globals,numVars,comm,subMaps); 

      TEST_EQUALITY(subMaps.size(),3,
            "\n   tInerlacedEpetra::test_buildSubMaps_num: " << toString(status) << ": "
         << "testing number of maps built ( " << subMaps.size() << " == " << 3 << "? ) ");

      bool cur = true;
      for(int i=0;i<3;i++)
         cur &= (subMaps[i].first==1);
      TEST_ASSERT(cur,
            "\n   tInerlacedEpetra::test_buildSubMaps_num: " << toString(status) << ": "
         << "testing that maps are associated with the correct numbers of variables");

      int * gids;

      // test the first of three used maps
      TEST_EQUALITY(subMaps[0].second->NumGlobalElements(),3,
            "\n   tInerlacedEpetra::test_buildSubMaps_num: " << toString(status) << ": "
         << "testing that first map has correct number of global elements ( "
         << subMaps[0].second->NumGlobalElements() << " == " << 3 << " ?)");
      gids = subMaps[0].second->MyGlobalElements();
      cur = (gids[0] == 0 && gids[1] == 3 && gids[2] == 6);
      TEST_ASSERT(cur,
            "\n   tInerlacedEpetra::test_buildSubMaps_num: " << toString(status) << ": "
         << "testing that first map is created correctly");

      // test the second of three used maps
      TEST_EQUALITY(subMaps[1].second->NumGlobalElements(),3,
            "\n   tInerlacedEpetra::test_buildSubMaps_num: " << toString(status) << ": "
         << "testing that second map has correct number of global elements ( "
         << subMaps[1].second->NumGlobalElements() << " == " << 3 << " ?)");
      gids = subMaps[1].second->MyGlobalElements();
      cur = (gids[0] == 1 && gids[1] == 4 && gids[2] == 7);
      TEST_ASSERT(cur,
            "\n   tInerlacedEpetra::test_buildSubMaps_num: " << toString(status) << ": "
         << "testing that second map is created correctly");

      // test the first of three used maps
      TEST_EQUALITY(subMaps[2].second->NumGlobalElements(),3,
            "\n   tInerlacedEpetra::test_buildSubMaps_num: " << toString(status) << ": "
         << "testing that the third map has correct number of global elements ( "
         << subMaps[2].second->NumGlobalElements() << " == " << 3 << " ?)");
      gids = subMaps[2].second->MyGlobalElements();
      cur = (gids[0] == 2 && gids[1] == 5 && gids[2] == 8);
      TEST_ASSERT(cur,
            "\n   tInerlacedEpetra::test_buildSubMaps_num: " << toString(status) << ": "
         << "testing that third map is created correctly");
   } 
   catch(...) {
      TEST_ASSERT(false,
            "\n   tInerlacedEpetra::test_buildSubMaps_num: " << toString(status) << ": "
            "buildSubMaps(int,vector<pair<int,RCP<Epetra_Map> > >) threw an unexpected "
            "exception");
   }

   return allPassed;
}

bool tInterlacedEpetra::test_buildSubMaps_vec(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   const Epetra_Comm & comm = *GetComm();

   try {
      std::vector<std::pair<int,RCP<Epetra_Map> > > subMaps;
      int globals = 15;

      std::vector<int> vars(3);
      vars[0] = 2;    
      vars[1] = 1;    
      vars[2] = 3;    

      // build a set of submaps: this should fail
      PB::Epetra::buildSubMaps(globals,vars,comm,subMaps); 

      TEST_ASSERT(false,
            "\n   tInerlacedEpetra::test_buildSubMaps_vec: " << toString(status) << ": "
            "buildSubMaps(int,vector<pair<int,RCP<Epetra_Map> > >) did not throw "
            "with incorrect parameters");
   } 
   catch(...) {
      TEST_MSG(
            "\n   tInerlacedEpetra::test_buildSubMaps_vec: "
         << "correctly threw an exception on incorrect parameters");
   }

   try {
      std::vector<std::pair<int,RCP<Epetra_Map> > > subMaps;
      int globals = 18;

      std::vector<int> vars(3);
      vars[0] = 2;    
      vars[1] = 1;    
      vars[2] = 3;    

      // build a set of submaps: this should fail
      PB::Epetra::buildSubMaps(globals,vars,comm,subMaps); 

      TEST_EQUALITY(subMaps.size(),3,
            "\n   tInerlacedEpetra::test_buildSubMaps_vec: " << toString(status) << ": "
         << "testing number of maps built ( " << subMaps.size() << " == " << 3 << "? ) ");

      bool cur = true;
      for(int i=0;i<3;i++)
         cur &= (subMaps[i].first==vars[i]);
      TEST_ASSERT(cur,
            "\n   tInerlacedEpetra::test_buildSubMaps_vec: " << toString(status) << ": "
         << "testing that maps are associated with the correct numbers of variables");

      int * gids;

      // test the first of three used maps
      TEST_EQUALITY(subMaps[0].second->NumGlobalElements(),6,
            "\n   tInerlacedEpetra::test_buildSubMaps_vec: " << toString(status) << ": "
         << "testing that first map has correct number of global elements ( "
         << subMaps[0].second->NumGlobalElements() << " == " << 6 << " ?)");
      gids = subMaps[0].second->MyGlobalElements();
      cur = (gids[0] == 0 && gids[1] == 1 
          && gids[2] == 6 && gids[3] == 7
          && gids[4] ==12 && gids[5] ==13);
      TEST_ASSERT(cur,
            "\n   tInerlacedEpetra::test_buildSubMaps_vec: " << toString(status) << ": "
         << "testing that first map is created correctly");

      // test the second of three used maps
      TEST_EQUALITY(subMaps[1].second->NumGlobalElements(),3,
            "\n   tInerlacedEpetra::test_buildSubMaps_vec: " << toString(status) << ": "
         << "testing that second map has correct number of global elements ( "
         << subMaps[1].second->NumGlobalElements() << " == " << 3 << " ?)");
      gids = subMaps[1].second->MyGlobalElements();
      cur = (gids[0] == 2 && gids[1] == 8 && gids[2] == 14);
      TEST_ASSERT(cur,
            "\n   tInerlacedEpetra::test_buildSubMaps_vec: " << toString(status) << ": "
         << "testing that second map is created correctly");

      // test the first of three used maps
      TEST_EQUALITY(subMaps[2].second->NumGlobalElements(),9,
            "\n   tInerlacedEpetra::test_buildSubMaps_vec: " << toString(status) << ": "
         << "testing that the third map has correct number of global elements ( "
         << subMaps[2].second->NumGlobalElements() << " == " << 9 << " ?)");
      gids = subMaps[2].second->MyGlobalElements();
      cur = (gids[0] == 3 && gids[1] == 4 && gids[2]== 5
          && gids[3] == 9 && gids[4] ==10 && gids[5]==11
          && gids[6] ==15 && gids[7] ==16 && gids[8]==17);
      TEST_ASSERT(cur,
            "\n   tInerlacedEpetra::test_buildSubMaps_vec: " << toString(status) << ": "
         << "testing that third map is created correctly");
   } 
   catch(...) {
      TEST_ASSERT(false,
            "\n   tInerlacedEpetra::test_buildSubMaps_vec: " << toString(status) << ": "
            "buildSubMaps(int,vector<pair<int,RCP<Epetra_Map> > >) threw an unexpected "
            "exception");
   }

   return allPassed;
}

} // end Test namespace
} // end PB namespace
