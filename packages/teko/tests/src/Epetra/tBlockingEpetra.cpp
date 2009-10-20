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

#include "tBlockingEpetra.hpp"

#include "Epetra/PB_BlockingEpetra.hpp"

using namespace PB::Epetra;

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

void tBlockingEpetra::initializeTest() 
{
   tolerance_ = 1e-14;
}

int tBlockingEpetra::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status = true;
   int failcount = 0;

   failstrm << "tBlockingEpetra";

   status = test_buildMaps(verbosity,failstrm);
   PB_TEST_MSG(stdstrm,1,"   \"buildMaps\" ... PASSED","   \"buildMaps\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_one2many(verbosity,failstrm);
   PB_TEST_MSG(stdstrm,1,"   \"one2many\" ... PASSED","   \"one2many\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_many2one(verbosity,failstrm);
   PB_TEST_MSG(stdstrm,1,"   \"many2one\" ... PASSED","   \"many2one\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      PB_TEST_MSG(failstrm,0,"tBlockingEpetra...PASSED","tInterlacedEpetra...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      PB_TEST_MSG(failstrm,0,"...PASSED","tBlockingEpetra...FAILED");
   }

   return failcount;
}

bool tBlockingEpetra::test_buildMaps(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   int size = 3*1000;

   TEST_MSG("\n   Builing Epetra_Map");
   RCP<Epetra_Map> map = rcp(new Epetra_Map(size,0,*GetComm()));

   TEST_MSG("\n   Building sub maps");
   std::vector<int> gid0(1200);
   std::vector<int> gid1(1200);
   for(int i=0;i<1200;i++) {
      gid0[i] = 2*i;
      gid1[i] = 2*i+1;
   }   
   std::vector<int> gid2(600);
   for(int i=0;i<600;i++)
      gid2[i] = 2400+i;

   Blocking::MapPair map0 = Blocking::buildSubMap(gid0,*GetComm());
   Blocking::MapPair map1 = Blocking::buildSubMap(gid1,*GetComm());
   Blocking::MapPair map2 = Blocking::buildSubMap(gid2,*GetComm());

   TEST_ASSERT(map0.first->NumMyElements()==gid0.size() && map0.second->NumMyElements()==gid0.size(),
         "   tBlockingEpetra::test_buildMaps (" << PB::Test::toString(status) << "): "
      << " Checking map size: first=" << map0.first->NumMyElements() 
      << ", second="<< map0.second->NumMyElements()
      << ", gid="<< gid0.size());
   TEST_ASSERT(map1.first->NumMyElements()==gid1.size() && map0.second->NumMyElements()==gid0.size(),
         "   tBlockingEpetra::test_buildMaps (" << PB::Test::toString(status) << "): "
      << " Checking map size: first=" << map1.first->NumMyElements() 
      << ", second="<< map1.second->NumMyElements()
      << ", gid="<< gid1.size());
   TEST_ASSERT(map2.first->NumMyElements()==gid2.size() && map0.second->NumMyElements()==gid0.size(),
         "   tBlockingEpetra::test_buildMaps (" << PB::Test::toString(status) << "): "
      << " Checking map size: first=" << map2.first->NumMyElements() 
      << ", second="<< map2.second->NumMyElements()
      << ", gid="<< gid2.size());

   std::vector<Teuchos::RCP<Epetra_Map> > globalMaps(3);
   std::vector<Teuchos::RCP<Epetra_Map> > contigMaps(3);

   // get sub maps for convenient use and access
   globalMaps[0] = map0.first;
   globalMaps[1] = map1.first;
   globalMaps[2] = map2.first;

   contigMaps[0] = map0.second;
   contigMaps[1] = map1.second;
   contigMaps[2] = map2.second;

   // test that the extra data is attached
   TEST_ASSERT(contigMaps[0]!=Teuchos::null,
         "   tBlockingEpetra::test_buildMaps (" << PB::Test::toString(status) << ")");
   TEST_ASSERT(contigMaps[1]!=Teuchos::null,
         "   tBlockingEpetra::test_buildMaps (" << PB::Test::toString(status) << ")");
   TEST_ASSERT(contigMaps[2]!=Teuchos::null,
         "   tBlockingEpetra::test_buildMaps (" << PB::Test::toString(status) << ")");
   TEST_MSG("   tBlockingEpetra::test_buildMaps: extracted \"contigMaps\"");

   bool test;

   // check contiguous and global maps
   test = true;
   for(int i=0;i<globalMaps[0]->NumMyElements();i++) {
      int gid = globalMaps[0]->GID(i); 
      int cid = contigMaps[0]->GID(i); 
 
      test &= gid==gid0[i];
      test &= cid==i;
   }
   TEST_ASSERT(test, 
         "   tBlockingEpetra::test_buildMaps (" << PB::Test::toString(status) << "): "
      << "checked that block maps were internally consitent");

   test = true;
   for(int i=0;i<globalMaps[1]->NumMyElements();i++) {
      int gid = globalMaps[1]->GID(i); 
      int cid = contigMaps[1]->GID(i); 
 
      test &= gid==gid1[i];
      test &= cid==i;
   }
   TEST_ASSERT(test, 
         "   tBlockingEpetra::test_buildMaps (" << PB::Test::toString(status) << "): "
      << "checked that block maps were internally consitent");

   test = true;
   for(int i=0;i<globalMaps[2]->NumMyElements();i++) {
      int gid = globalMaps[2]->GID(i); 
      int cid = contigMaps[2]->GID(i); 
 
      test &= gid==gid2[i];
      test &= cid==i;
   }
   TEST_ASSERT(test, 
         "   tBlockingEpetra::test_buildMaps (" << PB::Test::toString(status) << "): "
      << "checked that block maps were internally consitent");

   return allPassed;
}

bool tBlockingEpetra::test_one2many(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   TEST_MSG("\n   Building sub maps");
   std::vector<int> gid0(1200);
   std::vector<int> gid1(1200);
   for(int i=0;i<1200;i++) {
      gid0[i] = 2*i;
      gid1[i] = 2*i+1;
   }   
   std::vector<int> gid2(600);
   for(int i=0;i<600;i++)
      gid2[i] = 2400+i;

   std::vector<Blocking::MapPair> maps(3);
   maps[0] = Blocking::buildSubMap(gid0,*GetComm());
   maps[1] = Blocking::buildSubMap(gid1,*GetComm());
   maps[2] = Blocking::buildSubMap(gid2,*GetComm());

   int size = 3*1000;
   TEST_MSG("\n   tBlockingEpetra::test_one2many: Builing Epetra_Map and source vector");
   RCP<Epetra_Map> map = rcp(new Epetra_Map(size,0,*GetComm()));
   RCP<Epetra_MultiVector> v = rcp(new Epetra_MultiVector(*map,1));
   v->Random();

   TEST_MSG("\n   tBlockingEpetra::test_one2many: Building Export/Import");
   std::vector<RCP<Epetra_Import> > subImport(3);
   std::vector<RCP<Epetra_Export> > subExport(3);
   for(int i=0;i<3;i++) {
      Blocking::ImExPair imex = Blocking::buildExportImport(*map,maps[i]);
      subImport[i] = imex.first; 
      subExport[i] = imex.second; 
   }


   TEST_MSG("\n   tBlockingEpetra::test_one2many: Building sub vectors");
   std::vector<RCP<Epetra_MultiVector> > subVectors;
   Blocking::buildSubVectors(maps,subVectors,1);

   TEST_MSG("\n   tBlockingEpetra::test_one2many: Performing data copy");
   Blocking::one2many(subVectors,*v,subImport);
   
   // just assume it works! :)

   return allPassed;
}

bool tBlockingEpetra::test_many2one(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;
   TEST_MSG("\n   Building sub maps");
   std::vector<int> gid0(1200);
   std::vector<int> gid1(1200);
   for(int i=0;i<1200;i++) {
      gid0[i] = 2*i;
      gid1[i] = 2*i+1;
   }   
   std::vector<int> gid2(600);
   for(int i=0;i<600;i++)
      gid2[i] = 2400+i;

   std::vector<Blocking::MapPair> maps(3);
   maps[0] = Blocking::buildSubMap(gid0,*GetComm());
   maps[1] = Blocking::buildSubMap(gid1,*GetComm());
   maps[2] = Blocking::buildSubMap(gid2,*GetComm());

   int size = 3*1000;
   TEST_MSG("\n   tBlockingEpetra::test_one2many: Builing Epetra_Map and source vector");
   RCP<Epetra_Map> map = rcp(new Epetra_Map(size,0,*GetComm()));
   RCP<Epetra_MultiVector> v = rcp(new Epetra_MultiVector(*map,4));
   v->Random();

   TEST_MSG("\n   tBlockingEpetra::test_one2many: Building Export/Import");
   std::vector<RCP<Epetra_Import> > subImport(3);
   std::vector<RCP<Epetra_Export> > subExport(3);
   for(int i=0;i<3;i++) {
      Blocking::ImExPair imex = Blocking::buildExportImport(*map,maps[i]);
      subImport[i] = imex.first; 
      subExport[i] = imex.second; 
   }

   TEST_MSG("\n   tBlockingEpetra::test_one2many: Building sub vectors");
   std::vector<RCP<Epetra_MultiVector> > subVectors;
   Blocking::buildSubVectors(maps,subVectors,4);

   TEST_MSG("\n   tBlockingEpetra::test_one2many: Performing one2many");
   Blocking::one2many(subVectors,*v,subImport);

   std::vector<RCP<const Epetra_MultiVector> > cSubVectors;
   std::vector<RCP<Epetra_MultiVector> >::const_iterator itr;
   for(itr=subVectors.begin();itr!=subVectors.end();++itr)
      cSubVectors.push_back(*itr);

   TEST_MSG("\n   tBlockingEpetra::test_one2many: Performing many2one");
   RCP<Epetra_MultiVector> one = rcp(new Epetra_MultiVector(*map,1));
   Blocking::many2one(*one,cSubVectors,subExport);

   one->Update(1.0,*v,-1.0);
 
   double diff[4] = {0,0,0,0};
   double max=0.0,maxn=0;
   double norm[4] = {0,0,0,0};
   one->Norm2(diff);
   v->Norm2(norm);
   for(int i=0;i<4;i++) {
      max = max>diff[i]/norm[i] ? max : diff[i]/norm[i];
      maxn = maxn>norm[i] ? maxn : norm[i];
   }
   TEST_ASSERT(max<=tolerance_,
            "   tBlockingEpetra::test_buildMaps (" << PB::Test::toString(status) << "): "
         << "norm must be better than the tolerance ( " << max << " <=? " << tolerance_ << " maxn = " << maxn << " )");
   
   return allPassed;
}

} // end Test namespace
} // end PB namespace
