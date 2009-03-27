#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_VerboseObject.hpp"

// Thyra includes
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"
#include "Thyra_SpmdMultiVectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

// include basic Epetra information
#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
   #include "mpi.h"
#else
   #include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"

#include "Thyra_EpetraThyraWrappers.hpp"

#include "Epetra/PB_EpetraThyraConverter.hpp"

#include <iostream>
#include "tEpetraThyraConverter.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::rcp_dynamic_cast;

namespace {

bool compareEpetraMVToThyra(const Epetra_MultiVector & eX,
                            const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & tX,
                            int indexStart=-1,int indexEnd=-1);

const RCP<const Thyra::VectorSpaceBase<double> > buildCompositeSpace(const Epetra_Comm & Comm)
{
   const RCP<const Teuchos::Comm<int> > tComm = Thyra::create_Comm(rcpFromRef(Comm));

   // get process information
   int numProc = Comm.NumProc();
   int myPID   = Comm.MyPID();

   // how big is this vector
   int uElmts = 4;
   int vElmts = 8;
   int pElmts = 6;

   // build vector space
   Teuchos::Array<RCP<const Thyra::VectorSpaceBase<double> > > uvArray;
   uvArray.push_back(Thyra::defaultSpmdVectorSpace<double>(tComm,uElmts,uElmts*numProc)); 
   uvArray.push_back(Thyra::defaultSpmdVectorSpace<double>(tComm,vElmts,vElmts*numProc)); 

   Teuchos::Array<RCP<const Thyra::VectorSpaceBase<double> > > uvpArray;
   uvpArray.push_back(Thyra::productVectorSpace<double>(uvArray));
   uvpArray.push_back(Thyra::defaultSpmdVectorSpace<double>(tComm,pElmts,pElmts*numProc)); 

   const RCP<const Thyra::VectorSpaceBase<double> > uvpVS = Thyra::productVectorSpace<double>(uvpArray);

   return uvpVS;
}

bool compareEpetraMVToThyra(const Epetra_MultiVector & eX,
                            const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & tX,
                            int indexStart,int indexEnd)
{
   if(indexStart<0) {
      indexStart = 0;
      indexEnd = eX.GlobalLength();
   }

   // check the base case
   const RCP<const Thyra::ProductMultiVectorBase<double> > prodX
         = rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<double> > (tX);
   if(prodX==Teuchos::null) {
      // base case

      const Epetra_BlockMap & map = eX.Map();
      int vecs = eX.NumVectors();
      int * indicies = map.MyGlobalElements();

      // get vector view for comparing elements
      Thyra::ConstDetachedMultiVectorView<double> view(*tX);

      bool result = true;
      for(int i=0;i<map.NumMyElements();i++) {
         int gid = map.GID(i);

         // this is not in the range of vector elements we are interested in
         if(gid<indexStart || gid>=indexEnd) continue;

         // these values should be exactly equal
         for(int j=0;j<vecs;j++)
            result &= view(gid-indexStart,j) == eX[j][i];
      }

      return result;
   }

   const RCP<const Thyra::ProductVectorSpaceBase<double> > prodVS = prodX->productSpace(); 

   // loop over each subblock, comparing the thyra to epetra
   bool result = true;
   for(int i=0;i<prodVS->numBlocks();i++) {
      int size = prodVS->getBlock(i)->dim();

      // run comparison routine on relavant values
      result &= compareEpetraMVToThyra(eX,prodX->getMultiVectorBlock(i),indexStart,indexStart+size); 

      // shift starting index
      indexStart+= size;
   }

   return result;
}

} // end namespace


namespace PB {
namespace Test {

void tEpetraThyraConverter::initializeTest() {}

int tEpetraThyraConverter::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tEpetraThyraConverter";

   status = test_blockThyraToEpetra(verbosity,failstrm);
   allTests &= status;
   PB_TEST_MSG(stdstrm,1,"   \"blockThyraToEpetra\" ... PASSED","   \"blockThyraToEpetra\" ... FAILED");
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_single_blockThyraToEpetra(verbosity,failstrm);
   allTests &= status;
   PB_TEST_MSG(stdstrm,1,"   \"single_blockThyraToEpetra\" ... PASSED","   \"single_blockThyraToEpetra\" ... FAILED");
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_blockEpetraToThyra(verbosity,failstrm);
   allTests &= status;
   PB_TEST_MSG(stdstrm,1,"   \"blockEpetraToThyra\" ... PASSED","   \"blockEpetraToThyra\" ... FAILED");
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_single_blockEpetraToThyra(verbosity,failstrm);
   allTests &= status;
   PB_TEST_MSG(stdstrm,1,"   \"single_blockEpetraToThyra\" ... PASSED","   \"single_blockEpetraToThyra\" ... FAILED");
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      PB_TEST_MSG(failstrm,0,"tEpetraThyraConverter...PASSED","tEpetraThyraConverter...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      PB_TEST_MSG(failstrm,0,"...PASSED","tEpetraThyraConverter...FAILED");
   }

   return failcount;
}

bool tEpetraThyraConverter::test_blockThyraToEpetra(int verbosity,std::ostream & os)
{
   bool status;
   bool allPassed = true;

   const Epetra_Comm & Comm = *GetComm();
   const RCP<const Teuchos::Comm<int> > tComm = Thyra::create_Comm(rcpFromRef(Comm));

   // get process information
   int numProc = Comm.NumProc();
   int myPID   = Comm.MyPID();

   // how big is this vector
   int myElmts = 1000;
   int glElmts = myElmts*numProc;

   // build vector space
   const RCP<const Thyra::VectorSpaceBase<double> > vs
         = Thyra::defaultSpmdVectorSpace<double>(tComm,myElmts,glElmts); 
   const RCP<const Thyra::VectorSpaceBase<double> > prodVS = Thyra::productVectorSpace(vs,2); 

   // from the vector space build an epetra map
   const RCP<const Epetra_Map> map = PB::Epetra::thyraVSToEpetraMap(*prodVS,rcpFromRef(Comm));

   // create a vector
   const RCP<Thyra::MultiVectorBase<double> > tX = Thyra::createMembers<double>(prodVS,5);
   Thyra::randomize<double>(-10.0,10.0,tX); 

   const RCP<Epetra_MultiVector> eX = rcp(new Epetra_MultiVector(*map,5));
   PB::Epetra::blockThyraToEpetra(tX,*eX);

   TEST_ASSERT(eX!=Teuchos::null,
            "   tEpetraThyraConverter::test_blockThyraToEpetra: " << toString(status) 
         << ": blockThyraToEpetra returns not null");

   bool result = compareEpetraMVToThyra(*eX,tX);
   TEST_ASSERT(result,
            "   tEpetraThyraConverter::test_blockThyraToEpetra: " << toString(status) 
         << ": Epetra MV is compared to Thyra MV");

   return allPassed;
}

bool tEpetraThyraConverter::test_single_blockThyraToEpetra(int verbosity,std::ostream & os)
{
   bool status;
   bool allPassed = true;

   const Epetra_Comm & Comm = *GetComm();
   const RCP<const Teuchos::Comm<int> > tComm = Thyra::create_Comm(rcpFromRef(Comm));

   // get process information
   int numProc = Comm.NumProc();
   int myPID   = Comm.MyPID();

   // how big is this vector
   int myElmts = 1000;
   int glElmts = myElmts*numProc;

   // build vector space
   const RCP<const Thyra::VectorSpaceBase<double> > vs
         = Thyra::defaultSpmdVectorSpace<double>(tComm,myElmts,glElmts); 

   // from the vector space build an epetra map
   const RCP<const Epetra_Map> map = PB::Epetra::thyraVSToEpetraMap(*vs,rcpFromRef(Comm));

   // create a vector
   const RCP<Thyra::MultiVectorBase<double> > tX = Thyra::createMembers<double>(vs,5);
   Thyra::randomize<double>(-10.0,10.0,tX); 

   const RCP<Epetra_MultiVector> eX = rcp(new Epetra_MultiVector(*map,5));
   PB::Epetra::blockThyraToEpetra(tX,*eX);


   TEST_ASSERT(eX!=Teuchos::null,
            "   tEpetraThyraConverter::test_single_blockThyraToEpetra: " << toString(status) 
         << ": blockThyraToEpetra returns not null");

   bool result = compareEpetraMVToThyra(*eX,tX);
   TEST_ASSERT(result,
            "   tEpetraThyraConverter::test_single_blockThyraToEpetra: " << toString(status) 
         << ": Epetra MV is compared to Thyra MV");

   return allPassed;
}

bool tEpetraThyraConverter::test_blockEpetraToThyra(int verbosity,std::ostream & os)
{
   bool status;
   bool allPassed = true;

   const Epetra_Comm & Comm = *GetComm();
   const RCP<const Teuchos::Comm<int> > tComm = Thyra::create_Comm(rcpFromRef(Comm));
 
   // get process information
   int numProc = Comm.NumProc();
   int myPID   = Comm.MyPID();

   // how big is this vector
   int myElmts = 10;
   int glElmts = myElmts*numProc;

   // build vector space
   const RCP<const Thyra::VectorSpaceBase<double> > vs
         = Thyra::defaultSpmdVectorSpace<double>(tComm,myElmts,glElmts); 
   const RCP<const Thyra::VectorSpaceBase<double> > prodVS = Thyra::productVectorSpace(vs,2); 
 
   // from the vector space build an epetra map
   const RCP<const Epetra_Map> map = PB::Epetra::thyraVSToEpetraMap(*prodVS,rcpFromRef(Comm));
   
   // build an epetra multivector 
   Epetra_MultiVector eX(*map,3);
   eX.Random();

   // build a Thyra copy of this Epetra_MultiVector
   const RCP<Thyra::MultiVectorBase<double> >  tX = Thyra::createMembers(prodVS,eX.NumVectors());
   PB::Epetra::blockEpetraToThyra(eX,tX.ptr());

   bool result = true;
   result &= compareEpetraMVToThyra(eX,tX);
   TEST_ASSERT(result,
            "   tEpetraThyraConverter::test_blockEpetraToThyra: " << toString(status) 
         << ": Epetra MV is compared to Thyra MV");

   return allPassed;
}

bool tEpetraThyraConverter::test_single_blockEpetraToThyra(int verbosity, std::ostream & os)
{
   bool status;
   bool allPassed = true;

   const Epetra_Comm & Comm = *GetComm();
   const RCP<const Teuchos::Comm<int> > tComm = Thyra::create_Comm(rcpFromRef(Comm));
 
   // get process information
   int numProc = Comm.NumProc();
   int myPID   = Comm.MyPID();

   // how big is this vector
   int myElmts = 1000;
   int glElmts = myElmts*numProc;

   // build vector space
   const RCP<const Thyra::VectorSpaceBase<double> > vs
         = Thyra::defaultSpmdVectorSpace<double>(tComm,myElmts,glElmts); 
   const RCP<const Thyra::VectorSpaceBase<double> > prodVS = vs;
 
   // from the vector space build an epetra map
   const RCP<const Epetra_Map> map = PB::Epetra::thyraVSToEpetraMap(*prodVS,rcpFromRef(Comm));
   
   // build an epetra multivector 
   int vecs = 10;
   Epetra_MultiVector eX(*map,vecs);
   eX.Random();

   // build a Thyra copy of this Epetra_MultiVector
   const RCP<Thyra::MultiVectorBase<double> >  tX = Thyra::createMembers(prodVS,eX.NumVectors());
   PB::Epetra::blockEpetraToThyra(eX,tX.ptr());

   TEST_ASSERT(tX!=Teuchos::null,
            "   tEpetraThyraConverter::test_single_blockEpetraToThyra: " << toString(status) 
         << ": blockEpetraToThyra returns not null");

   bool result = compareEpetraMVToThyra(eX,tX);
   TEST_ASSERT(result,
            "   tEpetraThyraConverter::test_single_blockEpetraToThyra: " << toString(status) 
         << ": Epetra MV is compared to Thyra MV");

   return allPassed;
}

} // end Test namespace
} // end PB namespace
