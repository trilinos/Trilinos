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

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_VerboseObject.hpp"

// Thyra includes
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"
#include "Thyra_SpmdMultiVectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

// include basic Tpetra information
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"

#include "Thyra_TpetraThyraWrappers.hpp"

#include "Teko_TpetraThyraConverter.hpp"

#include <iostream>
#include "tTpetraThyraConverter.hpp"

typedef Teko::ST ST;
typedef Teko::LO LO;
typedef Teko::GO GO;
typedef Teko::NT NT;

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::rcp_dynamic_cast;

namespace {

/*
double compareTpetraMVToThyra(const Tpetra_MultiVector & eX,
                            const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & tX,
                            int indexStart=-1,int indexEnd=-1); */
double compareTpetraMVToThyra(const Tpetra::MultiVector<ST,LO,GO,NT> & eX,
                            const Teuchos::RCP<const Thyra::MultiVectorBase<ST> > & tX,
                            int verbosity,std::ostream & os,GO indexStart=-1,GO indexEnd=-1);

double compareTpetraMVToThyra(const Tpetra::MultiVector<ST,LO,GO,NT> & eX,
                            const Teuchos::RCP<const Thyra::MultiVectorBase<ST> > & tX,
                            int verbosity,std::ostream & os,GO indexStart,GO indexEnd)
{
   using Teuchos::outArg;
   if(indexStart<0) {
      indexStart = 0;
      indexEnd = eX.getGlobalLength();
   }

   ST maxerr = 0.0;

   // check the base case
   const RCP<const Thyra::ProductMultiVectorBase<ST> > prodX
         = rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<ST> > (tX);
   if(prodX==Teuchos::null) {
      // base case
      TEST_MSG("      compareTpetraMVToThyra - base case ( " << indexStart << ", " << indexEnd << " )" );

      const Tpetra::Map<LO,GO,NT> & map = *eX.getMap();
      int vecs = eX.getNumVectors();
/*
      // get vector view for comparing elements
      TEST_MSG("         " << "getting DetachedMultiVectorView");
      Thyra::ConstDetachedMultiVectorView<double> view(*tX);

      bool result = true;
      TEST_MSG("         " << "checking elements");
      for(int i=0;i<map.NumMyElements();i++) {
         int gid = map.GID(i);

         // this is not in the range of vector elements we are interested in
         if(gid<indexStart || gid>=indexEnd) continue;

         // these values should be exactly equal
         for(int j=0;j<vecs;j++) {
            bool local = view(gid-indexStart,j) == eX[j][i];
            result &= local;
            if(not local) { 
               double diff = std::fabs(view(gid-indexStart,j) - eX[j][i]);
               maxerr = maxerr > diff ? maxerr : diff;
            }
         }
      }
      TEST_MSG("         " << "check completed");

      TEST_MSG("      compareTpetraMVToThyra - finished base case");
*/
      const Teuchos::RCP<const Thyra::SpmdMultiVectorBase<ST> > spmd_tX
            = Teuchos::rcp_dynamic_cast<const Thyra::SpmdMultiVectorBase<ST> >(tX);

      Thyra::Ordinal stride = 0;
      Teuchos::ArrayRCP<const ST> localBuffer;
      spmd_tX->getLocalData(outArg(localBuffer), outArg(stride));

      TEST_MSG("         " << "stride = " << stride);
      TEST_MSG("         " << "checking elements");
      int thyraIndex = 0;      
      for(size_t i=0;i<map.getLocalNumElements();i++) {
         GO gid = map.getGlobalElement(i);

         // this is not in the range of vector elements we are interested in
         if(gid<indexStart || gid>=indexEnd) continue;

         // these values should be equal
         for(int j=0;j<vecs;j++) {
            ST diff = std::fabs(localBuffer[j*stride+thyraIndex]-eX.getData(j)[i]);
            maxerr = maxerr > diff ? maxerr : diff;
         }

         thyraIndex++;
      }
      TEST_MSG("         " << "check completed: maxerr = " << maxerr);
      TEST_MSG("      compareTpetraMVToThyra - finished base case");

      return maxerr;
   }

   const RCP<const Thyra::ProductVectorSpaceBase<ST> > prodVS = prodX->productSpace(); 
   TEST_MSG("      compareTpetraMVToThyra - recurse (" << indexStart << ", " << indexEnd << " )");

   // loop over each subblock, comparing the thyra to tpetra
   // bool result = true;
   for(int i=0;i<prodVS->numBlocks();i++) {
      GO size = prodVS->getBlock(i)->dim();

      // run comparison routine on relavant values
      ST val = compareTpetraMVToThyra(eX,prodX->getMultiVectorBlock(i),verbosity,os,indexStart,indexStart+size); 

      // shift starting index
      indexStart+= size;

      maxerr = maxerr > val ? maxerr : val;
   }

   TEST_MSG("      compareTpetraMVToThyra - finished recurse");
   return maxerr;
}

} // end namespace


namespace Teko {
namespace Test {

void tTpetraThyraConverter::initializeTest() {}

int tTpetraThyraConverter::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tTpetraThyraConverter";

   status = test_blockThyraToTpetra(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"blockThyraToTpetra\" ... PASSED","   \"blockThyraToTpetra\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_single_blockThyraToTpetra(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"single_blockThyraToTpetra\" ... PASSED","   \"single_blockThyraToTpetra\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_blockTpetraToThyra(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"blockTpetraToThyra\" ... PASSED","   \"blockTpetraToThyra\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_single_blockTpetraToThyra(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"single_blockTpetraToThyra\" ... PASSED","   \"single_blockTpetraToThyra\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      Teko_TEST_MSG(failstrm,0,"tTpetraThyraConverter...PASSED","tTpetraThyraConverter...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      Teko_TEST_MSG(failstrm,0,"...PASSED","tTpetraThyraConverter...FAILED");
   }

   return failcount;
}

bool tTpetraThyraConverter::test_blockThyraToTpetra(int verbosity,std::ostream & os)
{
   bool status;
   bool allPassed = true;

   const Teuchos::Comm<int> & Comm = *GetComm_tpetra();
   const RCP<const Teuchos::Comm<Teuchos::Ordinal> > tComm = Thyra::convertTpetraToThyraComm(rcpFromRef(Comm));

   // get process information
   int numProc = Comm.getSize();
   // int myPID   = Comm.MyPID();

   // how big is this vector
   LO myElmts = 1000;
   GO glElmts = myElmts*numProc;

   // build vector space
   const RCP<const Thyra::VectorSpaceBase<ST> > vs
         = Thyra::defaultSpmdVectorSpace<ST>(tComm,myElmts,glElmts); 
   const RCP<const Thyra::VectorSpaceBase<ST> > prodVS = Thyra::productVectorSpace(vs,2); 

   // from the vector space build an tpetra map
   TEST_MSG("\n   1. creating Map");
   const RCP<const Tpetra::Map<LO,GO,NT> > map = Teko::TpetraHelpers::thyraVSToTpetraMap(*prodVS,tComm);

   // create a vector
   const RCP<Thyra::MultiVectorBase<ST> > tX = Thyra::createMembers<ST>(prodVS,5);
   Thyra::randomize<ST>(-10.0,10.0,tX.ptr()); 

   TEST_MSG("   2. creating MultiVector");

   const RCP<Tpetra::MultiVector<ST,LO,GO,NT> > eX = rcp(new Tpetra::MultiVector<ST,LO,GO,NT>(map,5));
   TEST_MSG("   3. calling blockThyraToTpetra");
   Teko::TpetraHelpers::blockThyraToTpetra(tX,*eX);

   TEST_ASSERT(eX!=Teuchos::null,
            "\n   tTpetraThyraConverter::test_blockThyraToTpetra " << toString(status) 
         << ": blockThyraToTpetra returns not null");

   TEST_MSG("   4. comparing Tpetra to Thyra");
   ST result = compareTpetraMVToThyra(*eX,tX,verbosity,os);
   TEST_ASSERT(result==0.0,
            "\n   tTpetraThyraConverter::test_blockThyraToTpetra: " << toString(status) 
         << ": Tpetra MV is compared to Thyra MV (maxdiff = " << result << ")");

   return allPassed;
}

bool tTpetraThyraConverter::test_single_blockThyraToTpetra(int verbosity,std::ostream & os)
{
   bool status;
   bool allPassed = true;

   const Teuchos::Comm<int> & Comm = *GetComm_tpetra();
   const RCP<const Teuchos::Comm<Teuchos::Ordinal> > tComm = Thyra::convertTpetraToThyraComm(rcpFromRef(Comm));

   // get process information
   int numProc = Comm.getSize();
   // int myPID   = Comm.MyPID();

   // how big is this vector
   LO myElmts = 1000;
   GO glElmts = myElmts*numProc;

   // build vector space
   const RCP<const Thyra::VectorSpaceBase<ST> > vs
         = Thyra::defaultSpmdVectorSpace<ST>(tComm,myElmts,glElmts); 

   // from the vector space build an tpetra map
   const RCP<const Tpetra::Map<LO,GO,NT> > map = Teko::TpetraHelpers::thyraVSToTpetraMap(*vs,tComm);

   // create a vector
   const RCP<Thyra::MultiVectorBase<ST> > tX = Thyra::createMembers<ST>(vs,5);
   Thyra::randomize<ST>(-10.0,10.0,tX.ptr()); 

   const RCP<Tpetra::MultiVector<ST,LO,GO,NT> > eX = rcp(new Tpetra::MultiVector<ST,LO,GO,NT>(map,5));
   Teko::TpetraHelpers::blockThyraToTpetra(tX,*eX);

   TEST_ASSERT(eX!=Teuchos::null,
            "\n   tTpetraThyraConverter::test_single_blockThyraToTpetra: " << toString(status) 
         << ": blockThyraToTpetra returns not null");

   ST result = compareTpetraMVToThyra(*eX,tX,verbosity,os);
   TEST_ASSERT(result==0.0,
            "\n   tTpetraThyraConverter::test_single_blockThyraToTpetra: " << toString(status) 
         << ": Tpetra MV is compared to Thyra MV (maxdiff = " << result << ")");

   return allPassed;
}

bool tTpetraThyraConverter::test_blockTpetraToThyra(int verbosity,std::ostream & os)
{
   bool status;
   bool allPassed = true;

   const Teuchos::Comm<int> & Comm = *GetComm_tpetra();
   const RCP<const Teuchos::Comm<Teuchos::Ordinal> > tComm = Thyra::convertTpetraToThyraComm(rcpFromRef(Comm));
 
   // get process information
   int numProc = Comm.getSize();
   // int myPID   = Comm.MyPID();

   // how big is this vector
   LO myElmts = 1000;
   GO glElmts = myElmts*numProc;

   // build vector space
   const RCP<const Thyra::VectorSpaceBase<ST> > vs
         = Thyra::defaultSpmdVectorSpace<ST>(tComm,myElmts,glElmts); 
   const RCP<const Thyra::VectorSpaceBase<ST> > prodVS = Thyra::productVectorSpace(vs,2); 
 
   // from the vector space build an tpetra map
   const RCP<const Tpetra::Map<LO,GO,NT> > map = Teko::TpetraHelpers::thyraVSToTpetraMap(*prodVS,tComm);
   
   // build an tpetra multivector 
   Tpetra::MultiVector<ST,LO,GO,NT> eX(map,3);
   eX.randomize();

   // build a Thyra copy of this Tpetra_MultiVector
   const RCP<Thyra::MultiVectorBase<ST> >  tX = Thyra::createMembers(prodVS,eX.getNumVectors());
   Teko::TpetraHelpers::blockTpetraToThyra(eX,tX.ptr());

   ST result = compareTpetraMVToThyra(eX,tX,verbosity,os);
   TEST_ASSERT(result==0.0,
            "\n   tTpetraThyraConverter::test_blockTpetraToThyra: " << toString(status) 
         << ": Tpetra MV is compared to Thyra MV (maxdiff = " << result << ")");

   return allPassed;
}

bool tTpetraThyraConverter::test_single_blockTpetraToThyra(int verbosity, std::ostream & os)
{
   bool status;
   bool allPassed = true;

   const Teuchos::Comm<int> & Comm = *GetComm_tpetra();
   const RCP<const Teuchos::Comm<Teuchos::Ordinal> > tComm = Thyra::convertTpetraToThyraComm(rcpFromRef(Comm));
 
   // get process information
   int numProc = Comm.getSize();
   // int myPID   = Comm.MyPID();

   // how big is this vector
   LO myElmts = 1000;
   GO glElmts = myElmts*numProc;

   // build vector space
   const RCP<const Thyra::VectorSpaceBase<ST> > vs
         = Thyra::defaultSpmdVectorSpace<ST>(tComm,myElmts,glElmts); 
   const RCP<const Thyra::VectorSpaceBase<ST> > prodVS = vs;
 
   // from the vector space build an tpetra map
   const RCP<const Tpetra::Map<LO,GO,NT> > map = Teko::TpetraHelpers::thyraVSToTpetraMap(*prodVS,tComm);
   
   // build an tpetra multivector 
   int vecs = 10;
   Tpetra::MultiVector<ST,LO,GO,NT> eX(map,vecs);
   eX.randomize();

   // build a Thyra copy of this Tpetra::MultiVector
   const RCP<Thyra::MultiVectorBase<ST> >  tX = Thyra::createMembers(prodVS,eX.getNumVectors());
   Teko::TpetraHelpers::blockTpetraToThyra(eX,tX.ptr());

   TEST_ASSERT(tX!=Teuchos::null,
            "\n   tTpetraThyraConverter::test_single_blockTpetraToThyra: " << toString(status) 
         << ": blockTpetraToThyra returns not null");

   ST result = compareTpetraMVToThyra(eX,tX,verbosity,os);
   TEST_ASSERT(result==0.0,
            "\n   tTpetraThyraConverter::test_single_blockTpetraToThyra: " << toString(status) 
         << ": Tpetra MV is compared to Thyra MV (maxdiff = " << result << ")");

   return allPassed;
}

} // end Test namespace
} // end Teko namespace
