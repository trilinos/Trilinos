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

#include "Teko_EpetraThyraConverter.hpp"

#include <iostream>
#include "tEpetraThyraConverter.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::rcp_dynamic_cast;

namespace {

/*
double compareEpetraMVToThyra(const Epetra_MultiVector & eX,
                            const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & tX,
                            int indexStart=-1,int indexEnd=-1); */
double compareEpetraMVToThyra(const Epetra_MultiVector & eX,
                            const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & tX,
                            int verbosity,std::ostream & os,int indexStart=-1,int indexEnd=-1);

const RCP<const Thyra::VectorSpaceBase<double> > buildCompositeSpace(const Epetra_Comm & Comm)
{
   const RCP<const Teuchos::Comm<Teuchos::Ordinal> > tComm = Thyra::create_Comm(rcpFromRef(Comm));

   // get process information
   int numProc = Comm.NumProc();
   // int myPID   = Comm.MyPID();

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

double compareEpetraMVToThyra(const Epetra_MultiVector & eX,
                            const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & tX,
                            int verbosity,std::ostream & os,int indexStart,int indexEnd)
{
   if(indexStart<0) {
      indexStart = 0;
      indexEnd = eX.GlobalLength();
   }

   double maxerr = 0.0;

   // check the base case
   const RCP<const Thyra::ProductMultiVectorBase<double> > prodX
         = rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<double> > (tX);
   if(prodX==Teuchos::null) {
      // base case
      TEST_MSG("      compareEpetraMVToThyra - base case ( " << indexStart << ", " << indexEnd << " )" );

      const Epetra_BlockMap & map = eX.Map();
      int vecs = eX.NumVectors();
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

      TEST_MSG("      compareEpetraMVToThyra - finished base case");
*/
      const Teuchos::RCP<const Thyra::SpmdMultiVectorBase<double> > spmd_tX
            = Teuchos::rcp_dynamic_cast<const Thyra::SpmdMultiVectorBase<double> >(tX);
      Thyra::Ordinal stride = 0;
      const double * localBuffer = 0;
      spmd_tX->getLocalData(&localBuffer,&stride);

      TEST_MSG("         " << "stride = " << stride);
      TEST_MSG("         " << "checking elements");
      int thyraIndex = 0;      
      for(int i=0;i<map.NumMyElements();i++) {
         int gid = map.GID(i);

         // this is not in the range of vector elements we are interested in
         if(gid<indexStart || gid>=indexEnd) continue;

         // these values should be equal
         for(int j=0;j<vecs;j++) {
            double diff = std::fabs(localBuffer[j*stride+thyraIndex]-eX[j][i]);
            maxerr = maxerr > diff ? maxerr : diff;
         }

         thyraIndex++;
      }
      TEST_MSG("         " << "check completed: maxerr = " << maxerr);
      TEST_MSG("      compareEpetraMVToThyra - finished base case");

      return maxerr;
   }

   const RCP<const Thyra::ProductVectorSpaceBase<double> > prodVS = prodX->productSpace(); 
   TEST_MSG("      compareEpetraMVToThyra - recurse (" << indexStart << ", " << indexEnd << " )");

   // loop over each subblock, comparing the thyra to epetra
   // bool result = true;
   for(int i=0;i<prodVS->numBlocks();i++) {
      int size = prodVS->getBlock(i)->dim();

      // run comparison routine on relavant values
      double val = compareEpetraMVToThyra(eX,prodX->getMultiVectorBlock(i),verbosity,os,indexStart,indexStart+size); 

      // shift starting index
      indexStart+= size;

      maxerr = maxerr > val ? maxerr : val;
   }

   TEST_MSG("      compareEpetraMVToThyra - finished recurse");
   return maxerr;
}

} // end namespace


namespace Teko {
namespace Test {

void tEpetraThyraConverter::initializeTest() {}

int tEpetraThyraConverter::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tEpetraThyraConverter";

   status = test_blockThyraToEpetra(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"blockThyraToEpetra\" ... PASSED","   \"blockThyraToEpetra\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_single_blockThyraToEpetra(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"single_blockThyraToEpetra\" ... PASSED","   \"single_blockThyraToEpetra\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_blockEpetraToThyra(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"blockEpetraToThyra\" ... PASSED","   \"blockEpetraToThyra\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_single_blockEpetraToThyra(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"single_blockEpetraToThyra\" ... PASSED","   \"single_blockEpetraToThyra\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      Teko_TEST_MSG(failstrm,0,"tEpetraThyraConverter...PASSED","tEpetraThyraConverter...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      Teko_TEST_MSG(failstrm,0,"...PASSED","tEpetraThyraConverter...FAILED");
   }

   return failcount;
}

bool tEpetraThyraConverter::test_blockThyraToEpetra(int verbosity,std::ostream & os)
{
   bool status;
   bool allPassed = true;

   const Epetra_Comm & Comm = *GetComm();
   const RCP<const Teuchos::Comm<Teuchos::Ordinal> > tComm = Thyra::create_Comm(rcpFromRef(Comm));

   // get process information
   int numProc = Comm.NumProc();
   // int myPID   = Comm.MyPID();

   // how big is this vector
   int myElmts = 1000;
   int glElmts = myElmts*numProc;

   // build vector space
   const RCP<const Thyra::VectorSpaceBase<double> > vs
         = Thyra::defaultSpmdVectorSpace<double>(tComm,myElmts,glElmts); 
   const RCP<const Thyra::VectorSpaceBase<double> > prodVS = Thyra::productVectorSpace(vs,2); 

   // from the vector space build an epetra map
   TEST_MSG("\n   1. creating Map");
   const RCP<const Epetra_Map> map = Teko::Epetra::thyraVSToEpetraMap(*prodVS,rcpFromRef(Comm));

   // create a vector
   const RCP<Thyra::MultiVectorBase<double> > tX = Thyra::createMembers<double>(prodVS,5);
   Thyra::randomize<double>(-10.0,10.0,tX.ptr()); 

   TEST_MSG("   2. creating MultiVector");

   const RCP<Epetra_MultiVector> eX = rcp(new Epetra_MultiVector(*map,5));
   TEST_MSG("   3. calling blockThyraToEpetra");
   Teko::Epetra::blockThyraToEpetra(tX,*eX);

   TEST_ASSERT(eX!=Teuchos::null,
            "\n   tEpetraThyraConverter::test_blockThyraToEpetra " << toString(status) 
         << ": blockThyraToEpetra returns not null");

   TEST_MSG("   4. comparing Epetra to Thyra");
   double result = compareEpetraMVToThyra(*eX,tX,verbosity,os);
   TEST_ASSERT(result==0.0,
            "\n   tEpetraThyraConverter::test_blockThyraToEpetra: " << toString(status) 
         << ": Epetra MV is compared to Thyra MV (maxdiff = " << result << ")");

   return allPassed;
}

bool tEpetraThyraConverter::test_single_blockThyraToEpetra(int verbosity,std::ostream & os)
{
   bool status;
   bool allPassed = true;

   const Epetra_Comm & Comm = *GetComm();
   const RCP<const Teuchos::Comm<Teuchos::Ordinal> > tComm = Thyra::create_Comm(rcpFromRef(Comm));

   // get process information
   int numProc = Comm.NumProc();
   // int myPID   = Comm.MyPID();

   // how big is this vector
   int myElmts = 1000;
   int glElmts = myElmts*numProc;

   // build vector space
   const RCP<const Thyra::VectorSpaceBase<double> > vs
         = Thyra::defaultSpmdVectorSpace<double>(tComm,myElmts,glElmts); 

   // from the vector space build an epetra map
   const RCP<const Epetra_Map> map = Teko::Epetra::thyraVSToEpetraMap(*vs,rcpFromRef(Comm));

   // create a vector
   const RCP<Thyra::MultiVectorBase<double> > tX = Thyra::createMembers<double>(vs,5);
   Thyra::randomize<double>(-10.0,10.0,tX.ptr()); 

   const RCP<Epetra_MultiVector> eX = rcp(new Epetra_MultiVector(*map,5));
   Teko::Epetra::blockThyraToEpetra(tX,*eX);

   TEST_ASSERT(eX!=Teuchos::null,
            "\n   tEpetraThyraConverter::test_single_blockThyraToEpetra: " << toString(status) 
         << ": blockThyraToEpetra returns not null");

   double result = compareEpetraMVToThyra(*eX,tX,verbosity,os);
   TEST_ASSERT(result==0.0,
            "\n   tEpetraThyraConverter::test_single_blockThyraToEpetra: " << toString(status) 
         << ": Epetra MV is compared to Thyra MV (maxdiff = " << result << ")");

   return allPassed;
}

bool tEpetraThyraConverter::test_blockEpetraToThyra(int verbosity,std::ostream & os)
{
   bool status;
   bool allPassed = true;

   const Epetra_Comm & Comm = *GetComm();
   const RCP<const Teuchos::Comm<Teuchos::Ordinal> > tComm = Thyra::create_Comm(rcpFromRef(Comm));
 
   // get process information
   int numProc = Comm.NumProc();
   // int myPID   = Comm.MyPID();

   // how big is this vector
   int myElmts = 1000;
   int glElmts = myElmts*numProc;

   // build vector space
   const RCP<const Thyra::VectorSpaceBase<double> > vs
         = Thyra::defaultSpmdVectorSpace<double>(tComm,myElmts,glElmts); 
   const RCP<const Thyra::VectorSpaceBase<double> > prodVS = Thyra::productVectorSpace(vs,2); 
 
   // from the vector space build an epetra map
   const RCP<const Epetra_Map> map = Teko::Epetra::thyraVSToEpetraMap(*prodVS,rcpFromRef(Comm));
   
   // build an epetra multivector 
   Epetra_MultiVector eX(*map,3);
   eX.Random();

   // build a Thyra copy of this Epetra_MultiVector
   const RCP<Thyra::MultiVectorBase<double> >  tX = Thyra::createMembers(prodVS,eX.NumVectors());
   Teko::Epetra::blockEpetraToThyra(eX,tX.ptr());

   double result = compareEpetraMVToThyra(eX,tX,verbosity,os);
   TEST_ASSERT(result==0.0,
            "\n   tEpetraThyraConverter::test_blockEpetraToThyra: " << toString(status) 
         << ": Epetra MV is compared to Thyra MV (maxdiff = " << result << ")");

   return allPassed;
}

bool tEpetraThyraConverter::test_single_blockEpetraToThyra(int verbosity, std::ostream & os)
{
   bool status;
   bool allPassed = true;

   const Epetra_Comm & Comm = *GetComm();
   const RCP<const Teuchos::Comm<Teuchos::Ordinal> > tComm = Thyra::create_Comm(rcpFromRef(Comm));
 
   // get process information
   int numProc = Comm.NumProc();
   // int myPID   = Comm.MyPID();

   // how big is this vector
   int myElmts = 1000;
   int glElmts = myElmts*numProc;

   // build vector space
   const RCP<const Thyra::VectorSpaceBase<double> > vs
         = Thyra::defaultSpmdVectorSpace<double>(tComm,myElmts,glElmts); 
   const RCP<const Thyra::VectorSpaceBase<double> > prodVS = vs;
 
   // from the vector space build an epetra map
   const RCP<const Epetra_Map> map = Teko::Epetra::thyraVSToEpetraMap(*prodVS,rcpFromRef(Comm));
   
   // build an epetra multivector 
   int vecs = 10;
   Epetra_MultiVector eX(*map,vecs);
   eX.Random();

   // build a Thyra copy of this Epetra_MultiVector
   const RCP<Thyra::MultiVectorBase<double> >  tX = Thyra::createMembers(prodVS,eX.NumVectors());
   Teko::Epetra::blockEpetraToThyra(eX,tX.ptr());

   TEST_ASSERT(tX!=Teuchos::null,
            "\n   tEpetraThyraConverter::test_single_blockEpetraToThyra: " << toString(status) 
         << ": blockEpetraToThyra returns not null");

   double result = compareEpetraMVToThyra(eX,tX,verbosity,os);
   TEST_ASSERT(result==0.0,
            "\n   tEpetraThyraConverter::test_single_blockEpetraToThyra: " << toString(status) 
         << ": Epetra MV is compared to Thyra MV (maxdiff = " << result << ")");

   return allPassed;
}

} // end Test namespace
} // end Teko namespace
