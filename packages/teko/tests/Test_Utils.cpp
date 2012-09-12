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

#include "Test_Utils.hpp"

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_OperatorVectorTypes.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"

#include "Teko_Utilities.hpp"

// Epetra includes
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"

#include <iostream>

using namespace Teuchos;

namespace Teko {
namespace Test {

const RCP<const Thyra::LinearOpBase<double> > build2x2(const Epetra_Comm & comm,double a,double b,double c,double d)
{
   RCP<Epetra_Map> map = rcp(new Epetra_Map(2,0,comm));

   int indicies[2];
   double row0[2];
   double row1[2];

   indicies[0] = 0;
   indicies[1] = 1;

   // build a CrsMatrix
   RCP<Epetra_CrsMatrix> blk  = rcp(new Epetra_CrsMatrix(Copy,*map,2));
   row0[0] = a; row0[1] = b;  // do a transpose here!
   row1[0] = c; row1[1] = d;
   blk->InsertGlobalValues(0,2,&row0[0],&indicies[0]);
   blk->InsertGlobalValues(1,2,&row1[0],&indicies[0]);
   blk->FillComplete();

   return Thyra::epetraLinearOp(blk);
}

const RCP<const Thyra::MultiVectorBase<double> > BlockVector(const Epetra_Vector & eu, const Epetra_Vector & ev,
        const RCP<const Thyra::VectorSpaceBase<double> > & vs)
{
   typedef RCP<const Thyra::MultiVectorBase<double> > Vector;

   const RCP<const Thyra::ProductVectorSpaceBase<double> > pvs 
         = rcp_dynamic_cast<const Thyra::ProductVectorSpaceBase<double> >(vs);

   RCP<const Epetra_MultiVector> meu = rcpFromRef(eu);
   RCP<const Epetra_MultiVector> mev = rcpFromRef(ev);
   const Vector u = Thyra::create_MultiVector(meu,pvs->getBlock(0)); 
   const Vector v = Thyra::create_MultiVector(mev,pvs->getBlock(1)); 

   // build rhs: this is ugly...in 2 steps
   // (i). allocate space for rhs "Product" vector, this is the range of A
   // (j). Need to do a dynamic cast to set the blocks (VectorBase doesn't know about setBlock)
   // const RCP<Thyra::MultiVectorBase<double> > rhs = Thyra::createMembers(vs,1);
   // rcp_dynamic_cast<Thyra::DefaultProductMultiVector<double> >(rhs)->setBlock(0,u);
   // rcp_dynamic_cast<Thyra::DefaultProductMultiVector<double> >(rhs)->setBlock(1,v);
   std::vector<RCP<Thyra::MultiVectorBase<double> > > blocks;
   blocks.push_back(rcp_const_cast<Thyra::MultiVectorBase<double> >(u));
   blocks.push_back(rcp_const_cast<Thyra::MultiVectorBase<double> >(v));

   return buildBlockedMultiVector(blocks);
}

void Print(std::ostream & os,const std::string & s,const RCP<const Thyra::MultiVectorBase<double> > & v)
{
   Thyra::ConstDetachedVectorView<double> view(v->col(0));

   os << s << " = "; 
   for(int i=0;i<v->range()->dim();i++) os << view[i] << " ";
   os << std::endl;
}

RCP<Thyra::VectorBase<double> > BuildIVector(int j,const RCP<const Thyra::VectorSpaceBase<double> > & vs)
{
   RCP<Thyra::VectorBase<double> > v = Thyra::createMember(vs);
   Thyra::assign<double>(v.ptr(), 0.0);

   {
     Thyra::DetachedVectorView<double> view(v);
     view[j] = 1.0;
   }

   return v;
}

void HardExtract(std::ostream & os,const RCP<const Thyra::LinearOpBase<double> > & A)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   int ncols = A->domain()->dim(); 
   int nrows = A->range()->dim(); 

   // print the matrix the hard way
   for(int i=0;i<nrows;i++) {
      const RCP<Thyra::VectorBase<double> > u = Teko::Test::BuildIVector(i,A->range());
      for(int j=0;j<ncols;j++) {
         const RCP<const Thyra::VectorBase<double> > ej = Teko::Test::BuildIVector(j,A->domain());
         const RCP<Thyra::VectorBase<double> > v = Thyra::createMember(A->range()); 
         ::Thyra::apply(*A, Thyra::NOTRANS, *ej, v.ptr());
          
         // this for some reason messes up valgrind (2/9/09)
         double value = A->range()->scalarProd(*u,*v);

         os << value << "  "; 
      } 
      os << std::endl;
   }
}

// add two Thyra vectors
const Teuchos::RCP<Thyra::MultiVectorBase<double> > Add(const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & x,
                                                   const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & y)
{
   return Add(1.0,x,1.0,y);
}

const Teuchos::RCP<Thyra::MultiVectorBase<double> > Add(double ax,const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & x,
                                                   double ay,const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & y)
{
   const Teuchos::RCP<Thyra::MultiVectorBase<double> > z = Thyra::createMembers(x->range(),1);

   Thyra::linear_combination<double>(Teuchos::tuple<double>(ax,ay),
         Teuchos::tuple(Teuchos::ptrInArg(*x),Teuchos::ptrInArg(*y)),0.0,z.ptr());
 
   return z;
}

// compute ||x-y||_2
double Difference(const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & x,
                  const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & y)
{
   Teuchos::RCP<const Thyra::MultiVectorBase<double> > z = Add(1.0,x,-1.0,y);

   return Thyra::norm_2(*z->col(0));
}

// construct a diagonal matrix
const Teuchos::RCP<const Thyra::LinearOpBase<double> > DiagMatrix(int cnt,double * vec,std::string label)
{
   const RCP<Epetra_SerialComm> comm = rcp(new Epetra_SerialComm());
   const RCP<Epetra_Map> map = rcp(new Epetra_Map(cnt,0,*comm));
   const RCP<Epetra_CrsMatrix> ptrF  = rcp(new Epetra_CrsMatrix(Copy,*map,1));

   // construct a diagonal matrix
   for(int i=0;i<cnt;i++)
      ptrF->InsertGlobalValues(i,1,&vec[i],&i);
   ptrF->FillComplete();

   // return thyra object
   return Thyra::epetraLinearOp(ptrF,label);
}

// declare static allocation
std::list<std::pair<Teuchos::RCP<UnitTest>,std::string> > UnitTest::testList;
Teuchos::RCP<const Epetra_Comm> UnitTest::comm;

void UnitTest::AddTest(const Teuchos::RCP<UnitTest> & ut,const std::string & name)
{
   // add a unit test and string to the list
   testList.push_back(std::make_pair(ut,name));
}

bool UnitTest::RunTests(int verbosity, std::ostream & stdstrm,std::ostream & failstrm)
{
   bool allPassed = true;
   int count=0,numfailed=0;
   std::list<std::pair<Teuchos::RCP<UnitTest>,std::string> >::iterator itr;

   bool isParallel = GetComm()->NumProc()>1;

   // loop over the tests and run each
   for(itr=testList.begin();itr!=testList.end();++itr) {
      int localrun = 0;
      int localfail = 0;

      // skip a test if its not parallel
      if(isParallel && not (itr->first)->isParallel()) continue;

      stdstrm << "Running test \"" << itr->second << (verbosity>=1 ? "\"\n" : "\" ... ");
      
      // run the tests
      (itr->first)->initializeTest();
      bool status = 0==(localfail = (itr->first)->runTest(verbosity,stdstrm,failstrm,localrun));

      // output some stuff to the standard stream
      if(verbosity>=1) {
         stdstrm << "Test \"" << itr->second << "\" completed ... ";
         if(status)
            stdstrm << Teko::Test::toString(status) << " (" << localrun << ")" << std::endl;
         else
            stdstrm << Teko::Test::toString(status) << " (" << localfail << ")" << std::endl;
      } 
       
      allPassed &= status;
      if(not status) numfailed+=localfail;
      count += localrun;
   }

   // output status
   stdstrm << std::endl;
   stdstrm << "Tests Passed: " << count-numfailed 
           << ", Tests Failed: " << numfailed << std::endl;
   stdstrm << "(Incidently, you want no failures)" << std::endl;

   return allPassed;
}

Teuchos::RCP<const Epetra_Comm> UnitTest::GetComm()
{
   return comm;
}

void UnitTest::SetComm(const Teuchos::RCP<const Epetra_Comm> & c)
{
   comm = c;
}

bool UnitTest::CheckParallelBools(bool myBool,int & failPID)
{
   int myInt = myBool ? 1 : 0;
   std::vector<int> bools(GetComm()->NumProc());

   GetComm()->GatherAll(&myInt,&bools[0],1);
 
   failPID = -1;

   bool result = true;
   for(unsigned int i=0;i<bools.size();i++) {
      result &= bools[i]==1 ? true : false;
      if(bools[i]!=1) failPID = i; 
   }

   return result;
}

} // end Tests
} // end Teko
