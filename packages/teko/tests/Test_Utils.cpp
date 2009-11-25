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
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_DefaultLinearOpSourceDecl.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_SpmdLinearOpBase.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_OperatorVectorTypes.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"

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

const RCP<const Thyra::VectorBase<double> > BlockVector(const Epetra_Vector & eu, const Epetra_Vector & ev,
        const RCP<const Thyra::VectorSpaceBase<double> > & vs)
{
   typedef RCP<const Thyra::VectorBase<double> > Vector;

   const RCP<const Thyra::ProductVectorSpaceBase<double> > pvs 
         = rcp_dynamic_cast<const Thyra::ProductVectorSpaceBase<double> >(vs);
   const Vector u = Thyra::create_Vector(rcpFromRef(eu),pvs->getBlock(0)); 
   const Vector v = Thyra::create_Vector(rcpFromRef(ev),pvs->getBlock(1)); 

   // build rhs: this is ugly...in 2 steps
   // (i). allocate space for rhs "Product" vector, this is the range of A
   // (j). Need to do a dynamic cast to set the blocks (VectorBase doesn't know about setBlock)
   const RCP<Thyra::VectorBase<double> > rhs = Thyra::createMember(vs);
   rcp_dynamic_cast<Thyra::DefaultProductVector<double> >(rhs)->setBlock(0,u);
   rcp_dynamic_cast<Thyra::DefaultProductVector<double> >(rhs)->setBlock(1,v);

   return rhs;
}

void Print(std::ostream & os,const std::string & s,const RCP<const Thyra::VectorBase<double> > & v)
{
   Thyra::ConstDetachedVectorView<double> view(v);

   os << s << " = "; 
   for(int i=0;i<v->space()->dim();i++) os << view[i] << " ";
   os << std::endl;
}

RCP<Thyra::VectorBase<double> > BuildIVector(int j,const RCP<const Thyra::VectorSpaceBase<double> > & vs)
{
   RCP<Thyra::VectorBase<double> > v = Thyra::createMember(vs);
   Thyra::DetachedVectorView<double> view(v);
   Thyra::assign<double>(&*v,0.0);
   view[j] = 1.0;

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
         A->apply(Thyra::NONCONJ_ELE,*ej,&*v);
          
         // this for some reason messes up valgrind (2/9/09)
         double value = A->range()->scalarProd(*u,*v);

         os << value << "  "; 
      } 
      os << std::endl;
   }
}

// add two Thyra vectors
const Teuchos::RCP<Thyra::VectorBase<double> > Add(const Teuchos::RCP<const Thyra::VectorBase<double> > & x,
                                                   const Teuchos::RCP<const Thyra::VectorBase<double> > & y)
{
   return Add(1.0,x,1.0,y);
}

const Teuchos::RCP<Thyra::VectorBase<double> > Add(double ax,const Teuchos::RCP<const Thyra::VectorBase<double> > & x,
                                                   double ay,const Teuchos::RCP<const Thyra::VectorBase<double> > & y)
{
   const Teuchos::RCP<Thyra::VectorBase<double> > z = Thyra::createMember(x->space());

   Thyra::linear_combination<double>(Teuchos::tuple<double>(ax,ay),
         Teuchos::tuple(Teuchos::ptrInArg(*x),Teuchos::ptrInArg(*y)),0.0,z.ptr());
 
   return z;
}

// compute ||x-y||_2
double Difference(const Teuchos::RCP<const Thyra::VectorBase<double> > & x,
                  const Teuchos::RCP<const Thyra::VectorBase<double> > & y)
{
   Teuchos::RCP<const Thyra::VectorBase<double> > z = Add(1.0,x,-1.0,y);

   return Thyra::norm_2(*z);
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
