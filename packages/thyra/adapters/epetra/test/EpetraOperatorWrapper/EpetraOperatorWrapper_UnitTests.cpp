/*
// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#include "Thyra_EpetraOperatorWrapper.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Trilinos_Util_CrsMatrixGallery.h"
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

#include "Teuchos_UnitTestHarness.hpp"


namespace Thyra {


using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::inOutArg;
using Teuchos::as;


TEUCHOS_UNIT_TEST( EpetraOperatorWrapper, basic )
{

#ifdef HAVE_MPI
   Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
   Epetra_SerialComm comm;
#endif

   out << "\nRunning on " << comm.NumProc() << " processors\n";

   int nx = 39; // essentially random values
   int ny = 53;

   out << "Using Trilinos_Util to create test matrices\n";

   // create some big blocks to play with
   Trilinos_Util::CrsMatrixGallery FGallery("recirc_2d",comm,false); // CJ TODO FIXME: change for Epetra64
   FGallery.Set("nx",nx);
   FGallery.Set("ny",ny);
   RCP<Epetra_CrsMatrix> F = rcp(FGallery.GetMatrix(),false);

   Trilinos_Util::CrsMatrixGallery CGallery("laplace_2d",comm,false); // CJ TODO FIXME: change for Epetra64
   CGallery.Set("nx",nx);
   CGallery.Set("ny",ny);
   RCP<Epetra_CrsMatrix> C = rcp(CGallery.GetMatrix(),false);

   Trilinos_Util::CrsMatrixGallery BGallery("diag",comm,false); // CJ TODO FIXME: change for Epetra64
   BGallery.Set("nx",nx*ny);
   BGallery.Set("a",5.0);
   RCP<Epetra_CrsMatrix> B = rcp(BGallery.GetMatrix(),false);

   Trilinos_Util::CrsMatrixGallery BtGallery("diag",comm,false); // CJ TODO FIXME: change for Epetra64
   BtGallery.Set("nx",nx*ny);
   BtGallery.Set("a",3.0);
   RCP<Epetra_CrsMatrix> Bt = rcp(BtGallery.GetMatrix(),false);

   // load'em up in a thyra operator
   out << "Building block2x2 Thyra matrix ... wrapping in EpetraOperatorWrapper\n";
   const RCP<const LinearOpBase<double> > A =
     Thyra::block2x2<double>(
       Thyra::epetraLinearOp(F),
       Thyra::epetraLinearOp(Bt),
       Thyra::epetraLinearOp(B),
       Thyra::epetraLinearOp(C),
       "A"
       );

   const RCP<Thyra::EpetraOperatorWrapper> epetra_A =
     rcp(new Thyra::EpetraOperatorWrapper(A));

   // begin the tests!
   const Epetra_Map & rangeMap  = epetra_A->OperatorRangeMap();
   const Epetra_Map & domainMap = epetra_A->OperatorDomainMap();

   // check to see that the number of global elements is correct
   TEST_EQUALITY(rangeMap.NumGlobalElements(), 2*nx*ny);
   TEST_EQUALITY(domainMap.NumGlobalElements(), 2*nx*ny);

   // largest global ID should be one less then the # of elements
   TEST_EQUALITY(rangeMap.NumGlobalElements()-1, rangeMap.MaxAllGID());
   TEST_EQUALITY(domainMap.NumGlobalElements()-1, domainMap.MaxAllGID());

   // create a vector to test: copyThyraIntoEpetra
   {
      const RCP<VectorBase<double> > tv = Thyra::createMember(A->domain());
      Thyra::randomize(-100.0, 100.0, tv.ptr());
      const RCP<const VectorBase<double> > tv_0 =
        Thyra::productVectorBase<double>(tv)->getVectorBlock(0);
      const RCP<const VectorBase<double> > tv_1 =
        Thyra::productVectorBase<double>(tv)->getVectorBlock(1);
      const Thyra::ConstDetachedSpmdVectorView<double> vv_0(tv_0);
      const Thyra::ConstDetachedSpmdVectorView<double> vv_1(tv_1);

      int off_0 = vv_0.globalOffset();
      int off_1 = vv_1.globalOffset();
      
      // create its Epetra counter part
      Epetra_Vector ev(epetra_A->OperatorDomainMap());
      epetra_A->copyThyraIntoEpetra(*tv, ev);

      // compare handle_tv to ev!
      TEST_EQUALITY(tv->space()->dim(), as<Ordinal>(ev.GlobalLength()));
      const int numMyElements = domainMap.NumMyElements();
      double tval = 0.0;
      for(int i=0; i < numMyElements; i++) {
         int gid = domainMap.GID(i);
         if(gid<nx*ny)
            tval = vv_0[gid-off_0];
         else
            tval = vv_1[gid-off_1-nx*ny];
         TEST_EQUALITY(ev[i], tval);
      }
   }

   // create a vector to test: copyEpetraIntoThyra
   {
      // create an Epetra vector
     Epetra_Vector ev(epetra_A->OperatorDomainMap());
     ev.Random();

      // create its thyra counterpart
      const RCP<VectorBase<double> > tv = Thyra::createMember(A->domain());
      const RCP<const VectorBase<double> > tv_0 =
        Thyra::productVectorBase<double>(tv)->getVectorBlock(0);
      const RCP<const VectorBase<double> > tv_1 =
        Thyra::productVectorBase<double>(tv)->getVectorBlock(1);
      const Thyra::ConstDetachedSpmdVectorView<double> vv_0(tv_0);
      const Thyra::ConstDetachedSpmdVectorView<double> vv_1(tv_1);

      int off_0 = rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<double> >(
        tv_0->space())->localOffset();
      int off_1 = rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<double> >(
        tv_1->space())->localOffset();

      epetra_A->copyEpetraIntoThyra(ev, tv.ptr());
   
      // compare tv to ev!
      TEST_EQUALITY(tv->space()->dim(), as<Ordinal>(ev.GlobalLength()));
      int numMyElements = domainMap.NumMyElements();
      double tval = 0.0;
      for(int i=0;i<numMyElements;i++) {
         int gid = domainMap.GID(i);
         if(gid<nx*ny)
            tval = vv_0[gid-off_0];
         else
            tval = vv_1[gid-off_1-nx*ny];
         TEST_EQUALITY(ev[i], tval);
      }
   }

   // Test using Thyra::LinearOpTester
   const RCP<const LinearOpBase<double> > thyraEpetraOp = epetraLinearOp(epetra_A);
   LinearOpTester<double> linearOpTester;
   linearOpTester.show_all_tests(true);
   const bool checkResult = linearOpTester.check(*thyraEpetraOp, inOutArg(out));
   TEST_ASSERT(checkResult);

}


} // namespace Thyra
