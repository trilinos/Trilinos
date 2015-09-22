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

#include "tBlockUpperTriInverseOp.hpp"

#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"

#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "Teko_Utilities.hpp"
#include "Teko_BlockUpperTriInverseOp.hpp"

namespace Teko {
namespace Test {

using Teuchos::RCP;
using Teuchos::rcp;

static const RCP<const Thyra::LinearOpBase<double> > build2x2(const RCP<const Epetra_Map> & map,double a,double b,double c,double d)
{
   int indicies[2];
   double row0[2];
   double row1[2];

   indicies[0] = 0;
   indicies[1] = 1;

   // build a CrsMatrix
   RCP<Epetra_CrsMatrix> blk  = rcp(new Epetra_CrsMatrix(Copy,*map,2));
   row0[0] = a; row0[1] =  b; 
   row1[0] = c; row1[1] = d; 
   blk->InsertGlobalValues(0,2,&row0[0],&indicies[0]);
   blk->InsertGlobalValues(1,2,&row1[0],&indicies[0]);
   blk->FillComplete();

   return Thyra::epetraLinearOp(blk);
}

void tBlockUpperTriInverseOp::initializeTest()
{
   std::vector<int> indicies(2);
   std::vector<double> row0(2),row1(2);

   RCP<const Epetra_Comm> comm = GetComm();
   const RCP<Epetra_Map> map = rcp(new Epetra_Map(2,0,*comm));

   tolerance_ = 1.0e-11;
   RCP<const Thyra::LinearOpBase<double> > blk;

   // build forward operator
   A_ = rcp(new Thyra::DefaultBlockedLinearOp<double>());
   A_->beginBlockFill(3,3);

   // build 0,0 matrix
   blk = build2x2(map,  1.0,  3.0,
                        2.0, -1.0);
   A_->setBlock(0,0,blk);

   // build 0,1 matrix
   blk = build2x2(map,  2.0,  9.0,
                        8.0,  3.0);
   A_->setBlock(0,1,blk);

   // build 1,1 matrix
   blk = build2x2(map,  7.0,  8.0,
                       -2.0,  4.0);
   A_->setBlock(1,1,blk);

   // build 1,2 matrix
   blk = build2x2(map, -1.0,  6.0,
                        2.0,  1.0);
   A_->setBlock(1,2,blk);

   // build 2,2 matrix
   blk = build2x2(map,  3.0,  9.0,
                        7.0,  1.0);
   A_->setBlock(2,2,blk);

   A_->endBlockFill();

   // build inverse operator
   invA_ = rcp(new Thyra::DefaultBlockedLinearOp<double>());
   invA_->beginBlockFill(3,3);

   // build 0,0 matrix
   blk = build2x2(map,  0.142857142857143,   0.428571428571429,
                        0.285714285714286,  -0.142857142857143);
   invA_->setBlock(0,0,blk);
   invDiag_.push_back(blk);

   // build 0,1 matrix
   blk = build2x2(map, -0.454545454545455,   0.266233766233766,
                       -0.045454545454545,  -0.444805194805195);
   invA_->setBlock(0,1,blk);

   // build 0,2 matrix
   blk = build2x2(map,  0.303571428571429,  -0.271103896103896, 
                        0.069642857142857,   0.090746753246753);
   invA_->setBlock(0,2,blk);

   // build 1,1 matrix
   blk = build2x2(map,  0.090909090909091,  -0.181818181818182,
                        0.045454545454545,   0.159090909090909);
   invA_->setBlock(1,1,blk);
   invDiag_.push_back(blk);

   // build 1,2 matrix
   blk = build2x2(map, -0.050000000000000,  0.086363636363636,
                       -0.045833333333333, -0.019318181818182);
   invA_->setBlock(1,2,blk);

   // build 2,2 matrix
   blk = build2x2(map, -0.016666666666667,   0.150000000000000, 
                        0.116666666666667,  -0.050000000000000);
   invA_->setBlock(2,2,blk);
   invDiag_.push_back(blk);

   invA_->endBlockFill();
}

int tBlockUpperTriInverseOp::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tUpperTriInverseOp";

   status = test_apply(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"apply\" ... PASSED","   \"apply\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_alphabeta(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"alphabeta\" ... PASSED","   \"alphabeta\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      Teko_TEST_MSG(failstrm,0,"tUpperTriInverseOp...PASSED","tUpperTriInverseOp...FAILED");
   }
   else {// Normal Operatoring Procedures (NOP)
      Teko_TEST_MSG(failstrm,0,"...PASSED","tUpperTriInverseOp...FAILED");
   }

   return failcount;
}

bool tBlockUpperTriInverseOp::test_alphabeta(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;
   double diff;

   BlockedLinearOp U = getUpperTriBlocks(A_);
   LinearOp invTri = createBlockUpperTriInverseOp(U,invDiag_);

   RCP<Thyra::VectorBase<double> > src = Thyra::createMember(invA_->domain()); 
   RCP<Thyra::VectorBase<double> > dste = Thyra::createMember(invA_->range()); 
   
   Thyra::randomize<double>(-10,10,src.ptr());
   Thyra::randomize<double>(-10,10,dste.ptr());

   RCP<Thyra::MultiVectorBase<double> > dstn = dste->clone_v();

   diff = Teko::Test::Difference(dste,dstn);
   TEST_ASSERT(diff<=0.0,
          std::endl << "   tBlockUpperTriInverseOp::test_apply " << toString(status)
                    << ": exact copy failed (abserr=" << diff << " <= " << 0.0 << ")" );

   MultiVector dste_mv = dste;
   MultiVector dstn_mv = dstn;
   
   applyOp(invA_,src,dste_mv,3.2,-1.9);
   applyOp(invTri,src,dstn_mv,3.2,-1.9);

   diff = Teko::Test::Difference(dste,dstn)/Thyra::norm_2(*dste);
   TEST_ASSERT(diff<=tolerance_,
          std::endl << "   tBlockUpperTriInverseOp::test_apply " << toString(status)
                    << ": alpha/beta apply operation failed (relerr=" << diff << " <= " << tolerance_ << ")" );

   applyOp(invA_,src,dste_mv);
   applyOp(invTri,src,dstn_mv);

   diff = Teko::Test::Difference(dste,dstn)/Thyra::norm_2(*dste);
   TEST_ASSERT(diff<=tolerance_,
          std::endl << "   tBlockUpperTriInverseOp::test_apply " << toString(status)
                    << ": apply operation (relerr=" << diff << " <= " << tolerance_ << ")" );

   return allPassed;
}

bool tBlockUpperTriInverseOp::test_apply(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > U = getUpperTriBlocks(A_);
   RCP<const Thyra::LinearOpBase<double> > invTri = createBlockUpperTriInverseOp(U,invDiag_);

   Thyra::LinearOpTester<double> tester;
   tester.show_all_tests(true);
   std::stringstream ss;
   Teuchos::FancyOStream fos(Teuchos::rcpFromRef(ss),"      |||");
   const bool result = tester.compare( *invA_, *invTri, Teuchos::ptrFromRef(fos) );
   TEST_ASSERT(result,
          std::endl << "   tBlockUpperTriInverseOp::test_apply "
                    << ": Comparing implicitly generated operator to exact operator");
     if(not result || verbosity>=10) 
        os << ss.str(); 

   return allPassed;
}

} // end Test
} // end Teko
